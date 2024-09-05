#%%
from aiida import load_profile
from aiida import orm
from aiida_workgraph import task, WorkGraph
load_profile()

@task.calcfunction(outputs=[{"name": "structures"}, {"name": "volumes"}])
def scale_structure(structure: orm.StructureData, scales: list):
    """Scale the structure by the given scales."""
    atoms = structure.get_ase()
    volumes = {}
    structures = {}
    for i in range(len(scales)):
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scales[i], scale_atoms=True)
        structure = orm.StructureData(ase=atoms1)
        structures[f"s_{i}"] = structure
        volumes[f"s_{i}"] = structure.get_cell_volume()
    return {"structures": structures, "volumes": orm.Dict(volumes)}


@task.calcfunction()
# because this is a calcfunction, and the input scf_outputs are dynamic, we need use **scf_outputs.
def fit_eos(volumes: dict = None, **scf_outputs):
    """Fit the EOS of the data."""
    from ase.eos import EquationOfState
    from ase.units import kJ

    volumes_list = []
    energies = []
    for key, data in scf_outputs.items():
        unit = data.dict.energy_units
        energy = data.dict.energy
        if unit == "a.u.":  # convert to eV
            energy = energy * 27.21138602
        energies.append(energy)
        volumes_list.append(volumes.get_dict()[key])
    #
    eos = EquationOfState(volumes_list, energies)
    v0, e0, B = eos.fit()
    # convert B to GPa
    B = B / kJ * 1.0e24
    eos = orm.Dict({"energy unit": "eV", "v0": v0, "e0": e0, "B": B})
    return eos


# Output result from context to the output socket
@task.graph_builder(outputs=[{"name": "result", "from": "context.result"}])
def all_scf(structures, scf_inputs):
    """Run the scf calculation for each structure."""
    from aiida_workgraph import WorkGraph
    from aiida_quantumespresso.calculations.pw import PwCalculation

    wg = WorkGraph()
    for key, structure in structures.items():
        # scf = wg.tasks.new(PwCalculation, name=f"scf_{key}", structure=structure)
        # scf.set(scf_inputs)
        # save the output parameters to the context
        scf.set_context({"output_parameters": f"result.{key}"})
    return wg


@task.graph_builder(outputs=[{"name": "result", "from": "fit_eos.result"}])
def eos_workgraph(
    structure: orm.StructureData = None,
    code: orm.Code = None,
    scales: list = None,
    parameters: dict = None,
    kpoints: orm.KpointsData = None,
    pseudos: dict = None,
    metadata: dict = None,
):
    """Workgraph for EOS calculation.
    1. Get the scaled structures.
    2. Run the SCF calculation for each scaled structure.
    3. Fit the EOS.
    """
    wg = WorkGraph("EOS")
    scale_structure1 = wg.tasks.new(
        scale_structure, name="scale_structure", structure=structure, scales=scales
    )
    all_scf1 = wg.tasks.new(
        all_scf,
        name="all_scf",
        structures=scale_structure1.outputs["structures"],
        scf_inputs={
            "code": code,
            "parameters": orm.Dict(parameters),
            "kpoints": kpoints,
            "pseudos": pseudos,
            "metadata": metadata,
        },
    )
    wg.tasks.new(
        fit_eos,
        name="fit_eos",
        volumes=scale_structure1.outputs["volumes"],
        scf_outputs=all_scf1.outputs["result"],
    )
    return wg