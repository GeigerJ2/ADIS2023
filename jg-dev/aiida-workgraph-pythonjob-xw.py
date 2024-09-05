# # Aiida

import subprocess

# try:
#     subprocess.check_output(["verdi", "profile", "setup", "core.sqlite_dos", "-n", "--profile-name", "test", "--email", "no@email.com"])
# except:
#     pass

from pathlib import Path
from ase.build import bulk

from aiida import orm, engine, load_profile
from aiida.common.exceptions import NotExistent

load_profile('workgraph-dev')

try:
    localhost = orm.load_computer('localhost')
except NotExistent:
    localhost = orm.Computer(
        label='localhost',
        hostname='localhost',
        transport_type='core.local',
        scheduler_type='core.direct',
        workdir=Path('workdir').absolute().as_posix()
    ).store()
    localhost.configure()

try:
    pw_code = orm.load_code('pw@localhost')
except NotExistent:
    pw_code = orm.InstalledCode(
        label='pw',
        computer=localhost,
        filepath_executable='pw.x',
        default_calc_job_plugin='aiida_qe_basic.pw',
        prepend_text='export OMP_NUM_THREADS=1'
    ).store()

# ! rabbitmq-server -detached

# ! sleep 5

# results = engine.run(builder)

# results

# results['properties'].get_dict()

# ## Equation of State curve - WorkGraph from Xing

from aiida_workgraph import task, WorkGraph
from ase import Atoms
from ase.build import bulk
from ase.io.espresso import Namelist
from ase.calculators.espresso import Espresso, EspressoProfile
# ? from ase_quantumespresso.espresso import Espresso, EspressoProfile -> Renamed this. Why `ase_quantumespresso`?

# ### Creating the structure


atoms = bulk('Al', a=4.05, cubic=True)
structure = orm.StructureData(ase=atoms)

# ### Building blocks taken from `workgraph-collections`

# #### `generate_scaled_atoms` and `fit_eos`


@task(
    outputs=[
        {"name": "scaled_atoms", "identifier": "Namespace"},
        {"name": "volumes"},
    ]
)
def generate_scaled_atoms(atoms: Atoms, scales: list) -> dict:
    """Scale the structure by the given scales."""
    volumes = {}
    scaled_atoms = {}
    for i in range(len(scales)):
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scales[i], scale_atoms=True)
        scaled_atoms[f"s_{i}"] = atoms1
        volumes[f"s_{i}"] = atoms1.get_volume()
    return {"scaled_atoms": scaled_atoms, "volumes": volumes}


@task()
def fit_eos(volumes: dict, scf_results: dict) -> dict:
    """Fit the EOS of the data."""
    from ase.eos import EquationOfState
    from ase.units import kJ

    volumes_list = []
    energies = []
    for key, data in scf_results.items():
        energy = data["energy"]
        energies.append(energy)
        volumes_list.append(volumes[key])
    #
    eos = EquationOfState(volumes_list, energies)
    v0, e0, B = eos.fit()
    # convert B to GPa
    B = B / kJ * 1.0e24
    eos = {"energy unit": "eV", "v0": v0, "e0": e0, "B": B}
    return eos


# #### PwCalculator


@task(
    outputs=[
        {"name": "atoms"},
        {"name": "results"},
    ]
)
def pw_calculator(
    atoms: Atoms,
    pseudopotentials: dict,
    kpts: list = None,
    kspacing: float = None,
    command: str = "pw.x",
    input_data: dict = None,
    pseudo_dir: str = "./pseudopotentials",
    calculation: str = None,
) -> dict:
    """Run a Quantum Espresso calculation on the given atoms object."""

    input_data = {} if input_data is None else input_data

    from ase.io.espresso import Namelist
    from ase.calculators.espresso import Espresso, EspressoProfile

    profile = EspressoProfile(command=command, pseudo_dir=pseudo_dir)

    input_data = Namelist(input_data)
    input_data.to_nested(binary="pw")
    # set the calculation type
    if calculation:
        input_data.setdefault("CONTROL", {})
        input_data["CONTROL"]["calculation"] = calculation

    # Set the output directory
    input_data.setdefault("CONTROL", {})
    input_data["CONTROL"]["outdir"] = "out"

    calc = Espresso(
        profile=profile,
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=kpts,
        kspacing=kspacing,
    )

    atoms.calc = calc

    atoms.get_potential_energy()
    results = atoms.calc.results
    new_atoms = results.pop("atoms")
    # we only update the position and cell of the atoms object
    atoms.positions = new_atoms.positions
    atoms.cell = new_atoms.cell
    # Set atoms.calc to None to avoid pickling error
    atoms.calc = None
    return {"atoms": atoms, "results": results}


# #### `all_scf` function

@task.graph_builder(outputs=[{"name": "scf_results", "from": "context.results"}])
def all_scf(scaled_atoms, scf_inputs):
    """Run the scf calculation for each atoms."""

    wg = WorkGraph()
    for key, atoms in scaled_atoms.items():
        scf = wg.tasks.new(
            "PythonJob", function=pw_calculator, name=f"scf_{key}", atoms=atoms
        )
        scf.set(scf_inputs)
        # save the output parameters to the context
        scf.set_context({"results": f"results.{key}"})
    return wg

# #### EOS WorkGraph


@task.graph_builder(outputs=[{"name": "result", "from": "fit_eos.result"}])
def eos_workgraph(
    atoms: Atoms = None,
    command: str = "pw.x",
    computer: str = "localhost",
    scales: list = None,
    pseudopotentials: dict = None,
    pseudo_dir: str = None,
    kpts: list = None,
    input_data: dict = None,
    metadata: dict = None,
    run_relax: bool = True,
):
    """Workgraph for EOS calculation.
    1. Get the scaled atoms.
    2. Run the SCF calculation for each scaled atoms.
    3. Fit the EOS.
    """
    from copy import deepcopy

    input_data = input_data or {}

    wg = WorkGraph("EOS")
    # -------- relax -----------
    if run_relax:
        relax_task = wg.tasks.new(
            "PythonJob",
            function=pw_calculator,
            name="relax",
            atoms=atoms,
            metadata=metadata,
            computer=computer,
        )
        relax_input_data = deepcopy(input_data)
        relax_input_data.setdefault("CONTROL", {})
        relax_input_data["CONTROL"]["calculation"] = "vc-relax"
        relax_task.set(
            {
                "command": command,
                "input_data": relax_input_data,
                "kpts": kpts,
                "pseudopotentials": pseudopotentials,
                "pseudo_dir": pseudo_dir,
            }
        )
        atoms = relax_task.outputs["atoms"]
    # -------- scale_atoms -----------
    scale_atoms_task = wg.tasks.new(
        "PythonJob",
        function=generate_scaled_atoms,
        name="scale_atoms",
        atoms=atoms,
        scales=scales,
        computer=computer,
        metadata=metadata,
    )
    # -------- all_scf -----------
    all_scf1 = wg.tasks.new(
        all_scf,
        name="all_scf",
        scaled_atoms=scale_atoms_task.outputs["scaled_atoms"],
        scf_inputs={
            "command": command,
            "input_data": input_data,
            "kpts": kpts,
            "pseudopotentials": pseudopotentials,
            "pseudo_dir": pseudo_dir,
            "metadata": metadata,
            "computer": computer,
        },
    )
    # -------- fit_eos -----------
    wg.tasks.new(
        "PythonJob",
        function=fit_eos,
        name="fit_eos",
        volumes=scale_atoms_task.outputs["volumes"],
        scf_results=all_scf1.outputs["scf_results"],
        computer=computer,
        metadata=metadata,
    )
    return wg


# #### Actually submit the `WorkGraph`


# ? metadata dictionary has to be AiiDA data type

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
pseudo_dir = "/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/espresso/pseudo" # -> Must be string, why not Path?
input_data = {
    "system": {
        "occupations": "smearing",
        "degauss": 0.02,
        "smearing": "cold",
    },
}
# metadata = orm.Dict({"options": {"prepend_text": "export OMP_NUM_THREADS=1"}})

wg = eos_workgraph(
    # atoms=atoms,
    atoms=atoms,
    computer="localhost",
    scales=[0.95, 0.98, 1.0, 1.02, 1.05],
    command="mpirun -np 2 pw.x",
    pseudopotentials=pseudopotentials,
    pseudo_dir=pseudo_dir,
    input_data=input_data,
    kpts=(4, 4, 4),
    # metadata=metadata,
)
# ------------------------- Submit the calculation -------------------
# wg.run()
wg.submit(wait=True, timeout=200)




