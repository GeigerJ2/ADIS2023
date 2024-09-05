# Possible functions to generate the scaled structures

# @task(
#     outputs=[
#         {"name": "scaled_atoms", "identifier": "Namespace"},
#     ]
# )
# # ? Previously was using `Atoms`, but cannot call .get_ase() on the output structure from the relax step when the
# # ? WorkGraph was not submitted
# def generate_structures(atoms: Atoms, scales: list) -> dict:
#     """Scale the structure by the given scales."""
#     scaled_atoms_dict = {}
#     # atoms = structure.get_ase()
#     for i in range(len(scales)):
#         scaled_atoms = atoms.copy()
#         scaled_atoms.set_cell(atoms.cell * scales[i], scale_atoms=True)
#         scaled_atoms_dict[f"structure_{i}"] = scaled_atoms
#     # ? Do we need the nested dictionary here?
#     # return {"scaled_atoms": scaled_atoms_dict}
#     return scaled_atoms_dict

# generate_structures(atoms=atoms, scales=factor_lst)

# @engine.calcfunction
# def generate_structures(structure: orm.StructureData, factor_lst: orm.List):

#     scaled_structure_dict = {}
#     scaled_structure_list = []

#     for index, scaling_factor in enumerate(factor_lst.get_list()):

#         ase_structure = structure.get_ase()

#         new_cell = ase_structure.get_cell() * scaling_factor
#         ase_structure.set_cell(new_cell, scale_atoms=True)

#         scaled_structure_dict[f'structure_{index}'] = orm.StructureData(ase=ase_structure)
#         scaled_structure_list.append(orm.StructureData(ase=ase_structure))

#     return scaled_structure_dict

# import sys
# sys.path.append("/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/adis_tools")

# import importlib.util
# import sys
# spec = importlib.util.spec_from_file_location("adis_tools", "/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/adis_tools")
# foo = importlib.util.module_from_spec(spec)
# sys.modules["adis_tools"] = foo
# spec.loader.exec_module(foo)
# schemas = foo.schemas

# @task.graph_builder(outputs=[{"name": "scf_results", "from": "context.results"}])
# def all_scf(scaled_atoms, scf_inputs):
#     """Run the scf calculation for each atoms."""
#     from aiida_workgraph import WorkGraph
#     from .base import pw_calculator

#     wg = WorkGraph()
#     for key, atoms in scaled_atoms.items():
#         scf = wg.add_task(
#             "PythonJob", function=pw_calculator, name=f"scf_{key}", atoms=atoms
#         )
#         scf.set(scf_inputs)
#         # save the output parameters to the context
#         scf.set_context({"results": f"results.{key}"})
#     return wg

# @task.graph_builder(outputs=[{"name": "result", "from": "context.result"}])
# def all_scf(structures, scf_inputs):
#     """Run the scf calculation for each structure."""
#     from aiida_workgraph import WorkGraph
#     from aiida_quantumespresso.calculations.pw import PwCalculation

#     wg = WorkGraph()
#     for key, structure in structures.items():
#         scf = wg.add_task(PwCalculation, name=f"scf_{key}", structure=structure)
#         scf.set(scf_inputs)
#         # save the output parameters to the context
#         scf.set_context({"output_parameters": f"result.{key}"})
#     return wg

# ? Backup for using PythonJob -> Requires everything via PythonJob
# generate_structures_task = wg.add_task(
#     "PythonJob",
#     function=generate_structures,
#     name="scale_atoms",
#     structure=relax_task.outputs['structure'],
#     factor_lst=factor_lst,
#     computer='localhost',
#     # metadata=metadata,
# )

# print(relax_task.outputs.result.node['result'])
# print(relax_task) # .node.outputs.result.value)
# print(relax_task.outputs["result"]['structure'].get_ase())
# relaxed_structure = relax_task.outputs['result']['structure']
# relaxed_atoms = relaxed_structure.get_ase()

# Run QE using os.subprocess and PythonJob

import subprocess

scf_workdir = Path("/tmp/eos-wg/scf_0/")  # Example path
pwi_input = scf_workdir / "input.pwi"  # Input file path

# Construct the command
command = [
    "pw.x",  # Executable
    "-i", str(pwi_input)  # Arguments
]

# Environment variables
env = {
    **os.environ,  # Copy the existing environment
    "ESPRESSO_PSEUDO": PSEUDO_DIR  # Add/modify the ESPRESSO_PSEUDO variable
}

# Run the command
result = subprocess.run(command, env=env, capture_output=True, text=True, cwd=scf_workdir)

# Check the result
if result.returncode == 0:
    print("Quantum Espresso ran successfully.")
    print(result.stdout)
else:
    print("Error running Quantum Espresso.")
    print(result.stderr)

# run via `launch_shell_job`

def pw_launch_shell_job(input_dict, working_directory):
    # ? How to pass/copy the pseudopotential file? The actual working directory is only available when the shell job is
    # ? being submitted...

    from aiida_shell import launch_shell_job

    # Write input file to working_directory from input_dict
    write_input(input_dict=input_dict, working_directory=working_directory)

    # ? How to create the task only, rather than having to add it to the WG?

    # pw_task = dev_wg.add_task(
    results, node = launch_shell_job(
        # "ShellJob",
        # name="pw",
        command="pw.x",
        arguments="-i {pwi_input}",
        nodes={"pwi_input": orm.SinglefileData(working_directory / "input.pwi")},
        metadata={
            "options": {
                "prepend_text": f"export ESPRESSO_PSEUDO={PSEUDO_DIR}",
                # 'output_filename': 'pwscf.xml',
                # "additional_retrieve": ["pwscf.xml"],
            }
        },
        # parser_outputs=[{"name": "result"}],
        outputs=['pwscf.xml'],
        parser=parse_pw_wg,
    )

    # dev_wg.submit()
    return results, node


results, node = pw_launch_shell_job(
    input_dict=input_dict, working_directory=Path("/tmp/wg-dev")
)

# Original QE input file

# &CONTROL
#    calculation      = 'scf'
#    tstress          = .true.
#    tprnfor          = .true.
#    pseudo_dir       = '/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/espresso/pseudo'
# /
# &SYSTEM
#    occupations      = 'smearing'
#    degauss          = 0.02
#    smearing         = 'cold'
#    ntyp             = 1
#    nat              = 4
#    ibrav            = 0
# /
# &ELECTRONS
# /
# &IONS
# /
# &CELL
# /
# &FCP
# /
# &RISM
# /
# ATOMIC_SPECIES
# Al 26.9815385 Al.pbe-n-kjpaw_psl.1.0.0.UPF

# K_POINTS automatic
# 3 3 3  0 0 0

# CELL_PARAMETERS angstrom
# 4.45056315296600 0.00000000000000 0.00000000000000
# 0.00000000000000 4.45056315296600 0.00000000000000
# 0.00000000000000 0.00000000000000 4.45056315296600

# ATOMIC_POSITIONS angstrom
# Al 0.0000000000 0.0000000000 0.0000000000
# Al 0.0000000000 2.2252815765 2.2252815765
# Al 2.2252815765 0.0000000000 2.2252815765
# Al 2.2252815765 2.2252815765 0.0000000000

# scf_input_dict = MAIN_INPUT_DICT.copy()
# scf_input_dict["calculation"] = "scf"
# scf_input_dict["structure"] = ase_structure
# scf_workdir = main_workdir / f"scf_{ase_structure_key}"
# scf_workdir.mkdir(exist_ok=True)

# write_input(input_dict=scf_input_dict, working_directory=scf_workdir)