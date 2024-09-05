# # Aiida

import subprocess

try:
    subprocess.check_output(["verdi", "profile", "setup", "core.sqlite_dos", "-n", "--profile", "test", "--email", "no@email.com"])
except:
    pass

from pathlib import Path
from ase.build import bulk

from aiida import orm, engine, load_profile
from aiida.common.exceptions import NotExistent

load_profile()

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

from aiida_qe_basic.pw import PwCalculation

builder = PwCalculation.get_builder()

builder.code = pw_code
builder.structure = orm.StructureData(ase=bulk('Al', a=4.05, cubic=True))
builder.pseudopotentials = orm.Dict({"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"})
builder.parameters = orm.Dict(
    {
        'CONTROL': {
            'calculation': 'scf',
            # 'pseudo_dir': Path('files').absolute().as_posix(),
            'pseudo_dir': Path('espresso/pseudo/').absolute().as_posix(),
            # 'pseudo_dir': '/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/espresso/pseudo
        },
        'SYSTEM': {
            'occupations': 'smearing',
            'smearing': 'cold',
            'degauss': 0.02
        }
    }
)
builder.metadata.options.resources = {
    'num_machines': 1,
    'num_mpiprocs_per_machine': 1
}

# ! rabbitmq-server -detached

# ! sleep 5

results = engine.run(builder)

results

results['properties'].get_dict()

# ## Equation of State curve - basic QE
#
# Running an EOS without all the fancy features in the `aiida-quantumespresso` plugin.

from pathlib import Path

from aiida import orm, engine, load_profile

load_profile()

# ### Importing a structure

from ase.build import bulk

structure = orm.StructureData(ase=bulk('Al', a=4.05, cubic=True))

# ### Relaxing the geometry

resources = {
    'num_machines': 1,
    'num_mpiprocs_per_machine': 1
}

relax_params = {
    'CONTROL': {
        'calculation': 'vc-relax',
        # 'pseudo_dir': Path('files').absolute().as_posix(),
        'pseudo_dir': Path('espresso/pseudo/').absolute().as_posix(),
    },
    'SYSTEM': {
        'occupations': 'smearing',
        'smearing': 'cold',
        'degauss': 0.02
    }
}

from aiida_qe_basic.pw import PwCalculation

builder = PwCalculation.get_builder()

# pseudo_dir = Path(
#     "/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/espresso/pseudo"
# )
# al_pp = orm.UpfData(pseudo_dir / "Al.pbe-n-kjpaw_psl.1.0.0.UPF")

builder.code = orm.load_code("pw@localhost")
atoms = bulk("Al", a=4.05, cubic=True)
structure = orm.StructureData(ase=atoms)
builder.structure = structure

builder.pseudopotentials = orm.Dict({"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"})

# builder.pseudopotentials = orm.Dict({"Al": al_pp})
# builder.pseudopotentials = orm.Dict(
#     {
#         "Al": "/home/geiger_j/aiida_projects/workgraph-dev/git-repos/ADIS2023/espresso/pseudo/Al.pbe-n-kjpaw_psl.1.0.0.UPF"
#     }
# )

# pseudo_family = orm.load_group('SSSP/1.3/PBE/efficiency')
# pseudos = pseudo_family.get_pseudos(structure=structure)
# builder.pseudopotentials = pseudos

builder.parameters = orm.Dict(relax_params)
builder.metadata.options.resources = resources

results = engine.run(builder)
relaxed_structure = results['structure']
relaxed_structure

# ### Calc function to rescale structures
#
# The `calcfunction` below takes an input structure and rescales it to different volumes.

from aiida_qe_basic.pw import PwCalculation

# ? Has to return dict, rather than list

@engine.calcfunction
def rescale_list(structure: orm.StructureData, factor_list: orm.List):

    scaled_structure_dict = {}
    scaled_structure_list = []

    for index, scaling_factor in enumerate(factor_list.get_list()):

        ase_structure = structure.get_ase()

        new_cell = ase_structure.get_cell() * scaling_factor
        ase_structure.set_cell(new_cell, scale_atoms=True)

        scaled_structure_dict[f'structure_{index}'] = orm.StructureData(ase=ase_structure)
        scaled_structure_list.append(orm.StructureData(ase=ase_structure))

    return scaled_structure_dict

# Typically, you'd just run it by calling the function as you would a regular Python function:


# ? Little hack to use default bulk as `relaxed_structure`
relaxed_structure = orm.StructureData(ase=bulk('Al', a=4.05, cubic=True))
rescaled_structures = rescale_list(relaxed_structure, orm.List(list=[0.9, 0.95, 1.0, 1.05, 1.1]))

rescaled_structures

# ## EOS: Work function version

scf_inputs = {
    'CONTROL': {
        'calculation': 'scf',
        # 'pseudo_dir': Path('files').absolute().as_posix(),
    },
    'SYSTEM': {
        'occupations': 'smearing',
        'smearing': 'cold',
        'degauss': 0.02
    }
}

@engine.workfunction
def run_eos_wf(code: orm.Code, structure: orm.StructureData, scale_factors: orm.List):
    """Run an equation of state of a bulk crystal structure for the given element."""

    properties = {}

    for label, rescaled_structure in rescale_list(structure, scale_factors).items():

        builder = PwCalculation.get_builder()
        builder.code = code
        builder.structure = rescaled_structure
        builder.parameters = orm.Dict(scf_inputs)
        builder.pseudopotentials = orm.Dict({"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"})
        builder.metadata.options.resources = resources

        results = engine.run(builder)
        properties[label] = results['properties']

    return properties

results = run_eos_wf(
    code=orm.load_code('pw@localhost'),
    structure=relaxed_structure,
    scale_factors=[0.9, 0.95, 1.0, 1.05, 1.1]
)

results

volumes = []
energies = []

for result in results.values():
    volumes.append(result['volume'])
    energies.append(result['energy'])

import matplotlib.pyplot as plt

plt.plot(volumes, energies)

# ## Work chain version

@engine.calcfunction
def create_eos_dictionary(**kwargs) -> orm.Dict:
    eos = {
        label: (result['volume'], result['energy'])
        for label, result in kwargs.items()
    }
    return orm.Dict(eos)

create_eos_dictionary(**results).get_dict()

class EquationOfState(engine.WorkChain):
    """WorkChain to compute Equation of State using Quantum ESPRESSO."""

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input("code", valid_type=orm.Code)
        spec.input("structure", valid_type=orm.StructureData)
        spec.input("scale_factors", valid_type=orm.List)

        spec.outline(
            cls.run_eos,
            cls.results,
        )
        spec.output("eos_dict", valid_type=orm.Dict)

    def run_eos(self):

        calcjob_dict = {}

        for label, rescaled_structure in rescale_list(self.inputs.structure, self.inputs.scale_factors).items():

            builder = PwCalculation.get_builder()
            builder.code = self.inputs.code
            builder.structure = rescaled_structure
            builder.parameters = orm.Dict(scf_inputs)
            builder.pseudopotentials = orm.Dict({"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"})
            builder.metadata.options.resources = resources

            calcjob_dict[label] = self.submit(builder)

        self.ctx.labels = list(calcjob_dict.keys())

        return calcjob_dict

    def results(self):

        self.report(self.ctx)

        eos_results = {
            label: self.ctx[label].outputs['properties'] for label in self.ctx.labels
        }
        eos_dict = create_eos_dictionary(**eos_results)
        self.out('eos_dict', eos_dict)


engine.run(EquationOfState, code=orm.load_code('pw@localhost'),
           structure=relaxed_structure,
           scale_factors=orm.List([0.9, 0.95, 1.0, 1.05, 1.1]))

# ## Using the `builder`

builder = EquationOfState.get_builder()

builder.structure = relaxed_structure

builder

builder.scale_factors = orm.List([0.9, 0.95, 1.0, 1.05, 1.1])
builder.code = orm.load_code('pw@localhost')

results, node = engine.run_get_node(builder)

results['eos_dict'].get_dict()

eos = node.outputs.eos_dict.get_dict()

eos

plt.plot(
    [v[0] for v in eos.values()],
    [v[1] for v in eos.values()],
)


