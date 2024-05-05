import multiprocessing as mp

from src.coalescent_simulation import CoalescentSimulation
from src.read_simulation import ReadSimulation
from src.io import load_yaml, create_directory_structure

class Simulation:
    '''
    Class representing the main simulation object. This object will hold more specific simulations,
    such as the coalescent simulation or the read simulation.
    '''

    def __init__(self, config_fp):
        '''
        Create dictionary with overall simulation parameters, such as demography, populations, read
        paramters, etc.

        Args:
            config_fp (str): Filepath for a YAML configuration file.

        Raises:
            ValueError: On invalid configuration parameters.
        '''

        self.config_d = load_yaml(config_fp)

        # Create with None values so we can still check if they have been set without attribute errors
        self.coalescent_simulation = None
        self.read_simulation = None


    def run_coalescent_simulation(self, demography=None):
        '''
        Run the coalescent simulation. The coalescent simulation will simulate the "true" sequences
        for all individuals in the configuration file and write them to FASTA files. Read simulations
        depend on the coalescent simulation being run first.

        Args:
            demography (msprime.Demography, optional): msprime Demography object if none is specified in
                the configuration dictionary. Defaults to None.

        Raises:
            ValueError: On invalid configuration parameters.
        '''

        self.coalescent_simulation = CoalescentSimulation(self.config_d, demography)
        self.coalescent_simulation.run()

        # Write results of coalescent simulation
        create_directory_structure(self.config_d['output_directory'])
        self.focal_fps = self.coalescent_simulation.write_sequences()


    def run_read_simulation(self, num_procs=mp.cpu_count()):
        '''
        Run the read simulation. A read simulation will be created for every FASTA sequence generated
        in the coalescent simulation. The read simulations will output reads that simulate genotyping error,
        fragmentation, deamination, coverage and depth, etc.

        Args:
            num_procs (int, optional): Number of processes to spawn, will parallelize simulations
                to a per chromosome level. Defaults to number of virtual cores in the system.

        Raises:
            ValueError: If coalescent simulation has not been run first.
        '''

        print(mp.cpu_count())
        if not self.coalescent_simulation or not self.focal_fps:
            raise ValueError('Coalescent simulation must be run before read simulations!')

        self.read_simulation = ReadSimulation(self.config_d, num_procs)
        self.read_simulation.simulate_reads(self.focal_fps)
