import sys
import os
import msprime
import tskit
import demes
import numpy as np

from src.io import write_fasta_sequences, load_ancestral_sequence


class CoalescentSimulation:
    '''
    Class representing a coalescent simulation run. This can include:
        - Genotype matrices as generated by msprime.
        - Ancestral sequence.
        - Simulation parameters.
    '''

    def __init__(self, config_d, demography=None):
        '''
        Create the simulation parameters, such as demography, focal populations, number of individuals, etc.

        At the moment only configuration file loading is supported, with the exception of being able
        to specify an msprime `Demography` object in addition to the configuration file if no
        demes-compatible file is specified inside the configuration file.

        Args:
            config_d (str): Configuration dictionary containing coalescent simulation parameters.
            demography (Demography, optional): msprime Demography object if none is specified in
                the configuration dictionary. Defaults to None.

        Raises:
            ValueError: On invalid configuration parameters.
        '''
        if 'demography' not in config_d and not demography:
            raise ValueError('No demography specified')
        if 'demography' in config_d and demography:
            raise ValueError('Demography object provided even though configuration already specifies it!')

        # Required parameters
        self.dem_model      = msprime.Demography.from_demes(demes.load(config_d['demography'])) or demography
        self.foc_pops       = config_d.get('focal_populations', [])
        self.foc_pops_sizes = config_d.get('focal_population_sizes', [])
        self.ref_pops       = config_d.get('reference_populations', [])
        self.ref_pops_sizes = config_d.get('reference_population_sizes', [])

        # General simulation parameters
        self.anc_fp  = config_d.get('ancestral_sequence')
        self.seq_len = config_d.get('sequence_length', 10_000)
        self.mu      = config_d.get('mutation_rate', 1.5e-8)
        self.rho_map = config_d.get('recombination_map')
        self.rho     = config_d.get('recombination_rate', None if self.rho_map else 1.5e-8)
        self.ploidy  = config_d.get('ploidy', 2)
        self.out_dir = config_d.get('output_directory', '.')

        # Error handling
        if self.rho and self.rho_map:
            raise ValueError('"recombination_rate" and "recombination_map" parameters are both set!')
        if self.ploidy < 1:
            raise ValueError('"ploidy" must be larger than 0!')
        if self.seq_len < 1:
            raise ValueError('"sequence_length" must be larger than 0!')
        if self.mu < 0:
            raise ValueError('"mutation_rate" can\'t be negative!')
        if self.rho < 0:
            raise ValueError('"recombination_rate" can\'t be negative!')
        if len(self.ref_pops_sizes) != len(self.ref_pops):
            raise ValueError('"reference_populations" and "reference_population_sizes" should be equal in length!')
        if len(self.foc_pops_sizes) != len(self.foc_pops):
            raise ValueError('"focal_populations" and "focal_population_sizes" should be equal in length!')

        # Handle problem with tskit implementation of single-node tree sequences
        if self.ploidy == 1:
            for i in range(len(self.ref_pops_sizes)):
                if self.ref_pops_sizes[i] == 1:
                    print(f'WARNING: Adjusting size of population "{self.ref_pops[i]}" to 2' +
                          ' due to tskit issues with single-node tree sequences.')
                    self.reF_pops_sizes[i] = 2
            for i in range(len(self.foc_pops_sizes)):
                if self.foc_pops_sizes[i] == 1:
                    print(f'WARNING: Adjusting size of population "{self.foc_pops[i]}" to 2' +
                          ' due to tskit issues with single-node tree sequences.')
                    self.foc_pops_sizes[i] = 2

        # Load ancestral sequence if provided
        self.ancestral_sequence = None
        if self.anc_fp:
            self.ancestral_sequence = load_ancestral_sequence(self.anc_fp, self.seq_len)
            self.seq_len = len(self.ancestral_sequence)
            config_d['sequence_length'] = self.seq_len

        # Assume no contamination unless specified
        self.con_pop = config_d.get('contamination_population', None)

        # Assume sampling at present unless specified
        self.foc_pops_times = config_d.get('focal_population_times', [0] * len(self.foc_pops))
        self.ref_pops_times = config_d.get('reference_population_times', [0] * len(self.ref_pops))

        # Set to None so we can build them only if needed
        self.trees = None


    def run(self):
        '''
        Run the coalescent simulation, simulates coalescent trees and mutations.
        '''
        # Build sample dictionaries for sim_ancestry()
        sample_sets = []
        for pop, size, time in zip(self.foc_pops, self.foc_pops_sizes, self.foc_pops_times):
            sample_sets.append(msprime.SampleSet(size, population=pop, time=time, ploidy=self.ploidy))
        for pop, size, time in zip(self.ref_pops, self.ref_pops_sizes, self.ref_pops_times):
            sample_sets.append(msprime.SampleSet(size, population=pop, time=time, ploidy=self.ploidy))
        if self.con_pop:
            sample_sets.append(msprime.SampleSet(1, population=self.con_pop, time=0, ploidy=1))

        # Run coalescent simulation
        self.trees = msprime.sim_ancestry(sample_sets, sequence_length=self.seq_len,
                                          demography=self.dem_model,
                                          recombination_rate=self.rho or self.rho_map)
        self.trees = msprime.sim_mutations(self.trees, rate=self.mu)

        # Map population string IDs to integer IDs
        self.pop_id_d = {}
        for p in self.trees.populations():
            self.pop_id_d[p.metadata['name']] = p.id


    def write_sequences(self):
        '''
        Write fasta files for all of the simulated individuals.

        Raises:
            RuntimeError: If no coalescent simulation has taken place.

        Returns:
            Tuple of: (list of tuples (metadata, filepath) of the written focal sequences,
                       filepath for contamination sequence)
        '''
        if not self.trees:
            raise RuntimeError('Need to run coalescent simulation before writing VCFs!')

        if not self.ancestral_sequence:
            self.ancestral_sequence = tskit.random_nucleotides(self.seq_len)

        # Write samples in reference populations
        for pop in self.ref_pops:
            curr_samples = self.trees.samples(population=self.pop_id_d[pop])

            # Don't use contamination individual if it belongs to this population
            if self.con_pop == pop:
                curr_samples = np.delete(curr_samples, -1)

            curr_seqs = self.trees.alignments(reference_sequence=self.ancestral_sequence,
                                              samples=curr_samples)

            write_fasta_sequences(os.path.join(self.out_dir, 'reference_sequences'),
                                  pop, curr_seqs, self.ploidy)

        # Write contamination sequence if applicable
        cont_fp = None
        if self.con_pop:
            con_sample = self.trees.samples(population=self.pop_id_d[self.con_pop])[-1]
            con_seq    = self.trees.alignments(reference_sequence=self.ancestral_sequence,
                                               samples=[con_sample])

            cont_fp = write_fasta_sequences(os.path.join(self.out_dir, 'miscellaneous'),
                                            'CONTAMINATION_SEQUENCE', con_seq, 1)[0][1]

        # Write ancestral sequence
        write_fasta_sequences(os.path.join(self.out_dir, 'miscellaneous'),
                              'REFERENCE_SEQUENCE', [self.ancestral_sequence], 1)

        # Write samples in focal populations
        ret = []
        for pop in self.foc_pops:
            curr_samples = self.trees.samples(population=self.pop_id_d[pop])

            # Don't use contamination individual if it belongs to this population
            if self.con_pop == pop:
                curr_samples = np.delete(curr_samples, -1)

            curr_seqs = self.trees.alignments(reference_sequence=self.ancestral_sequence,
                                              samples=curr_samples)

            ret.extend(write_fasta_sequences(os.path.join(self.out_dir, 'focal_sequences'),
                                             pop, curr_seqs, self.ploidy))

        return ret, cont_fp
