import os
import multiprocessing as mp
import numpy as np

from src.io import parse_fragmentation_file, parse_damageprofiler_files, MisincorporationData


class ReadSimulation:
    '''
    Class representing a read simulation. This can include:
        - Read locations.
        - Error matrices
        - Simulation parameters.
    '''

    def __init__(self, config_d, num_procs=mp.cpu_count):
        '''
        Create the simulation parameters, such as coverage, fragment length distribution,
        damage parameters, etc.

        At the moment only configuration file loading is supported.

        Args:
            config_d (str): Configuration dictionary containing read simulation parameters.
            num_procs (int, optional): Number of processes to spawn, will parallelize simulations
                to a per chromosome level. Defaults to number of virtual cores in the system.
        '''

        # Required parameters
        self.foc_pops       = config_d['focal_populations']
        self.foc_pops_sizes = config_d['focal_population_sizes']

        # General simulation parameters
        self.mis_fps   = config_d.get('misincorporation_files', None)
        self.cont_p    = config_d.get('contamination_proportion', 0)
        self.coverage  = config_d.get('average_coverage', 5)
        self.frags     = config_d.get('fragmentation_distribution')
        self.seq_len   = config_d.get('sequence_length', 10_000)
        self.ploidy    = config_d.get('ploidy', 2)
        self.out_dir   = os.path.join(config_d.get('output_directory', '.'), 'focal_reads')
        self.num_procs = num_procs

        # Might be set later
        self.cont_sequence = None
        self.mis_5 = None
        self.mis_3 = None


    def simulate_reads_worker(self):
        '''
        Internal function that handles the read position and length calculations. This function
        will be called in parallel for all focal chromosomes

        Returns:
            List of tuples (length, start) describing all simulated reads.
        '''

        ret = []

        frags  = np.random.choice(self.lengths, self.num_reads, p=self.probs)
        starts = np.random.randint(0, self.seq_len - 1, self.num_reads)

        return [(self.seq_len - s if s + f - 1 >= self.seq_len else f, s) for f, s in zip(frags, starts)]


    def write_reads_with_errors_worker(self, sequence_data, reads):
        '''
        Write FASTQ formatted reads to storage, while simulating genotyping error and misincorporation.

        Args:
            sequence_data (tuple of metadata and filepath): Information for the true sequence.
            reads (list of length and start tuples): Information for the read lengths and placement.
        '''

        population, individual, chromosome = sequence_data[0]
        sequence_fp = sequence_data[1]
        out_fp = os.path.join(self.out_dir, f'{population}_{individual}.fastq')

        # Load true sequence into memory
        true_seq = []
        with open(sequence_fp) as in_f:
            for line in in_f:
                if line[0] == '>':
                    continue
                true_seq.extend(line.strip())

        # Random values that determine if read comes from contaminant or not
        cont_errs = np.random.random(len(reads)) < self.cont_p

        # Write each read to FASTQ file
        with open(out_fp, 'a') as out_f:
            for i, ((length, start), cont_err) in enumerate(zip(reads, cont_errs)):
                # Contaminated reads are written directly
                if cont_err:
                    to_write  = f'@SEQ_CONTAMINATED_{population}_{individual}_{chromosome}_{i}\n'
                    to_write += f'{self.cont_sequence[start:start + length]}\n+\n{"j" * length}\n'
                    out_f.write(to_write)
                    continue

                to_write  = f'@SEQ_{population}_{individual}_{chromosome}_{i}\n'
                true_read = true_seq[start:start + length]
                dmg_read  = []

                # Error simulation
                deam_errs = np.random.random(len(true_read))
                geno_errs = np.random.random(len(true_read))
                for i, (nuc, geno_err, deam_err) in enumerate(zip(true_read, geno_errs, deam_errs)):
                    # Deamination error
                    if self.mis_5 and i <= self.mis_5.max_pos:
                        nuc = np.random.choice(list('ATGC'), p=self.mis_5.weights[(i, nuc)])
                    elif self.mis_3 and (length - (i + 1)) <= self.mis_3.max_pos:
                        nuc = np.random.choice(list('ATGC'), p=self.mis_3.weights[(length - (i + 1), nuc)])

                    # Genotyping error
                    if geno_err < 0.00133333:
                        nuc = np.random.choice(list('ATGC'))

                    dmg_read.append(nuc)

                # Write to file
                to_write += f'{"".join(dmg_read)}\n+\n{"j" * length}\n'
                out_f.write(to_write)


    def simulate_reads(self, fasta_fps, cont_fp):
        '''
        Simulate the read lengths and positions for all individuals. If a fragmentation distribution
        is provided, this will be taken as the basis for the read lengths, otherwise a constant length
        of 70bp will be used. Reads / fragments are dispersed uniformly through the sequence.

        Args:
            fasta_fps (list of (metadata, filepath) tuples): Metadata and filepath for all true sequences.
        '''

        if self.frags:
            self.lengths, self.probs = parse_fragmentation_file(self.frags)
            avg_length = np.average(self.lengths, weights=self.probs)
            self.num_reads  = round((self.seq_len / avg_length) / self.ploidy)
        else:
            self.lengths, self.probs = [70], [1.0]
            self.num_reads = round((self.seq_len / 70) / self.ploidy)

        # Create contamination sequence if provided
        if cont_fp:
            self.cont_sequence = []
            with open(cont_fp) as in_f:
                for line in in_f:
                    if line[0] == '>':
                        continue
                    self.cont_sequence.extend(line.strip())

        # Get misincorporation information if provided
        if self.mis_fps:
            self.mis_5, self.mis_3 = parse_damageprofiler_files(self.mis_fps[0], self.mis_fps[1])

        # Simulate reads in parallel
        pool  = mp.Pool(self.num_procs)
        reads = pool.starmap(self.simulate_reads_worker,
                             [() for _ in range(sum(self.foc_pops_sizes) * self.ploidy)])

        pool.starmap(self.write_reads_with_errors_worker,
                     [(seq, read) for seq, read in zip(fasta_fps, reads)])
