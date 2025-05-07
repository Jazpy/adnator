import os
import sys
import subprocess
import random
import ctypes
import multiprocessing as mp
import numpy as np

from adnator.io import parse_fragmentation_file, parse_damageprofiler_files, MisincorporationData


class ReadSimulation:
    '''
    Class representing a read simulation. This can include:
        - Read locations.
        - Error matrices
        - Simulation parameters.
    '''

    def __init__(self, config_d, num_procs=mp.cpu_count, seed=None):
        '''
        Create the simulation parameters, such as coverage, fragment length distribution,
        damage parameters, etc.

        At the moment only configuration file loading is supported.

        Args:
            config_d (str): Configuration dictionary containing read simulation parameters.
            num_procs (int, optional): CURRENTLY IGNORED. Number of processes to spawn, will
                parallelize simulations to a per chromosome level. Defaults to number of
                virtual cores in the system.
            seed (int, optional): Random seed to use for read simulation (includes, damage, positions, etc.)
        '''

        # Required parameters
        self.foc_pops       = config_d['focal_populations']
        self.foc_pops_sizes = config_d['focal_population_sizes']

        # General simulation parameters
        self.mis_fps   = config_d.get('misincorporation_files', None)
        self.cont_p    = config_d.get('contamination_proportion', 0)
        self.coverage  = config_d.get('average_coverage', 5)
        self.frags     = config_d.get('fragmentation_distribution')
        self.frag_len  = int(config_d.get('fragment_length', 70))
        self.seq_len   = config_d.get('sequence_length', 10_000)
        self.gen_err   = config_d.get('genotyping_error', False)
        self.f_per_hap = config_d.get('file_per_haplotype', False)
        self.ploidy    = config_d.get('ploidy', 2)
        self.out_dir   = os.path.join(config_d.get('output_directory', '.'), 'focal_reads')
        self.num_procs = num_procs

        if self.coverage <= 0:
            raise ValueError('"coverage" must be larger than 0!')
        if self.cont_p < 0:
            raise ValueError('"contamination_proportion" can\'t be negative!')
        if seed == -1:
            raise ValueError('Please use a different random seed for read simulation!')
        self.seed = seed if seed else -1

        # Might be set later
        self.cont_sequence = None
        self.mis_5 = None
        self.mis_3 = None


    def simulate_reads(self, fasta_fps, cont_fp):
        '''
        Simulate the read lengths and positions for all individuals. If a fragmentation distribution
        is provided, this will be taken as the basis for the read lengths, otherwise a constant length
        will be used. Reads / fragments are dispersed uniformly through the sequence.

        Args:
            fasta_fps (list of (metadata, filepath) tuples): Metadata and filepath for all true sequences.
        '''
        # Get fragmentation distribution information if provided
        if self.frags:
            self.lengths, self.probs = parse_fragmentation_file(self.frags)
            avg_length = np.average(self.lengths, weights=self.probs)
            self.num_reads = round(((self.seq_len / avg_length) / self.ploidy) * self.coverage)
        else:
            self.lengths, self.probs = [self.frag_len], [1.0]
            self.num_reads = round(((self.seq_len / self.frag_len) / self.ploidy) * self.coverage)

        # Get misincorporation information if provided
        # TODO Handle this formatting for CPP in a better way, maybe JSON
        mis_5_pos_lst, mis_5_nuc_lst, mis_5_pro_lst = [], [], []
        mis_3_pos_lst, mis_3_nuc_lst, mis_3_pro_lst = [], [], []
        if self.mis_fps:
            self.mis_5, self.mis_3 = parse_damageprofiler_files(self.mis_fps[0], self.mis_fps[1])
            if self.mis_5:
                for (pos, nuc), mis_probs in self.mis_5.weights.items():
                    mis_5_pos_lst.append(pos)
                    mis_5_nuc_lst.append(nuc)
                    mis_5_pro_lst.extend(mis_probs)
            if self.mis_3:
                for (pos, nuc), mis_probs in self.mis_3.weights.items():
                    mis_3_pos_lst.append(pos)
                    mis_3_nuc_lst.append(nuc)
                    mis_3_pro_lst.extend(mis_probs)

        # Hand off to CPP
        src_dir, _ = os.path.split(__file__)
        so_path = os.path.join(src_dir, 'read_sim.so')

        # Check that .so exists, try creating it if not
        if not os.path.isfile(so_path):
            try:
                subprocess.run(['make'], cwd=src_dir, check=True)
            except:
                print('ERROR: Failed to build shared CPP library. Are Clang and OpenMP installed?')
                sys.exit(-1)

        cpp_library = ctypes.cdll.LoadLibrary(so_path)

        pop_lst, ind_lst, chr_lst, seq_lst, out_lst = [], [], [], [], []
        for (population, individual, chromosome), seq_fp in fasta_fps:
            pop_lst.append(population)
            ind_lst.append(str(individual))
            chr_lst.append(str(chromosome))
            seq_lst.append(seq_fp)
            if self.f_per_hap:
                out_lst.append(os.path.join(self.out_dir, f'{population}_{individual}_{chromosome}.fastq'))
            else:
                out_lst.append(os.path.join(self.out_dir, f'{population}_{individual}.fastq'))

        # General parameters
        pop_arr = (ctypes.c_char_p * len(pop_lst))()
        ind_arr = (ctypes.c_char_p * len(ind_lst))()
        chr_arr = (ctypes.c_char_p * len(chr_lst))()
        seq_arr = (ctypes.c_char_p * len(seq_lst))()
        out_arr = (ctypes.c_char_p * len(out_lst))()
        fra_arr = (ctypes.c_int    * len(self.lengths))()
        pro_arr = (ctypes.c_double * len(self.probs))()

        pop_arr[:] = [i.encode("utf-8") for i in pop_lst]
        ind_arr[:] = [i.encode("utf-8") for i in ind_lst]
        chr_arr[:] = [i.encode("utf-8") for i in chr_lst]
        seq_arr[:] = [i.encode("utf-8") for i in seq_lst]
        out_arr[:] = [i.encode("utf-8") for i in out_lst]
        fra_arr[:] = self.lengths
        pro_arr[:] = self.probs

        # Misincorporation parameters
        mis_5_pos_arr = (ctypes.c_size_t * len(mis_5_pos_lst))()
        mis_5_pro_arr = (ctypes.c_double * len(mis_5_pro_lst))()
        mis_3_pos_arr = (ctypes.c_size_t * len(mis_3_pos_lst))()
        mis_3_pro_arr = (ctypes.c_double * len(mis_3_pro_lst))()

        mis_5_nuc_arr = ''.join(mis_5_nuc_lst).encode('utf-8')
        mis_3_nuc_arr = ''.join(mis_3_nuc_lst).encode('utf-8')

        mis_5_pos_arr[:] = mis_5_pos_lst
        mis_5_pro_arr[:] = mis_5_pro_lst
        mis_3_pos_arr[:] = mis_3_pos_lst
        mis_3_pro_arr[:] = mis_3_pro_lst

        cpp_library.simulate_reads(ctypes.c_size_t(len(fasta_fps)), ctypes.c_size_t(self.seq_len),
                                   ctypes.c_size_t(self.num_reads), ctypes.c_size_t(len(self.lengths)), fra_arr,
                                   pro_arr, ctypes.c_double(self.cont_p), seq_arr, pop_arr, ind_arr, chr_arr,
                                   out_arr, ctypes.c_char_p((cont_fp if cont_fp else '').encode('utf-8')),
                                   ctypes.c_size_t(len(mis_5_pos_lst)), mis_5_pos_arr, mis_5_nuc_arr, mis_5_pro_arr,
                                   ctypes.c_size_t(len(mis_3_pos_lst)), mis_3_pos_arr, mis_3_nuc_arr, mis_3_pro_arr,
                                   ctypes.c_bool(self.gen_err), ctypes.c_int(self.seed), ctypes.c_size_t(self.ploidy))
