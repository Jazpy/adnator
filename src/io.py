import os
import yaml
import shutil

from collections import defaultdict, namedtuple


'''
Misincorporation objects are a tuple of transition probabilities and the maximum position considered by the
original misincorporation file.
'''
MisincorporationData = namedtuple('MisincorporationData', ['max_pos', 'weights'])


def load_yaml(fp):
    '''
    Load contents of file into configuration object while performing basic validation.

    Args:
        fp (str): Filepath for a YAML configuration file.

    Returns:
        Contents of YAML file in Python dictionary.

    Raises:
        ValueError: If required configuration parameters are not set.
    '''
    config = yaml.safe_load(open(fp))

    if 'focal_populations' not in config:
        raise ValueError('CONFIG ERROR: Need list of focal populations to simulate!')
    if 'focal_population_sizes' not in config:
        raise ValueError('CONFIG ERROR: Need number of individuals to simulate for focal populations!')
    if 'reference_population_sizes' not in config and 'reference_populations' in config:
        raise ValueError('CONFIG ERROR: Need number of individuals to simulate for reference populations!')

    return config


def create_directory_structure(main_dir_fp):
    '''
    Creates the directory structure needed to save simulated sequences.

    Args:
        main_dir_fp (str): Filepath for top-level directory where data will be saved.
    '''
    if os.path.exists(main_dir_fp):
        shutil.rmtree(main_dir_fp)
    os.makedirs(main_dir_fp)

    sub_dirs = ['reference_sequences', 'focal_sequences', 'focal_reads', 'miscellaneous']

    for dir in sub_dirs:
        os.makedirs(os.path.join(main_dir_fp, dir))


def write_fasta_sequences(out_dir, population, sequences, ploidy, chunking=60):
    '''
    Writes the FASTA sequences resulting from the coalescent simulation step. This are the "true"
    sequences before the read simulation step.

    Args:
        out_dir (str): Directory to write FASTA files into.
        population (str): Name of the population these sequences belong to.
        sequences (Iterable): Iterator over sequences to write.
        ploidy (int): Number of chromosomes that make up a single individual.
        chunking (int, optional): Number of base pairs to write per-line in FASTA file. Defaults to 60.

    Returns:
        List of tuples (metadata, filepath) of the written sequences.
    '''
    ret = []

    # Write sequences in chunks, assuming contguous sequences belong to the same individual
    for i, seq in enumerate(sequences):
        individual = i // ploidy
        chromosome = i %  ploidy
        out_fp = f'{out_dir}/{population}_{individual}_{chromosome}.fasta'
        ret.append(((population, individual, chromosome), out_fp))
        with open(out_fp, 'w') as out_f:
            out_f.write(f'>{population}_{individual}_{chromosome}\n')
            out_f.write('\n'.join([seq[x:x + chunking] for x in range(0, len(seq), chunking)]))

    return ret


def parse_fragmentation_file(frag_fp):
    '''
    Handles parsing of fragmentation length distribution file, should be in a two column format. The first
    column corresponds to fragment length, the second column corresponds to the proportion of fragments
    with the corresponding length. Second column should add up to 1.

    Args:
        frag_fp (str): Filepath for the fragmentation distribution file.

    Raises:
        ValueError: On invalid file format.

    Returns:
        Tuple of (lengths, probabilities) describing fragmentation length distribution.
    '''
    lengths = []
    probs   = []

    try:
        with open(frag_fp) as in_f:
            for line in in_f:
                toks = line.split()
                lengths.append(int(toks[0]))
                probs.append(float(toks[1]))

        assert len(lengths) > 0
    except:
        raise ValueError('Exception while parsing fragmentation file, is it present and well-formed?')

    # Truncate final probability if excess is small
    if sum(probs) != 1.0:
        if sum(probs) < 1.005 and max(probs) > 0.01:
            probs[probs.index(max(probs))] -= sum(probs) - 1.0
        else:
            raise ValueError('Probabilities in fragmentation file do not add up to 1!')

    return (lengths, probs)


def __parse_damageprofiler_file(fp):
    '''
    Auxiliary function to handle the parsing of a single damageprofiler file.

    Args:
        fp (str): Filepath for the damageprofiler file.

    Returns:
        Misincorporation object encoding probability of misincorporation given position and nucleotide,
        along with information for how much of the strand is considered.
    '''
    max_pos = 0
    weights = defaultdict(list)

    with open(fp) as in_f:
        line = in_f.readline().strip()

        while line[0] == '#' or not line:
            line = in_f.readline().strip()

        # Parse header
        idx_mapper = {}
        toks = line.split()
        for i, tok in enumerate(toks):
            if len(tok) == 3 and tok[1] == '>':
                idx_mapper[tok] = i

        # Parse position lines
        for line in in_f:
            if line[0] == '#' or not line.strip():
                continue

            toks = line.split()
            curr_pos = int(toks[0])
            max_pos  = max(curr_pos, max_pos)

            # Build weights for all possible transitions, starting with A>N
            a_t = float(toks[idx_mapper['A>T']])
            a_g = float(toks[idx_mapper['A>G']])
            a_c = float(toks[idx_mapper['A>C']])
            a_a = 1.0 - sum([a_t, a_g, a_c])
            weights[(curr_pos, 'A')] = [a_a, a_t, a_g, a_c]
            # T>N
            t_a = float(toks[idx_mapper['T>A']])
            t_g = float(toks[idx_mapper['T>G']])
            t_c = float(toks[idx_mapper['T>C']])
            t_t = 1.0 - sum([t_a, t_g, t_c])
            weights[(curr_pos, 'T')] = [t_a, t_t, t_g, t_c]
            # G>N
            g_a = float(toks[idx_mapper['G>A']])
            g_t = float(toks[idx_mapper['G>T']])
            g_c = float(toks[idx_mapper['G>C']])
            g_g = 1.0 - sum([g_a, g_t, g_c])
            weights[(curr_pos, 'G')] = [g_a, g_t, g_g, g_c]
            # C>N
            c_a = float(toks[idx_mapper['C>A']])
            c_t = float(toks[idx_mapper['C>T']])
            c_g = float(toks[idx_mapper['C>G']])
            c_c = 1.0 - sum([c_a, c_t, c_g])
            weights[(curr_pos, 'C')] = [c_a, c_t, c_g, c_c]

    return MisincorporationData(max_pos=max_pos, weights=weights)


def parse_damageprofiler_files(dmg_5_fp, dmg_3_fp):
    '''
    Handles parsing of damageprofiler misincorporation files into a dictionary encoding the probability
    of misincorporation given distance from read start and base in question.

    Args:
        dmg_5_fp (str): Filepath for the damageprofiler misincorporation file from the 5' end.
        dmg_3_fp (str): Filepath for the damageprofiler misincorporation file from the 3' end.

    Raises:
        ValueError: On invalid file format.

    Returns:
        Misincorporation objects encoding probability of misincorporation given position and nucleotide,
        along with information for how much of the strand is considered.
    '''
    try:
        mis_5 = __parse_damageprofiler_file(dmg_5_fp)
        mis_3 = __parse_damageprofiler_file(dmg_3_fp)
    except:
        raise ValueError('Exception while parsing damageprofiler file, is it present and well-formed?')

    return mis_5, mis_3


def load_ancestral_sequence(seq_fp, seq_len):
    '''
    Handles loading of a fasta file into memory to be used as an ancestral sequence for
    coalescent simulations.

    Args:
        seq_fp (str): Filepath for the sequence in FASTA format.
        seq_len (int): Maximum number of base pairs to load.

    Raises:
        ValueError: On invalid file format.

    Returns:
        FASTA sequence as a string.
    '''
    seqs = []
    curr_len = 0
    with open(seq_fp) as in_f:
        next(in_f)
        for line in in_f:
            line = line.strip()
            seqs.append(line)
            curr_len += len(line)

            if curr_len >= seq_len:
                break

    if curr_len < seq_len:
        print(f'WARNING: Provided ancestral sequence is shorter than requested length. New length = {curr_len}.')

    return ''.join(seqs)[:seq_len]
