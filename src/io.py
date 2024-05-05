import os
import yaml
import shutil


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
    if 'reference_populations' not in config:
        raise ValueError('CONFIG ERROR: Need list of reference populations to simulate')
    if 'reference_population_sizes' not in config:
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

    sub_dirs = ['reference_sequences', 'focal_sequences', 'focal_reads']

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
