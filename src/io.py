import os
import yaml


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

    if not os.path.exists(main_dir_fp):
        os.makedirs(main_dir_fp)

    ref_dir_fp = f'{main_dir_fp}/reference_sequences'
    if not os.path.exists(ref_dir_fp):
        os.makedirs(ref_dir_fp)

    foc_dir_fp = f'{main_dir_fp}/focal_sequences'
    if not os.path.exists(foc_dir_fp):
        os.makedirs(foc_dir_fp)
