# aDNAtor

![@gsusultimate](https://jazpy.github.io/links/CyberMamut.png)

Realistic-ish aDNA simulator.

# Installation

Package can be installed through `pip install adnator`.

`g++` and `OpenMP` need to be installed for read simulations.

# Overview

aDNAtor is a tool for the simulation of both complete sequences (FASTA files) and reads for these
sequences (FASTQ files). aDNAtor's primary use case is the study of damaged DNA; users can configure
parameters such as average read coverage, read fragmentation, misincorporation, genotyping error, etc.

Typical execution is split into two parts:

1. Coalescent simulation
2. Read simulation

In the coalescent simulation, `msprime` is used to simulate genealogies, mutations, recombination events,
and ground-truth nucleotide sequences. These sequences are used as the starting point for the read simulation.

The read simulation takes these ground-truth sequences, and randomly samples from them according
to user-provided parameters, introducing alterations such as genotyping error or deamination events.

aDNAtor's behavior is specified through a configuration file in `.yaml` format. The configuration
options available are detailed below.

## General parameters

`output_directory`: filepath for directory where all sequence and read files will be stored.
Four directories will be created inside `output_directory`:

  1. `focal_reads`: FASTQ files with simulated reads, with one file per individual.
  2. `focal_sequences`: FASTA files with ground-truth sequences, with one file per chromosome in the focal populations.
  3. `miscellaneous`: FASTA files for reference (ancestral) and contamination sequences.
  4. `reference_sequences`: FASTA files with ground-truth sequences, with one file per chromosome in the reference populations.

## Coalescent simulation parameters

`demography` (optional): filepath to a `demes` file specifying the demographic history for a set
of populations. If not present in the configuration file, an msprime `Demography` object needs
to be provided to aDNAtor's `Simulation` object's constructor.

`focal_populations`: list of strings corresponding to population IDs. aDNAtor will simulate both the
ground-truth sequences for these individuals, as well as FASTQ files resulting from read simulation.

`focal_population_sizes`: list of integers detailing how many individuals to simulate for each
population in `focal_populations`.

`focal_population_times` (optional): list of integers detailing how many generations in the past to
sample the individuals in `focal_populations`, defaults to sampling from the present (0 generations in the past).

`reference_populations` (optional): list of strings corresponding to population IDs. aDNAtor will
only simulate ground-truth FASTA sequences for these individuals, without introducing any kind of alterations.

`reference_population_sizes` (optional): list of integers detailing how many individuals to simulate for each
population in `reference_populations`.

`reference_population_times` (optional): list of integers detailing how many generations in the past to
sample the individuals in `reference_populations`, defaults to sampling from the present (0 generations in the past).

`ancestral_sequence` (optional): filepath to a FASTA file. This sequence will be used as the ancestral sequence
for all simulations. If not specified, a random string of nucleotides will be used for the ancestral sequence.

`sequence_length` (optional): length of the sequences to simulate, defaults to 10,000 base pairs.

`mutation_rate` (optional): mutation rate to use for coalescent simulations, defaults to 1.5e-8

`recombination_rate` (optional): recombination rate to use for coalescent simulations, defaults to 1.5e-8

`recombination_map` (optional): filepath to a recombination map in HapMap format. If specified, this recombination
map will be used for coalescent simulations.

`ploidy` (optional): ploidy of simulated individuals, defaults to 2.

## Read simulation parameters

`average_coverage` (optional): average coverage to simulate for FASTQ files, defaults to 5.

`fragmentation_distribution` (optional): filepath to a file detailing a read length distribution. This file is made up
of two columns without a header. The first column is the length of the read, and the second column is the probability
of a read having the corresponding length. Values in the second column should add up to 1.

`fragment_length` (optional): constant read length to simulate if no `fragmentation_distribution` argument is provided,
defaults to 70.

`misincorporation_files` (optional): list of two filepaths, corresponding to `5p_freq_misincorporations.txt` and
`3p_freq_misincorporations.txt` files as generated by `damageprofiler`. If provided, misincorporation will be simulated
for all reads following the specified distributions.

`genotyping_error` (optional): boolean value, used to enable or disable simulation of genotyping error. Defaults to False.

`contamination_population` (optional): string corresponding to a population ID. If provided, an extra chromosome
from this population will be simulated to serve as the source of contaminated reads.

`contamination_proportion` (optional): floating point value between 0 and 1, indicates the proportion of reads
that will be contaminated, defaults to 0.

`contamination_sequence` (optional): filepath to FASTA sequence to use as the source of contaminated reads.

# Example Usage

In order to run a simulation on the included demographic model `utilities/example_demography.yaml`, which specifies
two focal populations and two reference populations, with the following parameters:

1. Sequence length of 100kbp.
2. Sampling 5 individuals from focal population `FOC0`, 10 generations in the past.
3. Sampling 10 individuals from focal population `FOC1`, 50 generations in the past.
4. Sampling 5 individuals from reference population `REF0` in the present.
5. Sampling 10 individuals from reference population `REF1` in the present.
6. Providing the sequence in `utilities/ancestral_sequence.fasta` as the ancestral sequence.
7. With an average coverage of 1X.
8. With a contamination individual from population `REF0`, and a contamination proportion of 2%.
9. Simulating reads to follow the fragmentation distribution in `utilities/example_fragmentation_distribution.txt`.
10. Simulating the misincorporation rates detailed in `utilities/example_5p_misincorporations.txt` and `utilities/example_3p_misincorporations.txt`.
11. Placing all results in `example_data/`.

We would write the following configuration file (provided in `utilities/example_configuration.yaml`):

```yaml
# General simulation parameters
output_directory: './example_data/'

# Coalescent simulation parameters
demography: 'utilities/example_demography.yaml'
sequence_length: 100000
focal_populations: ['FOC0', 'FOC1']
focal_population_sizes: [5, 10]
focal_population_times: [10, 50]
reference_populations: ['REF0', 'REF1']
reference_population_sizes: [5, 10]
ancestral_sequence: 'utilities/ancestral_sequence.fasta'

# Read simulation parameters
average_coverage: 1
contamination_population: 'REF0'
contamination_proportion: 0.02
fragmentation_distribution: 'utilities/example_fragmentation_distribution.txt'
misincorporation_files: ['utilities/example_5p_misincorporations.txt', 'utilities/example_3p_misincorporations.txt']
```

We can then execute the coalescent and read simulations from Python:

```python
from adnator.simulation import Simulation


# Create simulation object with a configuration file
sim = Simulation('utilities/example_config.yaml')
# Run coalescent simulation (creates directories according to configuration file).
sim.run_coalescent_simulation()
# Run read and misincorporation simulation
sim.run_read_simulation()
```
