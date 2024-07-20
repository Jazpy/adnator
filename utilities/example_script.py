from adnator.simulation import Simulation


# Create simulation object with a configuration file
sim = Simulation('utilities/example_config.yaml')
# Run coalescent simulation (creates directories according to configuration file).
sim.run_coalescent_simulation()
# Run aDNA damage simulation
sim.run_read_simulation()
