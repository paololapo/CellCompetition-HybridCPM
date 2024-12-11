from cc3d import CompuCellSetup
from CompetitionSteppables import *
import os

name = "ScanJ23_"+str({{J23}})
print("Starting simulation: ", name)

path_to_save = "/home/lapo/Desktop/Uni/CellCompetition/TripleCompetition/Files"
path_to_save = "/mnt/CellCompetition/TripleCompetition/Files"
file_name = os.path.join(path_to_save, name)

#CompuCellSetup.register_steppable(steppable=countType1(frequency=1))
#CompuCellSetup.register_steppable(steppable=SillySteppable(frequency=1))


# Registering steppables in CompuCell3D for different simulation behaviors
# A "steppable" is a Python class that defines actions to be taken at each simulation step
# The "frequency" parameter determines how often the steppable runs during the simulation (in terms of time steps)

# Registering the ConstraintInitializerSteppableAdder, which initializes constraints for cells
# It will run every 1 simulation step (frequency=1)
CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppableAdder(frequency=1))

# Registering the GrowthSteppableLinear, which manages cell growth in a linear fashion
CompuCellSetup.register_steppable(steppable=GrowthSteppableLinear(frequency=10))

# Registering the MitosisSteppableAdder, which handles cell division (mitosis)
CompuCellSetup.register_steppable(steppable=MitosisSteppableAdder(frequency=10))

# Registering the CellMotilitySteppable, which defines cell movement during the simulation
CompuCellSetup.register_steppable(steppable=CellMotilitySteppable(frequency=10))

# Registering the DeathSteppable, which manages apoptosis (cell death) in the simulation
#CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=10))

# Registering steppeable for cell death due to biochemical competition
CompuCellSetup.register_steppable(steppable=DeathSteppablePerimiter(frequency=10))

# Registering the neighbourdata steppable, which tracks the number of neighboring cells
CompuCellSetup.register_steppable(steppable=neighbourdata(frequency=10))

# Registering the tracking steppable, which logs information about cell movement and state
CompuCellSetup.register_steppable(steppable=tracking(frequency=10, file_name=file_name))

# Registering the cleanup steppable, which removes cells based on certain conditions (like low volume)
CompuCellSetup.register_steppable(steppable=cleanup(frequency=1000))

# Stop the simulation after a given number of MCS
#CompuCellSetup.register_steppable(steppable=SimulationStopperSteppable(frequency=1, stop_mcs=stop_mcs))

# Start the simulation after all steppables have been registered
CompuCellSetup.run()

