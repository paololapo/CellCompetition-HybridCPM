from __future__ import division
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *

from datetime import datetime
import sys
import numpy as np
import random
from tempfile import TemporaryFile
import os

#from PlayerPython import *
#from cc3d import PlayerPython
#from PySteppables import *
#from PySteppablesExamples import MitosisSteppableBase
            
global relaxtime
global adderlist
growthratewt=6.3 #(default: 6.3)
growthratescrb=3.4 #(default: 3.4)
stiffness_kd=0.4 #(default: 0.4)
#stiffness_kd = {{stiffness_kd}}
stiffness_wt=2 #(default: 2.0)

growthrate3=4.8
stiffness_3=1
p_apo_3_coeff=1
volume_3_coeff=1

CI = {{CI}} # Sensitivity to Contact Inhibition (related to k, default: CI=0.1)
adderlist=[]
relaxtime=200
          
          
# Silly steppable to test a simple pipeline and debug
class SillySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """
    def finish(self):
        """
        Called after the last MCS to wrap up the simulation
        """
    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """         
  
  
      
#Initiate Target Volume, with some noise and Lambda
class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
    def start(self):
        for cell in self.cellList:
            # Randomly reposition cells within a predefined range
            x_shift=random.randint(0, 1200)-int(round(cell.xCOM));
            y_shift=random.randint(0, 1600)-int(round(cell.yCOM));
            self.moveCell(cell,[x_shift,y_shift,0])
            
            # If the cell type is 2, reset it to type 1
            if cell.type==2:  #cell.type actually wants the ID
                cell.type=1
                
        # Randomly select 20 cells
        cellarr=list(range(110))
        random.shuffle(cellarr)
        kdcells=cellarr[0:20]
        # Assign cell type 2 to the selected cells
        for cell in self.cellList:
            if cell.id in kdcells:
                cell.type=2             
        
        # Initialize plot window for visualizing initial volume distribution
        volList = [] # List to store initial volumes
        #self.pW=self.addNewPlotWindow(_title='initial seeding volumes' ,_xAxisTitle='initial volume',_yAxisTitle='N', _xScaleType='linear',_yScaleType='linear')
        #self.pW.addHistogramPlot(_plotName='initialVol',_color='green',_alpha=100)# _alpha is transparency 0 is transparent, 255 is opaque        

        # Assign initial target volume and Lambda based on cell type
        for cell in self.cellList:
            id = cell.id

            # Randomly assign target volume using a normal distribution
            # cell.targetVolume=random.randint(minvol-100,maxvol+100)    
            # cell.targetVolume=random.normalvariate(2100, 300)
            cell.targetVolume=random.normalvariate(1800, 500) # Mean: 1800, StdDev: 500
            round(cell.targetVolume)
            
            # Assign Lambda values based on cell type
            if cell.type == 1:  # For cell type 1
                cell.lambdaVolume = stiffness_wt
                add1 = [id, 1800, cell.type]
            elif cell.type == 2:  # For cell type 2
                cell.lambdaVolume = stiffness_kd
                add1 = [id, 1800, cell.type]
            elif cell.type == 3:
                cell.targetVolume = round(random.normalvariate(volume_3_coeff*1800, 500))
                cell.lambdaVolume = stiffness_3
                add1 = [id, volume_3_coeff*1800, cell.type]
                
            # Store cell information in a list for further use
            adderlist.append(add1)
            volList.append(cell.targetVolume)  # Append the volume to the volume list
  
        # Optionally save or plot the volume data
        # self.pW.addHistogram(plot_name='InitialVol', value_array=volList, number_of_bins=100)     
        # self.pW.savePlotAsData('initialVol.txt')    
   
   
   
# Class to simulate cell growth, adjusting the target volume over time based on cell type
class GrowthSteppableLinear(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # Create a plot window to track the average target volume added and total cell volume over time
        #self.pW = self.addNewPlotWindow(
        #    _title='AvgVtAdded',  # Title of the plot
        #    _xAxisTitle='Time (Frames)',  # X-axis label
        #    _yAxisTitle='AvgVtAdded',  # Y-axis label
        #    _xScaleType='linear',  # Linear scaling for the X-axis
        #    _yScaleType='linear'  # Linear scaling for the Y-axis
        #)
        
        # Add a plot to track the average target volume added, displayed as green dots
        #self.pW.addPlot('Avg Vt Added', _style='Dots', _color='green', _size=5)
        # Add a plot to track the average total volume, displayed as red dots
        #self.pW.addPlot('Avg V', _style='Dots', _color='red', _size=5)
        return
        
    def step(self, mcs):
        # Initialize counters for the number of cells, total volume added, and total cell volume
        cellcount = 0
        tvadded = 0  # Total target volume added
        totV = 0  # Total cell volume
        
        # Only execute growth logic after a relaxation time
        if mcs > relaxtime:
            for cell in self.cellList:
                cellcount += 1  # Increment the cell count
                totV += cell.volume  # Accumulate the cell's current volume
                
                # Logic for growing cells of different types
                if cell.type == 1:  # Growth for specific cell types (wild type)
                    cell.targetVolume += round(random.normalvariate(growthratewt, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))
                    tvadded += round(random.normalvariate(growthratewt, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))
                elif cell.type == 2:  # Growth for other specific cell types (scrambled)
                    cell.targetVolume += round(random.normalvariate(growthratescrb, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))
                    tvadded += round(random.normalvariate(growthratescrb, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))
                elif cell.type == 3:  
                    cell.targetVolume += round(random.normalvariate(growthrate3, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))
                    tvadded += round(random.normalvariate(growthrate3, 2.5)) * \
                        np.exp(-(CI / mcs) * ((cell.volume - cell.targetVolume) ** 2))   
                else:
                    # Cells of other types do not grow
                    cell.targetVolume = 0

        # Calculate averages if there are any cells; otherwise, set to zero
        if cellcount == 0:
            avgV = 0
            avgtvadded = 0
        else:
            avgV = totV / cellcount  # Average volume
            avgtvadded = tvadded / cellcount  # Average target volume added

        # Compute the time (adjusted for relaxation time)
        time = (mcs - relaxtime) / float(10)
        
        # Add data points to the plot for visualization and save to file
        #self.pW.addDataPoint('Avg Vt Added', time, avgtvadded)
        #self.pW.savePlotAsData('avg_vt_added.txt')
        #self.pW.addDataPoint('Avg V', time, avgV)
        #self.pW.savePlotAsData('avg_v.txt')

        # Optional: Growth based on a chemical concentration field (commented by default)
        # Uncomment this block to make growth a function of chemical concentration
        # field = CompuCell.getConcentrationField(self.simulator, "PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt = CompuCell.Point3D()
        # for cell in self.cellList:
        #     pt.x = int(cell.xCOM)  # Center of mass X coordinate
        #     pt.y = int(cell.yCOM)  # Center of mass Y coordinate
        #     pt.z = int(cell.zCOM)  # Center of mass Z coordinate
        #     concentrationAtCOM = field.get(pt)  # Get concentration at cell's center
        #     cell.targetVolume += 0.01 * concentrationAtCOM  # Adjust target volume based on concentration

              

# Class handling cell division (mitosis) during the simulation
class MitosisSteppableAdder(MitosisSteppableBase):  
    def __init__(self, frequency=10):
        # Initializes the MitosisSteppableAdder with a specified execution frequency
        MitosisSteppableBase.__init__(self, frequency)
        
    def start(self):
        """
        This function initializes global variables to track cell division information:
        - `parentlist`: Tracks parent and child relationships.
        - `Volumeanalysis`: Logs volume changes during division.
        """
        global parentlist 
        global parentid  
        global childid 
        global Volumeanalysis

        # Initialize variables for tracking parent-child relationships and volumes
        parentlist = []
        parentid = 0
        childid = 0
        Volumeanalysis = []

    def step(self, mcs):
        """
        The main step function:
        - Identifies cells eligible for division based on volume changes.
        - Triggers the division of those cells.
        - Logs division-related data such as parent-child IDs and volumes.
        """
        cells_to_divide = []

        # Check for eligible cells for division only after the relaxation time
        if mcs > relaxtime:
            for cell in self.cellList:
                for x in adderlist:
                    # Conditions for division:
                    # - Cell volume exceeds a threshold based on initial volume and simulation time.
                    # - Cell type is not 9 or 10 (reserved for other uses).
                    if cell.type == 3: f=volume_3_coeff 
                    else: f=1
                    if x[0] == cell.id and cell.volume - x[1] > (f*1800 - (mcs / 60)) and (cell.type not in [9, 10, 11]):
                        cells_to_divide.append(cell)  
                        # Log volume information for analysis
                        Volumeanalysis.append([x[0], (mcs - relaxtime) / float(10), x[1], cell.volume, cell.type])
        
        # Perform division for eligible cells
        for cell in cells_to_divide:
            # Different methods for mitosis; using division along the minor axis here
            self.divideCellAlongMinorAxis(cell)

            # Log parent and child relationships
            parentid = self.parentCell.id
            childid = self.childCell.id
            par1 = [parentid, childid, (mcs - relaxtime) / float(10), cell.type]
            parentlist.append(par1)

            # Log the initial volume of the child cell
            birthvolchild = self.childCell.targetVolume
            addchild = [childid, birthvolchild]
            adderlist.append(addchild)

    def updateAttributes(self):
        """
        Updates attributes of parent and child cells post-division:
        - Reduces the parent's target volume by half.
        - Copies attributes from parent to child.
        - Modifies cell types based on a predefined progression.
        """
        # Halve the parent's target volume post-division
        self.parentCell.targetVolume /= 2  
        self.cloneParent2Child()  # Copy all attributes from parent to child

    def finish(self):
        """
        This function is executed at the end of the simulation.
        Logs such as `parentlist` and `Volumeanalysis` can be saved here for analysis.
        """
        # Uncomment and modify paths to save logs if needed
        # np.savetxt('path_to_save/parentlist-%s' % datetime.now().strftime('%H-%M-%m-%d'), parentlist)
        # np.savetxt('path_to_save/Volumeanalysis-%s' % datetime.now().strftime('%H-%M-%m-%d'), Volumeanalysis)
        return



# Class for simulating random motility of cells in the simulation
class CellMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        # Initializes the steppable with a specified frequency of execution
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        This method is executed once at the beginning of the simulation.
        It can be used to initialize motility parameters for cells.
        """
        # Iterate over all cells in the simulation
        for cell in self.cellList:
            break  # This break stops the loop immediately, so no initialization occurs

            # The following code is commented out and would apply random motility forces to cells
            # Ensures the ExternalPotential plugin is loaded for applying forces to cells

            # Assigns random force components along the X and Y axes
            # `lambdaVecX` and `lambdaVecY` define the magnitude of forces
            cell.lambdaVecX = 10.1 * random.uniform(-1.0, 1.0)  # Force along the X-axis
            cell.lambdaVecY = 10.1 * random.uniform(-1.0, 1.0)  # Force along the Y-axis
            # Uncomment the following line to assign a Z-axis force
            # cell.lambdaVecZ = 0.0  # Force along the Z-axis

    def step(self, mcs):
        """
        This method is executed at each time step with the specified frequency.
        It updates the motility parameters of the cells, simulating random movement.
        """
        # Iterate over all cells in the simulation
        for cell in self.cellList:
            # Assign new random forces along the X and Y axes at each step
            cell.lambdaVecX = 10.1 * random.uniform(-1.0, 1.0)  # Random force along the X-axis
            cell.lambdaVecY = 10.1 * random.uniform(-1.0, 1.0)  # Random force along the Y-axis

            # Uncomment the line below to print the forces for debugging or analysis
            # print('cell.lambdaVecX=', cell.lambdaVecX, ' cell.lambdaVecY=', cell.lambdaVecY)




# Density-dependent apoptosis (cell death)
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        # Initialize the steppable with a specified execution frequency
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called once at the beginning of the simulation.
        Initializes global variables to track local density, apoptosis probabilities, and cumulative apoptosis events.
        """
        global LocalDensity  # Stores density information for individual cells
        global P_apo         # Tracks apoptosis probabilities and related data
        global cumulitiveapop  # Tracks cumulative apoptosis events over time
        cumulitiveapop = []
        LocalDensity = []
        P_apo = []

    def step(self, mcs):
        """
        Executes at every time step (defined by frequency).
        Calculates density for each cell, checks for density-dependent apoptosis, 
        and updates relevant properties and tracking variables.
        """
        apopwt = 0  # Counter for wild-type apoptosis events
        apopscrb = 0  # Counter for scrambled cell apoptosis events
        deathcount = 0  # Counter for total apoptosis events
        dens = 0  # Local density of a cell
        totdens = 0  # Accumulated density across neighbors
        cellcount = 0  # Counter for total cells processed

        if mcs > relaxtime:  # Skip apoptosis calculations during the relaxation period
            for cell in self.cellList:  # Iterate over all cells in the simulation
                if cell.type not in [9, 10, 11]:  # Skip already apoptotic cells
                    totdens = 0
                    cellneighbours = []  # IDs of neighboring cells
                    nieghbourvolume = []  # Volumes of neighboring cells
                    dens = (1 / cell.volume)  # Start with the density of the current cell

                    # Calculate contributions from neighbors
                    for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                        if neighbor:
                            cellneighbours.append(neighbor.id)
                            nieghbourvolume.append(neighbor.volume)
                            dens += (1 / neighbor.volume)  # Add neighbor's density contribution

                    # Total density for this cell
                    totdens += dens
                    cellcount += 1  # Increment cell count

                    # Check for density-dependent apoptosis for scrambled cells (type 2)
                    if cell.type == 2:
                        # Sigmoid-based probability function for density-dependent death
                        if 0.00072194 / (1 + np.exp(-509.4 * (dens * 3 - 0.0067))) > random.random():
                            cell.type = 10  # Mark cell as apoptotic
                            cell.targetVolume = 0
                            cell.lambdaVolume = 2  # Soft constraint for cell volume
                            deathcount += 1
                            apopscrb += 1

                    # Check for density-dependent apoptosis for wild-type cells (type 1)
                    if cell.type == 1:
                        # Sigmoid-based probability function for density-dependent death
                        if 0.00033074 / (1 + np.exp(-235.8 * (dens * 3 - 0.0152))) > random.random():
                            cell.type = 9  # Mark cell as apoptotic
                            cell.targetVolume = 0
                            cell.lambdaVolume = 2  # Soft constraint for cell volume
                            deathcount += 1
                            apopwt += 1

                    # Check for density-dependent apoptosis for type 3 cells
                    if cell.type == 3:
                        # Sigmoid-based probability function for density-dependent death
                        if p_apo_3_coeff*0.00033074 / (1 + np.exp(-235.8 * (dens * 3 - 0.0152))) > random.random():
                            cell.type = 11  # Mark cell as apoptotic
                            cell.targetVolume = 0
                            cell.lambdaVolume = 2  # Soft constraint for cell volume
                            deathcount += 1
                            apopwt += 1

                # Record local density and apoptosis probabilities over time
                time = (mcs - relaxtime) / float(10)  # Time in simulation units
                if cell.type == 2:  # Record for scrambled cells
                    LocalDensity.append([cell.id, cell.type, time, totdens])
                elif cell.type == 1:  # Record for wild-type cells
                    LocalDensity.append([cell.id, cell.type, time, totdens])
                P_apo.append([cell.id, cell.type, time, totdens, deathcount])
                cumulitiveapop.append([time, apopscrb, apopwt])  # Track cumulative apoptosis

    def finish(self):
        """
        Executed at the end of the simulation. Saves recorded data if necessary.
        Uncomment the lines below to save data to files.
        """
        # Save local density data
        # np.savetxt('', LocalDensity)

        # Save apoptosis probability data
        # np.savetxt('', P_apo)

        # Save cumulative apoptosis data
        # np.savetxt('', cumulitiveapop)
        return


        
# Apoptosis based on cell perimeter and shared contact with neighboring cells
class DeathSteppablePerimiter(SteppableBasePy):
    def __init__(self, frequency=10):
        # Initialize the steppable with a specified execution frequency
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called once at the beginning of the simulation.
        Initializes global variables to track perimeter-related data and apoptosis events.
        """
        global perimeterarray  # Stores perimeter-related information for cells
        global P_apo           # Tracks apoptosis probabilities and related data
        global cumulitiveapop  # Tracks cumulative apoptosis events over time
        cumulitiveapop = []
        perimeterarray = []
        P_apo = []

        # Uncomment below to create a plot window for apoptosis tracking
        # self.pW = self.addNewPlotWindow(
        #     _title='Cumulative Apoptosis',
        #     _xAxisTitle='Time (Frames)',
        #     _yAxisTitle='Apoptosis',
        #     _xScaleType='linear',
        #     _yScaleType='linear'
        # )
        # self.pW.addPlot('Apoptosis WT', _style='Dots', _color='green', _size=5)
        # self.pW.addPlot('Apoptosis Scrb', _style='Dots', _color='red', _size=5)

    def step(self, mcs):
        """
        Executes at every time step (defined by frequency).
        Calculates shared perimeter percentages and triggers apoptosis based on contact-dependent probabilities.
        """
        apopwt = 0  # Counter for wild-type apoptosis events
        apopscrb = 0  # Counter for scrambled cell apoptosis events
        deathcount = 0  # Counter for total apoptosis events

        if mcs > relaxtime:  # Skip apoptosis calculations during the relaxation period
            for cell in self.cellList:  # Iterate over all cells in the simulation
                # Skip cells of type 3 (no apoptosis for these cells)
                if cell.type == 3: continue
                
                # Determine whether the cell is scrambled or wild-type
                typecell = 2 if cell.type == 2 else 1
                neighbourpixel = []  # Stores types of neighboring pixels
                perimeterpercentage = 0  # Percentage of the perimeter shared with opposing cell types

                if cell.type not in [9, 10, 11]:  # Skip already apoptotic cells
                    # Iterate over all boundary pixels of the cell
                    for pixel in self.getCopyOfCellBoundaryPixels(cell):
                        a = self.point3DToNumpy(pixel)
                        centrecell = self.cellField[a[0], a[1], a[2]]  # Get the cell at the current boundary pixel

                        # Check neighboring pixels in all four cardinal directions
                        for dx, dy in [(-1, 0), (1, 0), (0, 1), (0, -1)]:
                            neighbor_pixel = self.cellField[a[0] + dx, a[1] + dy, a[2]]
                            if hasattr(neighbor_pixel, "id") and centrecell.id != neighbor_pixel.id:
                                neighbourpixel.append(neighbor_pixel.type)
                            elif not hasattr(neighbor_pixel, "id"):
                                neighbourpixel.append(0)  # Boundary touches the medium

                # Calculate shared perimeter percentage with opposing types
                if typecell == 2 and len(neighbourpixel) > 0:  # Scribble cells
                    perimeterpercentage = (
                        (neighbourpixel.count(1)) / len(neighbourpixel)
                    )
                    # Density-dependent apoptosis probability
                    if 3 * (0.000416 * perimeterpercentage) > random.random():
                        cell.type = 10  # Mark cell as apoptotic
                        cell.targetVolume = 0
                        cell.lambdaVolume = 2  # Enforce soft volume constraint
                        deathcount += 1
                        apopscrb += 1

                elif typecell == 1 and len(neighbourpixel) > 0:  # Wild-type cells
                    perimeterpercentage = (
                        (neighbourpixel.count(2)) / len(neighbourpixel)
                    )
                    # Contact-dependent apoptosis probability
                    if 10 * (0.0000416 * perimeterpercentage) > random.random():
                        cell.type = 9  # Mark cell as apoptotic
                        cell.targetVolume = 0
                        cell.lambdaVolume = 2  # Enforce soft volume constraint
                        deathcount += 1
                        apopwt += 1

                # Record data for analysis
                time = (mcs - relaxtime) / float(10)  # Convert simulation steps to time
                P_apo.append([cell.id, cell.type, time, perimeterpercentage, deathcount])
            # Append data for all cells
            perimeterarray.append([cell.id, cell.type, time, perimeterpercentage])
            cumulitiveapop.append([time, apopscrb, apopwt])

    def finish(self):
        """
        Executed at the end of the simulation. Saves recorded data if necessary.
        Uncomment the lines below to save data to files.
        """
        # Save perimeter-related data
        # np.savetxt('', perimeterarray)
        # Save apoptosis probability data
        # np.savetxt('', P_apo)
        # Save cumulative apoptosis data
        # np.savetxt('', cumulitiveapop)
        return
       
        

# Class to track and plot the count of wild-type and scrambled cells over time
class neighbourdata(SteppableBasePy):
    def __init__(self, frequency=10):
        """
        Initializes the steppable with a specified execution frequency.
        """
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called once at the beginning of the simulation.
        Sets up a plot window for tracking cell counts over time.
        """
        # Create a new plot window to display cell count over time
        #self.pW = self.addNewPlotWindow(
        #    _title='CellCount v Time',  # Title of the plot
        #    _xAxisTitle='Time(Frames)',  # X-axis label
        #    _yAxisTitle='Cell Count',   # Y-axis label
        #    _xScaleType='linear',       # Linear scaling for the x-axis
        #    _yScaleType='linear'        # Linear scaling for the y-axis
        #)

        # Add two plots to the window for wild-type and scrambled cells
        #self.pW.addPlot('cellcount WT', _style='Dots', _color='green', _size=5)  # Wild-type cell count
        #self.pW.addPlot('cellcount scrb', _style='Dots', _color='red', _size=5)  # Scrambled cell count

        # Optional: Uncomment the below section to track cell volume instead
        # self.pW = self.addNewPlotWindow(
        #     _title='Cellvolume v Time',
        #     _xAxisTitle='Time(Frames)',
        #     _yAxisTitle='Cell volume',
        #     _xScaleType='linear',
        #     _yScaleType='log'
        # )
        # self.pW.addPlot('cellvol WT', _style='Dots', _color='blue', _size=5)
        return

    def step(self, mcs):
        """
        Executes at every time step (defined by frequency).
        Counts the number of wild-type and scrambled cells and updates the plot.
        """
        cellcount_wt = 0  # Counter for wild-type cells
        cellcount_scrb = 0  # Counter for scrambled cells

        # Iterate through all cells in the simulation
        for cell in self.cellList:
            if cell.type == 1:  # Wild-type cells
                cellcount_wt += 1
            elif cell.type == 2:  # Scrambled cells
                cellcount_scrb += 1

        # Calculate simulation time based on relaxation time
        time = (mcs - relaxtime) / float(10)

        # Update the plot with the current cell counts
        #self.pW.addDataPoint('cellcount WT', time, cellcount_wt)
        #self.pW.addDataPoint('cellcount scrb', time, cellcount_scrb)

        # Save the plot data to files for analysis
        #self.pW.savePlotAsData('Cell_Count.txt')  # Wild-type cell counts
        #self.pW.savePlotAsData('Cell_Count_scrb.txt')  # Scrambled cell counts

        # Optional: Uncomment below to track average cell volume
        # self.pW.addDataPoint("cellvol WT", time, cellvolume / cellcount)

    def finish(self):
        """
        Called at the end of the simulation.
        Currently performs no additional tasks.
        """
        return

 

# Class to track cell positions and states over time
class tracking(SteppableBasePy):    
    def __init__(self, frequency=10, file_name="simulation"):
        """
        Initializes the steppable with a specified execution frequency.
        """
        SteppableBasePy.__init__(self, frequency)   
        self.file_name = file_name
        self.init_time = 0

    def start(self): 
        """
        Called once at the beginning of the simulation.
        Initializes the global variable `trackingfile` to store tracking data.
        """
        global trackingfile  # Global variable to store tracking data for all cells
        trackingfile = []
        self.init_time = datetime.now()
        print("Start time (h-m-s): ", self.init_time.strftime('%H-%M-%S'))

    def step(self, mcs): 
        """
        Executes at each time step (defined by frequency).
        Tracks cell positions, states, and types, and stores the data in `trackingfile`.
        """
        if mcs%100 == 0: print("mcs=", mcs, end="\r")
        for cell in self.cellList:
            # Calculate simulation time adjusted by relaxation time
            time = (mcs - relaxtime) / float(10)
            
            # Record relevant cell data: position, time, ID, type, and state
            ar1 = [cell.xCOM, cell.yCOM, int(time), int(cell.id), int(cell.type)]
            trackingfile.append(ar1)  # Append the data to the tracking list

    def finish(self):
        """
        Called at the end of the simulation.
        Finalizes the data collection process.
        Optionally saves the tracking data to a file (currently commented out).
        """
        time_taken = datetime.now() - self.init_time
        print("Time taken: ", time_taken)
        # Uncomment the line below to save `trackingfile` data to a file
        np.savetxt(self.file_name, trackingfile, delimiter=',', header="xCOM,yCOM,time,cell_id,cell_type", comments='')
        return
 
            


# Class to clean up cells during simulation based on specific criteria
# This cleanup steppable is designed to remove cells either due to density-dependent death 
# (low volume) or because they are apoptotic (cell types 9 and 10).
class cleanup(SteppableBasePy):     
    def __init__(self, frequency=100):
        """
        Initializes the cleanup steppable with a specified execution frequency.
        Also initializes a global list `extrusions` to track removed cells.
        """
        SteppableBasePy.__init__(self, frequency)     
        global extrusions  # Global variable to store information about extruded (removed) cells
        extrusions = []

    def step(self, mcs): 
        """
        Called at each simulation step based on the specified frequency.
        Handles the removal of cells that meet certain criteria.
        """
        # Only start cleanup after the relaxation time has passed
        if mcs > relaxtime:
            for cell in self.cellList:
                # If the cell's volume is below 500, log its removal
                if cell.volume < 500:
                    extrusions.append([cell.id, cell.type, (mcs - relaxtime) / float(10)]) 
                
                # Remove cells that are apoptotic (type 9 or 10) or have a volume below 500
                if cell.type in [9, 10, 11] or cell.volume < 500:
                    self.deleteCell(cell)  # Deletes the cell from the simulation
                    
    def finish(self):
        """
        Called at the end of the simulation.
        Optionally saves the `extrusions` data to a file for further analysis.
        """
        # Uncomment the line below to save `extrusions` data to a file
        # np.savetxt('', extrusions)   
        return
        


# class SimulationStopperSteppable(SteppableBasePy):
#     def __init__(self, frequency=1, stop_mcs=1000):
#         SteppableBasePy.__init__(self, frequency)
#         self.stop_mcs = stop_mcs

#     def step(self, mcs):
#         """
#         Check if the current MCS matches the stopping point and stop the simulation.
#         :param mcs: Current Monte Carlo Step
#         """
#         if mcs >= self.stop_mcs:
#             print(f"Stopping simulation at MCS {mcs}")
#             self.stop_simulation()