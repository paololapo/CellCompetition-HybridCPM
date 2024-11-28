# CompuCell3D: notes
This file contains a list of useful information and notes to start using CompuCell3D:
- [How to install CompuCell3D](#how-to-install-compucell3d)
- [CompuCell3D project organization](#compucell3d-project-organization)

## How to install CompuCell3d
Brief recap of what to do to install CompuCell3D based on the [official documentation](https://compucell3d.org/SrcBin). <br>
**0) Get Miniconda** <br>
*On linux:* install Miniconda from [https://docs.anaconda.com/miniconda/](https://docs.anaconda.com/miniconda/). <br>
*On Windows:* get Linux and install Miniconda. <br>

**1) Setup Miniconda env** <br>
Create a dedicated conda environment:
```
conda create -n cc3d_460_310 python=3.10
```
Activate the environment: <br>
```
conda activate cc3d_460_310 
```

**3) Install mamba** <br>
Install mamba to speed up package resolution:
```
conda install -c conda-forge mamba
```

**4) Install CC3D** <br>
Install CompuCell3D using mamba:
```
mamba install -c main -c conda-forge -c compucell3d compucell3d=4.6.0 
```

**5) Install Twedit++ [optional]** <br>
Install [Twedit](https://github.com/CompuCell3D/cc3d-twedit5) as IDE to write your code:
```
conda install -c compucell3d -c conda-forge cc3d-twedit5
```
Open the IDE and start writing your code:
```
python -m cc3d.twedit5
```

**6) Using CompuCell3D without GUI** <br>
You might want to run a simulation via command line (*batch mode*), without using the graphical interface. You can do it via:
```
python -m cc3d.run_script -i <full path to .cc3d file>
```
Moreover, you might want to make a *direct call* to CompuCell3D from a Python file. In this way, you can better run multiple simulation to compare different settings (e.g. to fine-tune some parameters). This is allowed and explained in the [official documentation](https://compucell3dreferencemanual.readthedocs.io/en/latest/calling_cc3d_directly_from_python.html)


## CompuCell3D project organization
A typical CompuCell3D project is organized in the following files: <br>

**CompuCell3D Project File (`.cc3d`)** <br>
The `.cc3d` file is the **main project file** that ties together all the components of the simulation. It is used to load the Python scripts, XML configuration, and any other simulation resources. When you run a CompuCell3D simulation, this file tells CompuCell3D how to initialize and run the simulation.
- **Purpose**: It links the simulation configuration (from the XML) and the Python scripts (steppables) and serves as the entry point for running the simulation in CompuCell3D.
- **Example**:
    ```xml
    <CompuCell3D>
        <Simulation>
            <Plugin Name="PythonSteppable" File="MySteppable.py" Frequency="10" />
            <Plugin Name="XMLConfig" File="simulation_config.xml" />
        </Simulation>
    </CompuCell3D>
    ```

**CompuCell3D XML File (`.xml`)** <br>
The `.xml` file is used to define the **simulation configuration** and **cell types**. This file contains all the **metadata** for the simulation, such as the number of cells, the model parameters, and the cell behaviors (growth, division, interaction rules). The XML file defines the **simulation grid**, **simulation time steps**, and how the steppables are connected to the overall simulation.
- **Purpose**: This file is the backbone for defining how the simulation environment is set up, the cell types used, and their interactions.
- **Example**:
    ```xml
    <CompuCell3D>
        <Simulation>
            <LengthUnits>microns</LengthUnits>
            <TimeUnits>minutes</TimeUnits>
            <Dimension x="100" y="100" z="100" />
        </Simulation>
        <Plugin Name="Potts" Type="Potts">
            <EnergyTerm Type="CellType" Weight="1.0" />
        </Plugin>
        <CellTypes>
            <CellType Type="1" Name="Type1" />
            <CellType Type="2" Name="Type2" />
        </CellTypes>
    </CompuCell3D>
    ```

**Steppable Python Files (`.py`)** <br>
These files contain the core logic for the **cellular behavior** during each simulation step. A "steppable" refers to a part of the simulation that is executed at each time step. Steppables control processes like cell growth, movement, apoptosis, and division. Each class within a steppable Python file inherits from `SteppableBasePy` and implements the `start()`, `step()`, and `finish()` methods.
- **Example**:
    ```python
    class MySteppable(SteppableBasePy):
        def __init__(self, frequency=1):
            SteppableBasePy.__init__(self, frequency)
        def start(self):
            pass
        def step(self, mcs):
            pass 
        def finish(self):
            pass
    ```

**Run Python file (`.py`)** <br>
The **Run** file contains the script that sets up and launches the entire simulation. It is the entry point for the simulation and registers all the steppables with the CompuCell3D environment. In the run file, you specify the simulation parameters, including which steppables to run, their frequencies, and any initial setup such as cell types or lattice dimensions.
- **Example**:
    ```python
    from Compucell3D import CompuCellSetup
    from my_steppable import MySteppable

    # Registering the steppables
    CompuCellSetup.register_steppable(steppable=MySteppable(frequency=10))

    # Running the simulation
    CompuCellSetup.run()
    ```

