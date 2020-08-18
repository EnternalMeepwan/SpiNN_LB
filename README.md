# SpiNN_LB

## TL;DR
SpiNN_LB is a Master thesis project. SpiNN_LB is a prototype of a basic D2Q9 lattice Boltzmann method implemented on the world-leading neural-inspired hardware archtechture, the SpiNNaker platform.

## Project information

**Student name:** Yuan Feng

**Student number:** s1909558

**EPCC supervisor name(s):** 
    * Dr. Kevin Stratford
**External supervisor(s):**
    * Dr. Alan Stokes (University of Manchester)


## Development Environment/ Dependency
The SpiNNaker is a dependency-hungry program. The whole development is under the following version of SpiNNaker software and do not guarantee to run at any other version.



| software               | tag                    | git hash  |
| ---------------------- | ---------------------- | --------- |
| spinn_common           | SpiNNaker 5.0.0 v5.0   | e954a29   |
| DataSpecification      | ~                      | dab4cf8   |
| SpiNNFrontEndCommon    |                        | 5dab2afe  |
| SpiNNMachine           |                        | dc6d29e   |
| SpiNNMan               |                        | f6d74bf   |
| PACMAN                 |                        | dae8e04e  |
| SpiNNStorageHandlers   |                        | 33494c6   |
| SpiNNUtils             |                        | 3ecf242   |
| sPyNNaker              |                        | 5c35967f0 |
| spalloc                |                        | acdec43   |
| sPyNNaker8             |                        | 27b45f5b  |
| spalloc_server         |                        | 9f16dfe   |
| SpiNNakerGraphFrontEnd |                        | 21d1885   |
| spinnaker_tools        | SpiNNaker 5.0.0 v3.2.5 | 12e29c5   |

If you have any trouble reflog to the git version and just want a check of the results, just use the executables in the `/SpiNN_LB/executable_backup`

## Compile and Run

Before compile, you need a SpiNNaker development installed, if not look here: https://spinnakermanchester.github.io/ .


### Compile:
```bash
cd SpiNN_LB
make
```

### Run 
And there would be two generated files: `lattice_cell.aplx` and `lattice_cell.dict`

Then if you have SpiNNaker board (if you do have a board you probably know how to run), run the python code in the `/SpiNN_LB/pythonCode`, especially`lattice_partitioned.py`. However, this project is not developed with a physical board on the hand. Therefore, the `lattice_partitioned.py` do not guarantee to work.

Most of you may not have a board. Take look at `spinn-20.cs.man.ac.uk/` you can connect to a SpiNNaker machine there in a Jupyter Notebook.
Then upload the `lattice_cell.aplx` and `lattice_cell.dict` and create a Jupyter Notebook with SpyNNaker kernel. Then run the code in `/SpiNN_LB/pythonCode/LB_method.ipynb`. 




## Run a different simulation
These code implement a 128 * 128 lattice Boltzmann method by default in Minion and Brown's paper (https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-96-4037).

If you want to run in a different scale you need to:
    1. Change the kinematic viscosity (nu) and relaxation time (tau) in C according to your physical settings
    2. Change the running time (in timestep) in Python according to your physical settings
    3. Change the simulation scale in Python according to your physical settings
    4. Probably adjust the maximum delay and time_scale_factor in Python and make it work.
    
## Visualization
In the `/SpiNN_LB/pythonCode/LB_method.ipynb`, we implement a contour visualizer of vorticity of the simulation. Use it.


## File Structure
```
.
├── build
├── C_implementation // A LBM implementation in C
│   └── source
│       ├── generate_pic.py // generate a contour graph for visualiztion
│       ├── Makefile
│       └──pure_c_lbm.c
├── evaluation_code
│   └── lbm_test.py
├── LICENSE
├── Makefile // Makefile compile the C 
├── pythonCode
│   ├── global_settings.py
│   ├── lattice_basic_cell.py
│   ├── lattice_edge.py
│   ├── lattice_partitioned.py
│   └── LB_method.ipynb // Jupyter Notebook that run the LBM on SpiNNaker
├── README.md
└── src
    └── lattice_cell.c // A basic lattice (as in lattice Boltzmann) on SpiNNaker in C
```


