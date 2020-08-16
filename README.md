# SpiNN_LB

## TL;DR
SpiNN_LB is a Master thesis project. SpiNN_LB is a prototype of a basic D2Q9 lattice Boltzmann method implemented on the world-leading neural-inspired hardware archtechture, the SpiNNaker platform.


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

Or take look at `spinn-20.cs.man.ac.uk/` you can connect to a SpiNNaker machine ther.
Then upload the `lattice_cell.aplx` and `lattice_cell.dict` and create a Jupyter Notebook with SpyNNaker kernel. Then run the code in `/SpiNN_LB/pythonCode/LB_method.ipynb`. 

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


