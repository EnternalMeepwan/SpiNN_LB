# 
#   This is a basic implementation of a lattice Boltzmann method on SpiNNaker.
#   This project is currently being released under the GPL 3.0 license. Use it Free as in
#   Freedom. 
#  
#   The up-to-date information of SpiNNaker Project can be found here:
#   https://spinnakermanchester.github.io/
#  

import math
import os
from pacman.model.graphs.machine import MachineEdge
import spinnaker_graph_front_end as front_end
from lattice_basic_cell import LatticeBasicCell
from lattice_edge import LatticeEdge
import global_settings




n_chips = (MAX_X_SIZE_OF_FABRIC * MAX_Y_SIZE_OF_FABRIC) // 10

# set up the front end and ask for the detected machines dimensions
front_end.setup(
    n_chips_required=n_chips, model_binary_folder=os.path.dirname(os.path.abspath("__file__")), machine_time_step=time_step,time_scale_factor=time_scale_factor)

# figure out if machine can handle simulation
cores = front_end.get_number_of_available_cores_on_machine()

# chech if there is enough cores
if cores <= (MAX_X_SIZE_OF_FABRIC * MAX_Y_SIZE_OF_FABRIC):
    raise KeyError("Don't have enough cores to run simulation")

# contain the vertices for the connection aspect
vertices = [
    [None for _ in range(MAX_X_SIZE_OF_FABRIC)]
    for _ in range(MAX_Y_SIZE_OF_FABRIC)]

# build vertices
for x in range(0, MAX_X_SIZE_OF_FABRIC):
    for y in range(0, MAX_Y_SIZE_OF_FABRIC):
        u_x, u_y = initVelocity(x, y)
        vert = LatticeBasicCell(
            "cell{}".format((x * MAX_X_SIZE_OF_FABRIC) + y),
            x, y, u_x, u_y)
        vertices[x][y] = vert
        front_end.add_machine_vertex_instance(vert)

# build edges
for x in range(0, MAX_X_SIZE_OF_FABRIC):
    for y in range(0, MAX_Y_SIZE_OF_FABRIC):
        #   diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
        positions = [
            (x, (y + 1) % MAX_Y_SIZE_OF_FABRIC, "E"),
            ((x + 1) % MAX_X_SIZE_OF_FABRIC,
             (y + 1) % MAX_Y_SIZE_OF_FABRIC, "SE"),
            ((x + 1) % MAX_X_SIZE_OF_FABRIC, y, "S"),
            ((x + 1) % MAX_X_SIZE_OF_FABRIC,
             (y - 1) % MAX_Y_SIZE_OF_FABRIC, "SW"),
            (x, (y - 1) % MAX_Y_SIZE_OF_FABRIC, "W"),
            ((x - 1) % MAX_X_SIZE_OF_FABRIC,
             (y - 1) % MAX_Y_SIZE_OF_FABRIC, "NW"),
            ((x - 1) % MAX_X_SIZE_OF_FABRIC, y, "N"),
            ((x - 1) % MAX_X_SIZE_OF_FABRIC,
             (y + 1) % MAX_Y_SIZE_OF_FABRIC, "NE")]

         # build edges for each direction for this vertex
        for (dest_x, dest_y, compass) in positions:
            front_end.add_machine_edge_instance(LatticeEdge( vertices[x][y], vertices[dest_x][dest_y],compass, "edge between {} and {}".format(vertices[x][y], vertices[dest_x][dest_y])), LatticeBasicCell.PARTITION_ID)
            vertices[x][y].set_direction_vertex(direction=compass, vertex=vertices[dest_x][dest_y])

# run the simulation
front_end.run(runtime)

# get recorded data
recorded_data = dict()

# if not front_end.use_virtual_machine():
buffer_manager = front_end.buffer_manager()

# get the data per vertex
for x in range(0, MAX_X_SIZE_OF_FABRIC):
    for y in range(0, MAX_Y_SIZE_OF_FABRIC):
        recorded_data[x, y] = vertices[x][y].get_data(
            front_end.buffer_manager(),
            front_end.placements().get_placement_of_vertex(
                vertices[x][y]))
# clear the machine
front_end.stop()
