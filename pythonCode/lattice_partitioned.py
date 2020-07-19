# Copyright (c) 2017-2019 The University of Manchester
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import math
import os
# from pacman.model.graphs.machine import MachineEdge

import spinnaker_graph_front_end as front_end

# 12000
runtime = 120000
time_step = 10000
MAX_X_SIZE_OF_FABRIC = 128
MAX_Y_SIZE_OF_FABRIC = 128
n_chips = (MAX_X_SIZE_OF_FABRIC * MAX_Y_SIZE_OF_FABRIC) // 10

ex = [0, 1, 0, -1, 0, 1, -1, -1, 1]
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1]



def initVelocity(x_pos, y_pos):
    """
    Init the velocity u, v = (u_x, u_y)
  
    ydim = xdim = N = 128
    K = 30 / N
    delta = 0.05
  
    u = tanh( K (y - 0.25 * ydim) ) for y <= 0.5 * ydim 
    u = U_0 tanh( K (0.75 * ydim - y) ) for y >  0.5 * ydim
    v = delta sin( 2pi *x / xdim )
    """
    U_0 = 0.01
    K = 30.0
    delta = 0.05
    x_temp = 1.0 * (x_pos) / MAX_X_SIZE_OF_FABRIC
    y_temp = 1.0 * (y_pos) / MAX_Y_SIZE_OF_FABRIC
    if y_temp <= 0.5:
        u_x = U_0 * math.tanh(K * (y_temp - 0.25))
    else:
        u_x = U_0 * math.tanh(K * (0.75 - y_temp))
    u_y = U_0 * delta * math.sin(2 * math.pi * (x_temp + 0.25))
    return u_x, u_y


# set up the front end and ask for the detected machines dimensions
front_end.setup(
    n_chips_required=n_chips, model_binary_folder=os.path.dirname(os.path.abspath("__file__")), machine_time_step=time_step,time_scale_factor=10)

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
