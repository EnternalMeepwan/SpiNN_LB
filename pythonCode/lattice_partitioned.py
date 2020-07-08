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
from pacman.model.graphs.machine import MachineEdge
import spinnaker_graph_front_end as front_end

runtime = 500
# machine_time_step = 100
MAX_X_SIZE_OF_FABRIC = 128
MAX_Y_SIZE_OF_FABRIC = 128
n_chips = (MAX_X_SIZE_OF_FABRIC * MAX_Y_SIZE_OF_FABRIC) // 15

ex = [0, 1, 0, -1, 0, 1, -1, -1, 1]
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1]

# Calculate the position of  8 neighbours for a givin Lattice
def generatePositionMap(x_dim, y_dim):
    positionMap = {}
    #   diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
    for i in range(9):
        x_temp = int(x_dim - ex[i]) % MAX_X_SIZE_OF_FABRIC
        y_temp = int(y_dim - ey[i]) % MAX_Y_SIZE_OF_FABRIC
        positionMap[i] = (x_temp, y_temp)
    return positionMap

# Calculate the key for a givin position
def generateKeys(position):
    return position[0] * MAX_X_SIZE_OF_FABRIC + position[1]

#! LB code: initialize the velocity
def initVelocity(x_pos, y_pos):
    U_0 = 0.01
    K = 30.0
    delta = 0.05
    x_temp = 1.0 * (x_pos) / MAX_X_SIZE_OF_FABRIC
    y_temp = 1.0 * (y_pos) / MAX_Y_SIZE_OF_FABRIC
    if y_temp <= 0.5:
        u_x = U_0 * math.tanh(K * (y_temp - 0.25))
    else:
        u_x = U_0 * math.tanh(K * (0.75 - y_temp))
#     print("the u_x of ({}, {}) is {}".format(x_pos, y_pos, u_x))
    u_y = U_0 * delta * math.sin(2.0 * math.pi * (x_temp + 0.25))
#     print("the u_y of ({}, {}) is {}".format(x_pos, y_pos, u_y))    
    return u_x, u_y

# set up the front end and ask for the detected machines dimensions
front_end.setup(
    n_chips_required=1024, model_binary_folder=os.path.dirname(os.path.abspath("__file__")), machine_time_step=10000)

# figure out if machine can handle simulation
cores = front_end.get_number_of_available_cores_on_machine()
print(cores)
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
        keymap = generatePositionMap(x, y)
        vert = LatticeBasicCell(
            "cell{}".format((x * MAX_X_SIZE_OF_FABRIC) + y),
            x, y, generateKeys(keymap[1]), generateKeys(keymap[2]), generateKeys(keymap[3]), generateKeys(keymap[4]),
            generateKeys(keymap[5]), generateKeys(keymap[6]), generateKeys(keymap[7]),
            generateKeys(keymap[8]), u_x, u_y)
        vertices[x][y] = vert
        front_end.add_machine_vertex_instance(vert)

# verify the initial state
output = ""
for x in range(0, MAX_X_SIZE_OF_FABRIC):
    for y in range(0, MAX_Y_SIZE_OF_FABRIC):
        output += "({0}, {1}) ".format(vertices[x][y].x_position, vertices[x][y].y_position)
    output += "\n"
print(output)
print("\n\n")

# build edges
for x in range(0, MAX_X_SIZE_OF_FABRIC):
    for y in range(0, MAX_Y_SIZE_OF_FABRIC):
        keymap = generatePositionMap(x, y)
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

        for (dest_x, dest_y, compass) in positions:
            front_end.add_machine_edge_instance(
                MachineEdge(
                    vertices[x][y], vertices[dest_x][dest_y],
                    label=compass), LatticeBasicCell.PARTITION_ID)

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

# visualise it in text form (bad but no vis this time)
for time in range(0, runtime, 2):
    print("at time {}".format(time))
    output = ""
    for x in range(0, MAX_Y_SIZE_OF_FABRIC):
        for y in range(0, MAX_Y_SIZE_OF_FABRIC):
            output += "{} {}".format(recorded_data[x, y][time], recorded_data[x, y][time+1])
        output += "\n"
    print(output)
    print("\n\n")

# clear the machine
front_end.stop()
