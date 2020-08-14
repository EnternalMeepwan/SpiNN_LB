# 
#   This is a basic implementation of a lattice Boltzmann method on SpiNNaker.
#   This project is currently being released under the GPL 3.0 license. Use it Free as in
#   Freedom. 
#  
#   The up-to-date information of SpiNNaker Project can be found here:
#   https://spinnakermanchester.github.io/
#  


import math
import random
random.seed(100)
runtime = 12000
time_step = 1000
time_scale_factor = 5
cores_per_chip = 10

max_offset = 1500

# generate a random timer offser for a givin core.
def generate_offset(processor):
    delay = math.ceil((processor - random.random()) / (cores_per_chip+4) * max_offset)
    if delay > max_offset:
        delay = max_offset - random.random() / (cores_per_chip) * max_offset
    return int(delay)



MAX_X_SIZE_OF_FABRIC = 128
MAX_Y_SIZE_OF_FABRIC = 128


ex = [0, 1, 0, -1, 0, 1, -1, -1, 1]
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1]

# init the velocity for a lattice in position x_pos, y_pos
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