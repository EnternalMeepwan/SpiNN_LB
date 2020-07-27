import math
import random
random.seed(100)
runtime = 12000
time_step = 1000
time_scale_factor = 5
cores_per_chip = 10

max_offset = 1500
def generate_offset(processor):
    return min(math.ceil((processor - random.random()) / (cores_per_chip+4) * max_offset), max_offset)
