runtime = 12000
time_step = 1000
time_scale_factor = 10
cores_per_chip = 10
_max_offset_factor = 0.2
max_offset = int(time_step * time_scale_factor * _max_offset_factor)

# generate the timer offset for a processor
def generate_offset(processor):
    return int(
        math.ceil(max_offset / cores_per_chip) * processor)
