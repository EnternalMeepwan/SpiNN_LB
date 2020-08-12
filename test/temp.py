import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math






data_spinnaker = pd.read_csv("data/Vorticity_128*128.csv").dropna(axis = 1).iloc[1:-1,1:-1].to_numpy()
data_C_double = pd.read_csv("C_development/source/Vorticity_double.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()
data_C_float = pd.read_csv("C_development/source/Vorticity_float.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()


#L_1_norm_spinn_double = abs(data_C_double - data_spinnaker).sum()
#L_1_norm_C_float_double = abs(data_C_double - data_C_float).sum()

#L_1_norm_spinn_vs_c_float = abs()


def calL1Num(data1, data2):
    return abs(data1-data2).sum()
    
def calL2Num(data1, data2):
    return math.sqrt(data1-data2)



