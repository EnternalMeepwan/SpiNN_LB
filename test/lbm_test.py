import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math



#C double
c_double_256=pd.read_csv("data/Vorticity_double256.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()
c_double_128=pd.read_csv("data/Vorticity_double128.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()
c_double_64=pd.read_csv("data/Vorticity_double64.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()

#C float
c_float_256=pd.read_csv("data/Vorticity_float256.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()
c_float_128=pd.read_csv("data/Vorticity_float128.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()
c_float_64=pd.read_csv("data/Vorticity_float64.csv", header=None).dropna(axis=1).iloc[1:-1,1:-1].to_numpy()

#spinn
spi_256=pd.read_csv("data/Vorticity_spinn_256.csv").dropna(axis = 1)
spi_128=pd.read_csv("data/Vorticity_spinn_128.csv").dropna(axis = 1).iloc[1:-1,1:-1].to_numpy()
spi_64=pd.read_csv("data/Vorticity_spinn_64.csv").dropna(axis = 1).iloc[1:-1,1:-1].to_numpy()


def calL1Num(data1, data2):
    return abs(data1-data2).sum()
    
def calL2Num(data1, data2):
    return math.sqrt(((data1-data2)**2).sum())

def calLmaxNum(data1, data2):
    return (abs(data1-data2)).max()



spi_256

