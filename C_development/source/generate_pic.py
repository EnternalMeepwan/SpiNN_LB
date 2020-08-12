import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

X = np.linspace(0,63,64)
Y = np.linspace(0,63,64)
data = pd.read_csv("Vorticity_double64.csv", header=None).dropna(axis = 1)
def f(X,Y):
    return (data.iloc[X,Y])

fig = plt.figure()
plt.title('64 * 64 contour of Vorticity')
ax = plt.subplot(1,1,1)
plt.contour(X, Y, f(X,Y),63, colors='k', linewidths=(0.7,))
plt.contourf(X, Y, f(X,Y),63, alpha=0.25, colors='black', linewidth=(0.7))
ax.set_aspect('equal')
plt.tight_layout()
plt.show()