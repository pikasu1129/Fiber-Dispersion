from cmath import sin
import sympy as sy
import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0, 0.5 * 1e-9, 100000)
fc = 28 * 1e9
y = np.sin(2 * np.pi * fc * (t + 10 * 1e-12))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(t, y, color = 'Red')
ax.set_ylabel("Power")
ax.set_xlabel("Time[s]")
plt.show()