import sympy as sy
import numpy as np
import math
import matplotlib.pyplot as plt
import random

#8bitで1バイト　→　1GBpsは125MB/s
#1Gbpsの1周期　→　1ns


L = 0.1 #(km)
bit = 25
Lam = 1.55 * 10 ** -6 #(m)
d = 16 #(ps/km*nm)
c = 3 * 10 ** 8
y = 38.25 * 10 ** -3 #(mm)
b2 = (d / (2 * math.pi * c)) * (Lam ** 2) * (10 ** -3)
nm = 3.96 #(電気信号の実行屈折率)
ng = 2.19 #(光波の群屈折率)
y = 38.25 * 10 ** -3 #(mm)

total  = (y / c) * (nm + ng) #(s)
initial = 1000
pitch = 50 * 10 ** -6 #(um)
pitchmm =  pitch * 10 ** 3
dt = pitch * (nm + ng) / c
sumw = (total + dt * initial) / dt

polnumber = 1 + int(sumw) - initial
electrodelength = pitch * polnumber
elecctrodelengthmm = electrodelength * 10 ** 3

#28GHzの信号
#時間tを生成

t = np.arange(0, bit * 25, 24)

digital = []
stepram = []
for n in np.arange(0, bit * 25, 24):
    rm = random.randint(0,1)
    digital.append(rm)

print(digital)

#x = np.arange(0, 24)
y = digital

plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "Bold"

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.step(t, y, color = 'Red')

ax.set_ylabel("Power")
ax.set_xlabel("Time[ps]")
#ax.set_xticklabels([0, 100, 200, 300, 400, 500, 600])
plt.show()

f = 28 * 10 ** 3 #(Hz)
signal = math.cos(2 * math.pi * f)
