import sympy as sy
import numpy as np
import math
import matplotlib.pyplot as plt
import random

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

#28GHzの信号を変調を生成

digital = []

for n in range(0,bit):
    digital[n] = random.randint(0,1)

print(digital[0])


