#!/usr/bin/python3
''' Script to calculate '''

import newton
from constants import *

F_Sun_Min = newton.F_G(M_Sun, 1, AE + r_Earth)
F_Sun_Max = newton.F_G(M_Sun, 1, AE - r_Earth)
F_Moon_Min = newton.F_G(M_Moon, 1, r_Moon_Earth + r_Earth)
F_Moon_Max = newton.F_G(M_Moon, 1, r_Moon_Earth - r_Earth)
F_Earth = newton.F_G(M_Earth, 1, r_Earth)

A_F_Sun = 1.0 - (F_Earth - F_Sun_Max) / F_Earth
A_F_Moon = 1.0 - (F_Earth - F_Moon_Max) / F_Earth

print(F_Sun_Max - F_Sun_Min)
print(F_Moon_Max - F_Moon_Min)
