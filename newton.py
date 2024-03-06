''' Newton laws '''

import math
from constants import *

def F(m, a):
	return m * a

def m(F, a):
	return F / a

def a(F, m):
	return F / m

def F_G(m_1, m_2, r):
	return G * m_1 * m_2 / math.pow(r, 2)
