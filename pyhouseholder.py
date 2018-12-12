#!/usr/bin/env python

import numpy as np
from math import copysign, hypot

a = np.matrix([
    [1, 0, 0, 0],
    [2, 1, 0, 0],
    [3, 2, 1, 0],
    [4, 3, 2, 1]])
q, r = np.linalg.qr(a)
print 'Q'
print q
print 'R'
print r
