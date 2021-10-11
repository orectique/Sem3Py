from src.Sem3Py.Numerical import Eigenvalue as egn
from src.Sem3Py.Numerical import Roots
import math

import numpy as np

k = lambda x: x - math.cos(x)
l = lambda x: 1 + math.sin(x)

k = Roots(k, l)
print(k.newton_raphson(1, delta = 0.0000000001))
