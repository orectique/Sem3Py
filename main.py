from src.Sem3Py.Numerical import roots
import math

k = roots('x**2 + sqrt(x+1)')
arr = k.bisection(-1.0, 2.0, 10)
print(arr)