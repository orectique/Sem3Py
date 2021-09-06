from src.Sem3Py.Numerical import roots
import math

func = lambda x : x**2 + x*5 - 8
k = roots(func)
arr = k.bisection(1.0, 2.0, delta = 0.5)
print(arr)