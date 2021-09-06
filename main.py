from src.Sem3Py.Numerical import roots
import math

func = lambda x : x**2 + x*5 - 8
k = roots(func)
arr = k.false_position(1.0, 2.0, 10)
print(arr)