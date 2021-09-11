from src.Sem3Py.Numerical import solutions as sln
import numpy as np

a = np.array([[1, 1, 1], [2, 3, 7], [1, 3, -2]])
b = np.array([3, 0, 17])

obj = sln()

x = obj.gauss(a, b)
print(x)

