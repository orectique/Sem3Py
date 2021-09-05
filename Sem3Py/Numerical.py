import numpy as np

class roots:

    """
    Instantiate a function.

    :param func: The function to evaluate.
    :type func: int
    """
    def __init__(self, func: str):
        self.func = func

    def bisection(self, lbound: float, ubound: float, itn: int):
        """
        Finding the root of a function using bisection method.

        :param lbound: The lower bound of the range to evaluate.
        :type lbound: int

        :param ubound: The upper bound of the range to evaluate.
        :type ubound: int

        :param itn: The number of iterations to evaluate.
        :type itn: int

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        arr = []
        a = lbound
        b = ubound
        for i in range(itn):
            x = (a + b)/2
            arr = np.vstack([arr, [a, b, x, eval(self.func)]])
            if eval(self.func) > 0:
                b = x
            elif eval(self.func) < 0:
                a = x
            else:
                return arr

        return arr

        


