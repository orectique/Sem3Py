import numpy as np

class roots:

    """
    Instantiate a function.

    :param func: The function to evaluate.
    :type func: obj
    """
    def __init__(self, func):
        self.func = func

    def bisection(self, lbound: float, ubound: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using bisection method.

        :param lbound: The lower bound of the range to evaluate.
        :type lbound: int

        :param ubound: The upper bound of the range to evaluate.
        :type ubound: int

        :param itn: The number of iterations to evaluate.
        :type itn: int
        
        :param itn: The number of iterations to evaluate.
        :type itn: int

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        arr = []
        a = lbound
        b = ubound
        itn = 0
        while itn <= max_itn and abs(a-b) >= delta:
            x = (a + b)/2
            if itn == 0:
                arr = [a, b, x, self.func(x)]
            else:
                arr = np.vstack([arr, [a, b, x, self.func(x)]])
            if self.func(x) > 0:
                b = x
            elif self.func(x) < 0:
                a = x
            else:
                return arr
            itn += 1

        return arr

        

    def incremental():
        return
    
    def false_position():
        return
    
    def secant():
        return
    
    def riddler():
        return
    
    def newton_rhapson():
        return

