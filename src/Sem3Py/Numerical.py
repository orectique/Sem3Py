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

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        arr = []
        a = lbound
        b = ubound
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
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
    
    def false_position(self, x0: float, x1: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using bisection method.

        :param x0: The lower bound of the range to evaluate.
        :type lbound: int

        :param x1: The upper bound of the range to evaluate.
        :type ubound: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        arr = []
        a = x0
        b = x1
        lst = [a, b]
        i = 2
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a - b) >= delta:
                lst.append(lst[0] - (lst[i - 1] - lst[0])*(self.func(lst[0]))/(self.func(lst[i-1]) - self.func(lst[0])))
                if itn == 0:
                    arr = [itn, i, lst[i-1], lst[i], self.func(lst[i])]
                else:
                    arr = np.vstack([arr, [itn, i, lst[i-1], lst[i], self.func(lst[i])]])
                a = lst[i - 2]
                b = lst[i - 1]
                itn += 1
                i += 1
            return arr
    
    def secant(self, x0: float, x1: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using bisection method.

        :param x0: The lower bound of the range to evaluate.
        :type lbound: int

        :param x1: The upper bound of the range to evaluate.
        :type ubound: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        arr = []
        a = x0
        b = x1
        lst = [a, b]
        i = 2
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a - b) >= delta:
                lst.append(lst[i-1] - (self.func(lst[i - 1])*(lst[i-1] - lst[i-2]))/(self.func(lst[i-1]) - self.func(lst[i-2])))
                if itn == 0:
                    arr = [itn, i, lst[i-1], lst[i], self.func(lst[i])]
                else:
                    arr = np.vstack([arr, [itn, i, lst[i-1], lst[i], self.func(lst[i])]])
                a = lst[i - 2]
                b = lst[i - 1]
                itn += 1
                i += 1
            return arr
        return
    
    def ridder():
        return
    
    def newton_raphson():
        return

