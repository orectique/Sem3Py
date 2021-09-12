import numpy as np

class optimize:
    
    """
    Instantiate a function.

    :param func: The lambda function to evaluate.
    :type func: Lambda function object
    
    :param dx: The derivative of the function passed to 'func'. Optional - Used only in Newton-Raphson method.
    :type dx: Lambda function object
    """
    def __init__(self, func, dx, ddx):
        self.func = func
        self.dx = dx
        self.ddx = ddx
        
    def golden_search(self, lbound: float, ubound: float, max_itn = 100, delta = 0.000001):
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
        
        phi = (1 + 5**0.5)/2
        a = lbound
        b = ubound
        c = b + (a - b)/phi
        d = a + (b - a)/phi
        fc = self.func(c)
        fd = self.func(d)
        itn = 0
        arr = [itn, a, c, d, b]
        
        while abs(b - a) >= delta and itn <= max_itn:
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b + (a - b)/phi
                fc = self.func(c)
            else:
                a = c
                c = d
                fc = fd
                d = (b -  a)/phi
                fd = self.func(d)
            
            arr = np.vstack([arr, [itn, a, c, d, b]])
            itn += 1
        
        return arr
    
    def fibonacci(self, lbound: float, ubound: float, max_itn = 100, delta = 0.000001):
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
        
        a = lbound
        b = ubound
        farr = [1, 1]
        n = 1
        
        
        
        
        return
    
    def newton(self, x0: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using Newton-Raphson method.

        :param x0: The initial value to evaluate.
        :type x0: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        
        x = x0
        itn = 0
        
        while itn <= max_itn and abs(self.dx(x)/self.ddx(x)) >= delta:
            x = x - (self.dx(x)/self.ddx(x))
            
            if itn == 0:
                arr = [itn, x, self.dx(x), self.ddx(x), abs(self.dx(x)/self.ddx(x))]
            
            else: 
                arr = np.vstack([arr, [itn, x, self.dx(x), self.ddx(x), abs(self.dx(x)/self.ddx(x))]])
            itn += 1
            
        return arr
    
    def secant(self, x0: float, x1: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using secant method.

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

        a = x0
        b = x1
        lst = [a, b]
        i = 2
        itn = 0
        
        if self.dx(a)*self.dx(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a - b) >= delta:
                lst.append(lst[i-1] - (self.dx(lst[i - 1])*(lst[i-1] - lst[i-2]))/(self.dx(lst[i-1]) - self.dx(lst[i-2])))
                if itn == 0:
                    arr = [itn, i, lst[i-1], lst[i], self.dx(lst[i])]
                else:
                    arr = np.vstack([arr, [itn, i, lst[i-1], lst[i], self.dx(lst[i])]])
                a = lst[i - 2]
                b = lst[i - 1]
                itn += 1
                i += 1
            return arr
        return
    
    def line_search(self):
        return
    
    def hooke_jeeves(self):
        return
    
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
  
        a = lbound
        b = ubound
        itn = 0
        
        if self.dx(a)*self.dx(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a-b) >= delta:
                x = (a + b)/2
                if itn == 0:
                    arr = [a, b, x, self.dx(x)]
                else:
                    arr = np.vstack([arr, [a, b, x, self.dx(x)]])
                if self.dx(x) > 0:
                    b = x
                elif self.dx(x) < 0:
                    a = x
                else:
                    return arr
                itn += 1

            return arr
        
    def steepest_decent():
        return