
import numpy as np
from numpy import array


class pvalues:
    #create a class of operations on pvalues
    #attribute length, order, order of order, array format
    def __init__(self, p):
        self.p = p
        self.n = self.p.size
        self.o = np.argsort(self.p)[::-1]
        self.ro = np.argsort(self.o)
        self.i = np.arange(1,self.n+1,1)[::-1]
    
    def BH(self):    
        return(np.minimum(1, np.minimum.accumulate(self.n/self.i*self.p[self.o]))[self.ro])
    
    def BY(self):
        self.q = sum(1/np.arange(1,self.n+1,1))
        return(np.minimum(1, np.minimum.accumulate((self.q *self.n)/self.i*self.p[self.o]))[self.ro])
    
        