
# find positions of elements in the list
def pozicija(testlist, cond):
    '''
    This function returns a list of indices of elements
    in `testlist` that satisfy the condition `cond`.
    
    Parameters:
        testlist : The list to search.
        cond : The condition to satisfy.
    
    Returns:
        [i for i,x in enumerate(testlist) if cond(x)] : The list of indices of elements that satisfy the condition.
    '''
    return [i for i,x in enumerate(testlist) if cond(x)]

def polynom(c, x):
    '''
    Calculate a polynomial.
    
    Parameters:
        c : The coefficients of the polynomial.
        x : The variable of the polynomial.
    
    Returns:
        npsum(out, axis=0) : The value of the polynomial.
    '''
    from numpy import sum as npsum
    
    out = [k*x**j for j,k in enumerate(c)]
    return npsum(out, axis=0)

def is_number(x):
    return isinstance(x, int) # Test this
    try:
        int(eval(x))
        return True
    except:
        return False

def mode(l): # Why not use np mode?
    
    from collections import Counter
    
    return Counter(l).most_common(1)[0][0]

def moving_average(a, n=3) :
    '''
    Calculate the moving average of a time series.
    
    Parameters:
        a : The time series.
        n : The number of points to average.
    
    Returns:
        ret[n - 1:] / n : The moving average of the time series.
    '''
    from numpy import cumsum
    
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_sum(a, n=3) :
    '''
    Calculate the moving sum of a time series.
    
    Parameters:
        a : The time series.
        n : The number of points to sum.
    
    Returns:
        ret[n - 1:] / n : The moving sum of the time series.
    '''
    from numpy import cumsum
    
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:]

def autocorr(sett, dtrange, nsplits = 1):
    '''
    Autocorrelation function.
    
    Parameters:
        sett : The time series.
        dtrange : The range of time lags.
        nsplits : The number of splits.
    
    Returns:
        out : The autocorrelation function.
    '''
    from numpy import zeros, corrcoef, array, mean, std
    
    if nsplits == 1:
        ret = zeros(len(dtrange))
        for k,i in enumerate(dtrange):
            if i==0:
                ret[k] = 1.
            else:
                ret[k] = corrcoef(sett[:len(sett)-i],sett[i:])[0,1]
        return ret
    else:
        out = []        
        for j in range(nsplits):
            ret = zeros(len(dtrange))
            ss = sett[j*len(sett)//nsplits : (j+1)*len(sett)//nsplits]
            for k,i in enumerate(dtrange):
                if i==0:
                    ret[k] = 1.
                else:
                    ret[k] = corrcoef(ss[:len(ss)-i],ss[i:])[0,1]
            out += [ret]
        out = array(out)
        return ( mean(out,axis=0), std(out,axis=0) )

def OU(theta,mu,sigma,tmax,x0,dt):
    '''
    Simulate an Ornstein-Uhlenbeck process with parameters theta, mu, sigma
    and initial condition x0 for a time tmax with time step dt.
    
    Parameters:
    
    Returns:
    
    '''
    from random import randn
    from numpy import empty
    
    maxindex = int(float(tmax)/dt)
    x = empty(maxindex)
    x[0] = x0
    w  = randn(maxindex)
    a1 = 1.-theta*dt
    a2 = mu*theta*dt
    b  = sigma*dt**.5*w
    for t in range(maxindex-1):
        x[t+1] = a1*x[t] - a2 + b[t]
    
    return x


def order(testlist):
    '''
    Return the order of elements in a list.
    
    Parameters:
        testlist : The list to order.
    
    Returns:
        out : The order of the elements in the list.
    '''
    import numpy as np
    
    tmp = sorted([[i,el] for i,el in enumerate(testlist)], key=lambda xi: xi[1])
    return np.array([el[0] for el in tmp])
    
def tally(mylist):
    ''' 
    Tally elements in a list
    
    Parameters:
        mylist : The list to tally.
    
    Returns:
        out : The tally of the elements in the list.
    '''
    from collections import Counter
    
    return sorted(Counter(mylist).most_common(), key=lambda duple: duple[0])


def multi_map(some_function, iterable, processes=1):
    '''
    Apply a function to an iterable using multiple processes.
    
    Parameters:
        some_function : The function to apply to the iterable.
        iterable : The iterable to apply the function to.
        processes : The number of processes to use.
    
    Returns:
        out : The output of the function applied to the iterable.
    '''
    assert type(processes) == int
    if processes==1:
        out = map(some_function, iterable)
    elif processes>1:
        from multiprocessing import Pool
        pool = Pool(processes)
        out  = pool.map(some_function, iterable)
        pool.close()
        pool.join()
    else:
        print ("invalid number of processes", processes)
        quit()
    return out


from contextlib import contextmanager
@contextmanager
def suppress_stdout():
    ''' Suppress the output of a block of code. '''
    import sys, os
    
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def stochasticMaximize(fun,x0,steps = 10000, temp = 1., step = .1):
    '''
    Stochastic optimization of a function.
    
    Parameters:
        fun : The function to optimize.
        x0 : The initial guess.
        steps : The number of steps to take.
        temp : The temperature of the system.
        step : The step size.
    
    Returns:
        outPars : The optimized parameters.
        fun(outPars) : The value of the function at the optimized parameters.
    '''
    import numpy as np
    from MCMC import MCMC
    
    global exponent, mcmc
    exec(open('MCMCworker_RNApOnly_exclusions.py').read())
    from os.path import expanduser
    nPars = len(x0)
    exponent = fun
    mcmc = MCMC(x0, Nsave=10*nPars, filename=expanduser('~/tmp/mcmc'), step = step, temp = temp, exclude=np.array([],dtype=int))
    mcmc.cycle(steps,adjust=True)
    outPars = np.loadtxt(mcmc.filename+".out", skiprows=steps//10//nPars*9//10)[:,1:-1].mean(axis=0)
    return outPars, fun(outPars)

 
def deep_getsizeof(o, ids):
    """Find the memory footprint of a Python object
 
    This is a recursive function that drills down a Python object graph
    like a dictionary holding nested dictionaries with lists of lists
    and tuples and sets.
 
    The sys.getsizeof function does a shallow size of only. It counts each
    object inside a container as pointer only regardless of how big it
    really is.
 
    :param o: the object
    :param ids:
    :return:
    """
    from collections import Mapping, Container
    from sys import getsizeof

    d = deep_getsizeof
    if id(o) in ids:
        return 0
 
    r = getsizeof(o)
    ids.add(id(o))
 
    if isinstance(o, str): # Removed isinstance(o, unicode) as it is not defined in Python 3
        return r
 
    if isinstance(o, Mapping):
        return r + sum(d(k, ids) + d(v, ids) for k, v in o.iteritems())
 
    if isinstance(o, Container):
        return r + sum(d(x, ids) for x in o)
 
    return r 
