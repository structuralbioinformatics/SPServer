from scipy.stats import hypergeom as hyp
import numpy     as np


def accumulative_hypergeometric(k, n, K, N):
    '''
    [k]: SUCCESS IN THE CLUSTER
    [n]: SIZE OF THE CLUSTER
    [K]: SUCCESS IN POPULATION
    [N]: SIZE OF THE POPULATION
    '''
    k, n, K, N = int(k), int(n), int(K), int(N)
    sf = hyp.sf(k, N, K, n)
    if sf < 1:
        return sf + hyp.pmf(k, N, K, n)
    else:
        return 1 - hyp.cdf(k, N, K, n) + hyp.pmf(k, N, K, n)


def frequency(k, n):
    '''
    [k]: SELECTION
    [n]: POPULATION
    '''
    return float(k)/float(n)


def logodds(k, n, K, N):
    '''
    [k]: SUCCESS IN THE CLUSTER
    [n]: SIZE OF THE CLUSTER
    [K]: SUCCESS IN POPULATION
    [N]: SIZE OF THE POPULATION
    '''
    return np.log(frequency(k, n)/frequency(K, N))


def mutual_information(k, n, K, N):
    '''
    [k]: SUCCESS IN THE CLUSTER
    [n]: SIZE OF THE CLUSTER
    [K]: SUCCESS IN POPULATION
    [N]: SIZE OF THE POPULATION
    '''
    return frequency(k, N)*logodds(k, n, K, N)
