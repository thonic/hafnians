import numpy as np
from scipy.special import binom
import numdifftools as nd
from collections.abc import Sequence, Callable
from copy import copy


def partial_derivative(func: Callable, var: int, order: int, point: Sequence[float]):
    """

    :param func:
    :param var:
    :param order:
    :param point:
    :return:
    """
    args = point[:]

    def wraps(x):
        args[var] = x
        return func(*args)

    pdf = nd.Derivative(wraps, n=order)
    return pdf(point[var])[0]


def next_n(n):
    max_nonzero = np.max(np.nonzero(n))
    n[max_nonzero] -= 1
    return n


class PreviousIndices(object):
    def __init__(self, n):
        self.n = n

    def __iter__(self):
        self.next_ind = np.zeros(len(self.n))
        return self

    def __next__(self):
        current_index = copy(self.next_ind)
        nz = np.array(np.nonzero(current_index))
        if len(nz.flatten()) == 0:
            self.next_ind[0] = 1
            return current_index

        max_nonzero = np.max(nz)
        if current_index[max_nonzero] < self.n[max_nonzero]:
            self.next_ind[max_nonzero] += 1
            return current_index

        if max_nonzero + 1 < len(self.n):
            self.next_ind[max_nonzero + 1] += 1
            return current_index

        if current_index[-1] == self.n[-1]:
            self.next_ind[-1] += 1
            return current_index

        raise StopIteration()


def recursive_derivatives_exp(f: Callable, x: Sequence[float], n: Sequence[int]):
    """

    :param x:
    :param f:
    :param n:
    :return:
    """
    if np.all(n[1:] == 0) and n[0] == 1:
        return np.array([n]), [np.array([1]), np.array([partial_derivative(f, 0, 1, x)])]

    deeper_n = next_n(n)
    deeper_ns, deeper_bells, deeper_fs = recursive_derivatives_exp(f, x, deeper_n)
    bns = [[binom(deeper_n[i], j) for i, j in enumerate(m)] for m in PreviousIndices()]
    bns = [np.prod(bm) for bm in bns]

    next_y = np.sum(
        bns
        * deeper_bells
        * partial_derivative(
            f,
        )
    )
