import numpy as np
from sdeint import itoint
import matplotlib.pyplot as plt

from numba import jit, njit

from timeit import timeit

m = np.vectorize(lambda x: max(0, x))


class LotkaVolterra:

    def __init__(self, S, mu, sigma, l, K, T, symm, scaled = True):
        self.S = S
        self.mu = mu
        self.sigma = sigma
        self.l = l
        self.K = K
        self.T = T
        if scaled:
            self.A = np.random.normal(mu / S, sigma / np.sqrt(S), (S, S))
        else:
            self.A = np.random.normal(mu, sigma, (S, S))
        np.fill_diagonal(self.A, 0)

        if symm:
            self.A = (self.A + self.A.T) / 2

    def fG(self, k=0.75):
        def f(x, t):
            y = m(x)
            return y * (y**(k - 1) - (y / self.K)**3 - (self.A @ y)) + self.l

        def G(x, t):
            y = m(x)
            return np.diag(np.sqrt(2 * self.T * y))
        return f, G

    def fG_numba(self, k=0.75):
        A = np.ascontiguousarray(self.A)
        K = self.K
        T = self.T

        @jit("f8[:](f8[:], f8)")
        def f(x, t):
            #x = np.ascontiguousarray(x)
            return x * (x**(k - 1) - x / K - A @ x)

        @jit("f8[:, :](f8[:], f8)")
        def G(x, t):
            #x = np.ascontiguousarray(x)
            return np.diag(np.sqrt(T * x))
        return f, G

    def solve(self, k, x0, max_time=100, dt=1, numba=False):
        if numba:
            f, G = self.fG_numba(k)
        else:
            f, G = self.fG(k)
        return itoint(f, G, y0=x0, tspan=np.linspace(0, max_time, round(max_time / dt)))

    def plot(self, x0, k=0.75, max_time=100, dt=1, numba=False):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))

        np.random.seed(0)
        axes[0].plot(self.solve(1, x0, max_time, dt, numba))
        np.random.seed(0)
        axes[1].plot(self.solve(k, x0, max_time, dt, numba))
        plt.tight_layout()


if __name__ == '__main__':
    lv = LotkaVolterra(S=50, mu=1, sigma=0.1, l=.1, K=1e2, T=0.1, symm=True)
    lv.solve(0.75, x0=np.ones(lv.S), max_time=50, dt=.01)
