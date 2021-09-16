from sublinearLV import *
%matplotlib inline


### May rescaled


lv = LotkaVolterra(S = 20, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = True)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .1)



lv = LotkaVolterra(S = 50, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = True)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .1)

lv = LotkaVolterra(S = 50, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = False)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .01)


lv = LotkaVolterra(S = 50, mu = 0, sigma = 0, l=.1, K = 100, T = 0.01, symm = False)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .01)


### not rescaled


lv = LotkaVolterra(S = 5, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = True, scaled = False)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .01)

lv = LotkaVolterra(S = 50, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = True, scaled = False)

X0 = np.random.random((5, lv.S))
for i in range(5):
    lv.plot(x0 = X0[i], max_time = 50, dt = .01)


# test


lv = LotkaVolterra(S = 100, mu = 1, sigma = 0.1, l=.1, K = 100, T = 0.01, symm = True)
