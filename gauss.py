

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

data = np.array([[-71.6, 0.2],
[-70.0, 0.3],
[-69.0, 0.1],
[-67.7, 0.5],
[-66.5, 0.2],
[-65.3, 0.6],
[-64.1, 0.4],
[-62.9, 0.7],
[-61.7, 1.1],
[-60.5, 1.5],
[-59.3, 2.2],
[-58.1, 3.3],
[-56.9, 5.1],
[-55.7, 6.9],
[-54.5, 8.5],
[-53.3, 8.8],
[-52.1, 8.4],
[-50.9, 7.7],
[-49.7, 5.6],
[-48.5, 3.8],
[-47.2, 2.2],
[-46.0, 1.7],
[-44.8, 1.3],
[-43.6, 0.8],
[-42.4, 0.6],
[-41.2, 0.8],
[-40.0, 0.4],
[-38.8, 0.9],
[-37.6, 0.1],
[-36.4, 0.8]])

x = data[:,0]
y = data[:,1]


def func(params):
    alpha, mu, sigma = tuple(params)

    return alpha * 1.0 / np.sqrt(2*np.pi) / sigma * np.exp( -(x-mu)**2 / sigma**2 )

def func_residual(params):
    return func(params) - y

fit_result = least_squares(
        func_residual,
        [10.0, -30, 10],
        jac='3-point',
        bounds=[[0.001, -100, 0.1], [10000, -30, 200]],
        loss='linear',
        tr_solver='exact',
        verbose=1
        )

print(fit_result)


def final(x, params):
    alpha, mu, sigma = tuple(params)

    return alpha * 1.0 / np.sqrt(2*np.pi) / sigma * np.exp( -(x-mu)**2 / sigma**2 )

plt.figure()
plt.plot(x, y, 'ro')


xfine = np.linspace(x.min(), x.max(), 1001)
plt.plot(xfine, final(xfine, fit_result.x), 'b-')
print('mu = ', fit_result.x[1] )
plt.axvline( fit_result.x[1] , color='k', linestyle='--')

plt.show()
