import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

def fit_parabola(x1, y1, x2, y2, x3, y3):
    denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    b = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom
    c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
    return a, b, c

xfull = np.linspace(1, 264, 264)

# INFINITY
x1 = np.array([1, 6, 100, 264])
y1 = np.array([1, 0.4, 0.2, 0.2])
tck = scipy.interpolate.splrep(x1, y1)
f = scipy.interpolate.interp1d(x1, y1, kind='linear')

a, b, c = fit_parabola(x1[0], y1[0], x1[1], y1[1], x1[3], y1[3])
y1full = a + b*xfull + c*xfull**2

#y1full = f(xfull)
#y1full = scipy.interpolate.splev(xfull, tck)

# X
x2 = np.array([1, 6, 100, 264])
y2 = np.array([1, 0.6, 0.8, 1.])
tck = scipy.interpolate.splrep(x2, y2)
f = scipy.interpolate.interp1d(x2, y2, kind='linear')
y2full = f(xfull)
#y2full = scipy.interpolate.splev(xfull, tck)

# 4X
x3 = np.array([1, 6, 10, 100, 264])
y3 = np.array([1, 0.5, 0.6, 0.7, 0.8])
tck = scipy.interpolate.splrep(x3, y3)
f = scipy.interpolate.interp1d(x3, y3, kind='linear')
y3full = f(xfull)
#y3full = scipy.interpolate.splev(xfull, tck)

# 8X
x4 = np.array([1, 6, 10, 20, 100, 264])
y4 = np.array([1, 0.4, 0.5, 0.6, 0.65, 0.7])
tck = scipy.interpolate.splrep(x4, y4)
f = scipy.interpolate.interp1d(x4, y4, kind='linear')
y4full = f(xfull)
#y4full = scipy.interpolate.splev(xfull, tck)


fig = plt.figure()
plt.plot(xfull, y1full)
plt.plot(xfull, y2full)
plt.plot(xfull, y3full)
plt.plot(xfull, y4full)
plt.legend([r'$\infty$ Histories', 'X Histories', '4X Histories', '8X Histories'])
plt.show()
