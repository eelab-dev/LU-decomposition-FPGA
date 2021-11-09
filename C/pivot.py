import numpy as np

x = np.array(([0.003, 59.14, 59.17], [5.291, -6.13, 46.78]))
y = np.array(([5.291, -6.13, 46.78], [0.003, 59.14, 59.17]))
x1 = np.copy(x)
y1 = np.copy(y)

x1[1, 0:] = x1[1, 0:] - x[0, 0:] * x[1, 0] / x[0, 0]
y1[1, 0:] = y1[1, 0:] - y[0, 0:] * y[1, 0] / y[0, 0]

print(x1)
xd = np.zeros(2)
yd = np.zeros(2)

xd[1] = x1[1, 2] / x1[1, 1]
xd[0] = (x[0, 2] - x[0, 1] * xd[1]) / x[0, 0]

yd[1] = y1[1, 2] / y1[1, 1]
yd[0] = (y[0, 2] - y[0, 1] * yd[1]) / y[0, 0]
np.set_printoptions(precision=54)
print(xd, yd, sep='\n')
