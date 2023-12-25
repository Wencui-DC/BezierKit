import math
import numpy as np
import matplotlib.pyplot as plt

# initiate a bezier curve
class createBezier:
    def __init__(self, ctrlpts):
        self.ctrlpts = np.matrix(ctrlpts)
        self.p = len(ctrlpts)-1
        self.space = len(ctrlpts[0])
        self.sampleSize = 50


# calculate the bernstein polynomial
def bernstein(i, p, u):
    if (i >= 0) and (i <= p):
        return math.comb(p, i) * math.pow(u, i) * math.pow(1-u, p-i)
    else:
        return 0


# calculate the derivative of bernstein polynomial
def bernsteinDer(i, p, u, order):
    derB = 0
    factor = 1
    degree = p - order
    for m in range(order+1):
        j = i - order + m
        if (j >= 0) and (j <= degree):
            derB += math.pow(-1, m) * math.comb(order, m) * bernstein(j, degree, u)

    if (derB != 0) and (order >= 1):
        for j in range(order):
            factor *= (p - j)

    return derB * factor


# calculate the order-th derivative of a bezier curve at designated u
def bezierDer(bezier, u, order):
    if order < 0:
        raise ValueError('Derivative order must be >= 0')

    BezierDer = 0
    for i in range(bezier.p + 1):
        BezierDer += bernsteinDer(i, bezier.p, u, order) * bezier.ctrlpts[i, :]

    return BezierDer


# evaluate a Bezier curve a single u
def evaluateBezier(bezier, u):
    curvePt = 0
    for i in range(bezier.p + 1):
        curvePt += bernstein(i, bezier.p, u) * bezier.ctrlpts[i, :]

    return curvePt


# calculate the trace of a bezier curve
def bezierTrace(bezier):
    step = 1/(bezier.sampleSize-1)
    U = np.arange(0, 1+step, step)
    trace = np.zeros((bezier.sampleSize, bezier.space))
    for i in range(bezier.sampleSize):
        trace[i, :] = evaluateBezier(bezier, U[i])

    return trace


# plot a bezier
def plotBezier(bezier):
    trace = bezierTrace(bezier)
    fig = plt.figure()
    if bezier.space == 2:
        x = np.array(bezier.ctrlpts[:, 0])
        y = np.array(bezier.ctrlpts[:, 1])
        plt.plot(trace[:, 0], trace[:, 1], color='black')
        plt.scatter(x, y, s=50, marker='o', color='blue')
        plt.plot(x, y, linestyle='dashed', color='blue', linewidth=1)
    elif bezier.space == 3:
        x = np.array(bezier.ctrlpts[:, 0])
        y = np.array(bezier.ctrlpts[:, 1])
        z = np.array(bezier.ctrlpts[:, 2])
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(trace[:, 0], trace[:, 1], trace[:, 2], color='black')
        ax.scatter(x, y, z, s=50, marker='o', color='blue')
        ax.plot(x, y, z, linestyle='dashed', color='blue', linewidth=1)
    else:
        raise ValueError('Bezier.space is neither 2 nor 3!')

    plt.title(f'a Bezier curve with degree of {bezier.p}')
    plt.show()

    return 0
