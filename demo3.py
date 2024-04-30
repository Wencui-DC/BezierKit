import BezierKit
import numpy as np


ctrlpts = [0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0], [9, 9, 9]
knots = [0, 0, 0, 0, .3, 1, 1, 1, 1]
weights = [2, 1, 1, 1, 3]
degree = 3

u = .33
nurbs = BezierKit.nurbs(ctrlpts, knots, weights, degree)

# pt = nurbs.evaluate(u)
# print(pt)

order = 10
# basisDers = BezierKit.basis.derivatives(3, degree, u, knots, order)
# print(basisDers)
nurbsDers = nurbs.derivative(u, order)
for i in range(order+1):
    print("C(%.2f)的%d阶导数 = [%.4f, %.4f, %.4f]" % (u, i, nurbsDers[i][0], nurbsDers[i][1], nurbsDers[i][2]))

