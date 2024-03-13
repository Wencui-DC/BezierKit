import BezierKit
import numpy as np


ctrlpts = [0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0], [9, 9, 9]
knots = [0, 0, 0, 0, .3, 1, 1, 1, 1]
weights = [2, 1, 1, 1, 3]
degree = 3


u = 0.05
nurbs = BezierKit.nurbs(ctrlpts, knots, weights, degree)

pt = nurbs.evaluate(u)
print(pt)

order = 10
nurbsDers = nurbs.derivative(u, order)
print(nurbsDers)


# i = BezierKit.basis.findSpan(nurbs.n, nurbs.p, u, nurbs.U)
# basisDers = BezierKit.basis.derivatives(i, nurbs.p, u, nurbs.U)

# wDers = BezierKit.basis.wFunDers(i, nurbs.p, nurbs.weights, basisDers)
# print(wDers)

# pStarDers = BezierKit.basis.pStarFunDers(i, nurbs.p, nurbs.dimension, nurbs.weights, nurbs.ctrlpts, basisDers)
# print(pStarDers)
# print(type(pStarDers[0]))


# span = 3
# p =degree
# # a = basisDers[0,0] * ctrlpts[span-p+0,:] * weights[span-p+0, 0]
# j=1
# print( nurbs.weights[span-p+j,:] *  basisDers[0,j] * nurbs.ctrlpts[span-p+j, :])
# # print(basisDers[0,0] * weights[span-p+0, 0]* ctrlpts[span-p+0,:] )