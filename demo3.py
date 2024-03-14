import BezierKit
import numpy as np


ctrlpts = [0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0], [9, 9, 9]
knots = [0, 0, 0, 0, .3, 1, 1, 1, 1]
weights = [2, 1, 1, 1, 3]
degree = 3


u = 0.59
nurbs = BezierKit.nurbs(ctrlpts, knots, weights, degree)

pt = nurbs.evaluate(u)
print(pt)

order = 10
nurbsDers = nurbs.derivative(u, order)
print(nurbsDers)
