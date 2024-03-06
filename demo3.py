import BezierKit
import numpy as np


ctrlpts = np.array([
    [0, 1, 2], [1, 3, 3], [4, 2, 4], [6, 0, 5]
])
knots = np.array([0, 0, 0, 0, 1, 1, 1, 1])
weights = [1, 1, 1, 2]
degree = 3

u = 0.02
nurbs = BezierKit.nurbs(ctrlpts, knots, weights, degree)
# i = BezierKit.basis.findSpan(nurbs.n, degree, u, knots)
# print(i)

# Basis = BezierKit.basis.basisFun(i, u, degree, knots)
# print(Basis)

pt = nurbs.evaluate(u)
print(pt)
