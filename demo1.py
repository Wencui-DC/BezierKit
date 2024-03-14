import BezierKit
import math
import numpy as np
## sample1: create a bezier curve
# 3D bezier
ctrlpts = [0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0]
knots = [0, 0, 0, 0, 1, 1, 1, 1]

rationB = BezierKit.bezier(ctrlpts)
rationB.weights = [2, 1, 1, 1]
u = 0.59
# pt = rationB.evaluate(u)
# print(pt)

# order = 8
# wDerRecord = np.zeros([order, 1])
# pStarDerRecrod = np.zeros([order, rationB.dimension])
# for j in range(order):
#     i = j + 1
#     wDerRecord[j, 0] = BezierKit.bernstein.wDer(rationB.weights, rationB.p, u, i)
#     pStarDerRecrod[j, :] = BezierKit.bernstein.pStarDer(rationB.weights, rationB.ctrlpts, rationB.p, u, i)
#
# print(pStarDerRecrod)
# print(wDerRecord)

# # show the order-th derivative
order = 10
der = rationB.derivative(u, order)
for i in range(order+1):
    print('The %d-order derivative of C(%.2f) =[%.4f, %.4f, %.4f]' % (i, u, der[i, 0], der[i, 1], der[i, 2]))
