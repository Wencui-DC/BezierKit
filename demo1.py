import BezierKit
import math

## sample1: create a bezier curve
# 3D bezier
ctrlpts = ([0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0])
weights = [1, 1/math.sqrt(2), 1, 2]

rationB = BezierKit.bezier(ctrlpts)
rationB.weights = weights
u = 0.1
# pt = rationB.evaluate(u)
# print(pt)

# show the order-th derivative
order = 8
for i in range(order+1):
    der = rationB.derivative(u, i)
    print('The %d-order derivative of C(%.2f) =[%.4f, %.4f, %.4f]' % (i, u, der[0,0], der[0,1], der[0,2]))
