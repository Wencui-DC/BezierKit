import BezierKit

## sample1: create a bezier curve
# 3D bezier
ctrlpts = [0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0]
knots = [0, 0, 0, 0, 1, 1, 1, 1]

bezier = BezierKit.bezier(ctrlpts)
bezier.weights = [2, 1, 1, 1]
u = 0.55
pt = bezier.evaluate(u)
print(pt)


# # show the order-th derivative
order = 2
der = bezier.derivative(u, order)
for i in range(order+1):
    print('The %d-order derivative of C(%.2f) =[%.4f, %.4f, %.4f]' % (i, u, der[i, 0], der[i, 1], der[i, 2]))
