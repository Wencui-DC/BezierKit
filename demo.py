import BezierKit

## sample1: create a bezier curve
# 3D bezier
ctrlpts = ([0, 1.00005519, 1], [0.55342686, 0.99873585, 2], [0.99873585, 0.55342686, 3], [1.00005519, 0, 4])

# 2D bezier (approximately a circle)
# ctrlpts = ([0, 1.00005519], [0.55342686, 0.99873585], [0.99873585, 0.55342686], [1.00005519, 0])

# create a bezier
bezier = BezierKit.Bezier(ctrlpts)


# ## sample2: plot a bezier
# bezier.Plot()
#
#
# ## sample3: calculate the trace of the bezier
# trace = bezier.Trace()
#
#
# ## sample4: evaluate the curve at single u
# u = 0.23
# pt = bezier.Evaluate(u)
# if bezier.space == 3:
#     print("C(%.2f) = [%.2f, %.2f, %.2f]" % (u, pt[0, 0], pt[0, 1], pt[0, 2]))
# elif bezier.space == 2:
#     print("C(%.2f) = [%.2f, %.2f]" % (u, pt[0, 0], pt[0, 1]))
# else:
#     raise ValueError('The dimension of ctrlpts is neither 2 nor 3!')
#
#
# ## sample5: calculate the order-th derivative of a bezier curve
# order = 5
# for i in range(order + 1):
#     der = bezier.Derivative(u, i)
#     if bezier.space == 3:
#         print("C(%.2f)的%d阶导数 = [%.4f, %.4f, %.4f]" % (u, i, der[0, 0], der[0, 1], der[0, 2]))
#     elif bezier.space == 2:
#         print("C(%.2f)的%d阶导数 = [%.4f, %.4f]" % (u, i, der[0, 0], der[0, 1]))
#     else:
#         raise ValueError("The dimension of ctrlpts is neither 2 nor 3!")


## sample6: calculate the length of the curve
Len = bezier.Length(0,1)
print("The arc length of the bezier curve is %.10f" % Len)
