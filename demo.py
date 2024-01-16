import BezierKit

## sample1: create a bezier curve
# 3D bezier
ctrlpts = ([0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0])

# 2D bezier (approximately a unit circle)
# ctrlpts = ([0, 1.00005519], [0.55342686, 0.99873585], [0.99873585, 0.55342686], [1.00005519, 0])


# create a bezier
bezier = BezierKit.bezier(ctrlpts)


## sample2: plot a bezier
bezier.vis()


## sample3: calculate the trace of the bezier
trace = bezier.trace()


## sample4: evaluate the curve at single u
u = 0.23
pt = bezier.evaluate(u)
if bezier.dimension == 3:
    print("C(%.2f) = [%.4f, %.4f, %.4f]" % (u, pt[0, 0], pt[0, 1], pt[0, 2]))
elif bezier.dimension == 2:
    print("C(%.2f) = [%.4f, %.4f]" % (u, pt[0, 0], pt[0, 1]))
else:
    raise ValueError('The dimension of ctrlpts is neither 2 nor 3!')


## sample5: calculate the order-th derivative of a bezier curve
u = 0.23
order = 5
for i in range(order + 1):
    der = bezier.derivative(u, i)
    if bezier.dimension == 3:
        print("The %d-order derivative of C(%.2f) = [%.4f, %.4f, %.4f]" % (i, u, der[0, 0], der[0, 1], der[0, 2]))
    elif bezier.dimension == 2:
        print("The %d-order derivative of C(%.2f) = [%.4f, %.4f]" % (i, u, der[0, 0], der[0, 1]))
    else:
        raise ValueError("The dimension of ctrlpts is neither 2 nor 3!")


## sample6: calculate the length of the curve
Len = bezier.length()
print("The arc length is %.10f" % Len)
# The default parameter a and b are 0 and 1. You can adjust them within [0, 1],
# for example: Len = bezier.length(0.5, 1)


## sample7: calculate the curve's curvature
u = 0.23
k = bezier.curvature(u)
print("The curvature at %.2f is %.4f" % (u, k))
