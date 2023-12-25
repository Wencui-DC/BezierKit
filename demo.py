import myBezier

## sample1: create a bezier curve
# 3D bezier
ctrlpts = ([0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0])
# 2D bezier
#ctrlpts = ([0, 2], [2, 2], [2, 0])
bezier = myBezier.createBezier(ctrlpts)

## sample2: plot a bezier
myBezier.plotBezier(bezier)

## sample3: evaluate the curve at single u
u = 0.1
pt = myBezier.evaluateBezier(bezier, u)
if bezier.space == 3:
    print("C(%.2f) = [%.2f, %.2f, %.2f]" % (u, pt[0, 0], pt[0, 1], pt[0, 2]))
elif bezier.space == 2:
    print("C(%.2f) = [%.2f, %.2f]" % (u, pt[0, 0], pt[0, 1]))
else:
    raise ValueError('The dimension of ctrlpts is neither 2 nor 3!')

## sample4: calculate the trace of the bezier
trace = myBezier.bezierTrace(bezier)

## sample5: calculate the order-th derivative of a bezier curve
order = 5
for i in range(order + 1):
    der = myBezier.bezierDer(bezier, u, i)
    if bezier.space == 3:
        print("C(%.2f)的%d阶导数 = [%.4f, %.4f, %.4f]" % (u, i, der[0, 0], der[0, 1], der[0, 2]))
    elif bezier.space == 2:
        print("C(%.2f)的%d阶导数 = [%.4f, %.4f]" % (u, i, der[0, 0], der[0, 1]))
    else:
        raise ValueError("The dimension of ctrlpts is neither 2 nor 3!")
