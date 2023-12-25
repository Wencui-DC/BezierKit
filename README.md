# Examples: more details please see demo.py

## sample1: create a bezier curve
import myBezier

ctrlpts = ([0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0])

bezier = myBezier.createBezier(ctrlpts)

## sample2: plot a bezier
myBezier.plotBezier(bezier)

## sample3: evaluate the curve at single u
u = 0.1

pt = myBezier.evaluateBezier(bezier, u)

## sample4: calculate the trace of the bezier
trace = myBezier.bezierTrace(bezier)

## sample5: calculate the order-th derivative of a bezier curve
order = 5

der = myBezier.bezierDer(bezier, u, order)

