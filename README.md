# Examples: more details please see demo.py

## sample1: create a bezier curve
import BezierKit

ctrlpts = ([0, 2, 15], [2, 2, 10], [2, 0, 5], [0, 0, 0])

bezier = BezierKit.bezier(ctrlpts)

## sample2: plot a bezier
bezier.plot()

![img.png](demoPic.png)

## sample3: evaluate the curve at single parameter u
u = 0.1

pt = bezier.evaluate(u)

print(pt)

## sample4: calculate the trace of the bezier
trace = bezier.trace()

## sample5: calculate the order-th derivative of a bezier curve
order = 5

der = bezier.derivative(u, order)

print(der)

## sample6: calculate the arc length
Len = bezier.length()

print("The arc length is %.10f" % Len)

