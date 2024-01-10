import math
import numpy as np
import matplotlib.pyplot as plt

# global variables
abscissae = [0,	-0.201194093997435,	0.201194093997435,
             -0.394151347077563, 0.394151347077563,	-0.570972172608539,
             0.570972172608539,	-0.724417731360170,	0.724417731360170,
             -0.848206583410427, 0.848206583410427, -0.937273392400706,
             0.937273392400706,-0.987992518020485, 0.987992518020485]

weights = [0.202578241925561, 0.198431485327112, 0.198431485327112,
           0.186161000015562, 0.186161000015562, 0.166269205816994,
           0.166269205816994, 0.139570677926154, 0.139570677926154,
           0.107159220467172, 0.107159220467172, 0.0703660474881081,
           0.0703660474881081, 0.0307532419961173, 0.0307532419961173]


# Bernstein class
class Bernstein:
    @staticmethod
    # calculate the bernstein polynomial 
    def BernsteinPoly(i, p, u):
        if 0 <= i <= p:
            return math.comb(p, i) * math.pow(u, i) * math.pow(1 - u, p - i)
        else:
            return 0

    @staticmethod
    # calculate the derivative of bernstein polynomial
    def Derivative(i, p, u, order):
        derB = 0
        factor = 1
        degree = p - order
        for m in range(order + 1):
            j = i - order + m
            if (j >= 0) and (j <= degree):
                derB += math.pow(-1, m) * math.comb(order, m) * Bernstein.BernsteinPoly(j, degree, u)

        if (derB != 0) and (order >= 1):
            for j in range(order):
                factor *= (p - j)

        return derB * factor


# Bezier class
class Bezier(Bernstein):
    def __init__(self, ctrlpts):
        super().__init__()
        self.ctrlpts = np.matrix(ctrlpts)
        self.p = len(ctrlpts)-1
        self.space = len(ctrlpts[0])
        self.sampleSize = 50


    # evaluate a Bezier curve a single u
    def Evaluate(self, u):
        curvePt = 0
        for i in range(self.p + 1):
            curvePt += Bernstein.BernsteinPoly(i, self.p, u) * self.ctrlpts[i, :]

        return curvePt


    # calculate the order-th derivative of a bezier curve at designated u
    def Derivative(self, u, order):
        if order < 0:
            raise ValueError('Derivative order must be >= 0')

        BezierDer = 0
        for i in range(self.p + 1):
            BezierDer += Bernstein.Derivative(i, self.p, u, order) * self.ctrlpts[i, :]

        return BezierDer

    # calculate the trace of a bezier curve
    def Trace(self):
        step = 1 / (self.sampleSize - 1)
        U = np.arange(0, 1 + step, step)
        trace = np.zeros((self.sampleSize, self.space))
        for i in range(self.sampleSize):
            trace[i, :] = self.Evaluate(U[i])

        return trace


    # @brief: calculate the arc length of the bezier curve
    # @a: the lower bound of parameter u \in [0,1]
    # @b: the higher bound of parameter u \in [0,1]
    # @return: the arc length
    def Length(self, a, b):
        if a<0 and b>1:
            raise ValueError('Interval of U is not within [0, 1]')

        coef_1 = (b - a) / 2
        coef_2 = (b + a) / 2
        Len = 0
        abscissaeLen = 15
        for i in range(abscissaeLen):
            u = coef_1 * abscissae[i] + coef_2
            firstDer = self.Derivative(u, 1)
            normSquare = 0
            for j in range(self.space):
                normSquare += firstDer[0, j] ** 2

            Len += weights[i] * math.sqrt(normSquare)

        return Len * coef_1


    # plot a bezier
    def Plot(self):
        trace = self.Trace()
        fig = plt.figure()
        if self.space == 2:
            x = np.array(self.ctrlpts[:, 0])
            y = np.array(self.ctrlpts[:, 1])
            plt.plot(trace[:, 0], trace[:, 1], color='black')
            plt.scatter(x, y, s=50, marker='o', color='blue')
            plt.plot(x, y, linestyle='dashed', color='blue', linewidth=1)
        elif self.space == 3:
            x = np.array(self.ctrlpts[:, 0])
            y = np.array(self.ctrlpts[:, 1])
            z = np.array(self.ctrlpts[:, 2])
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(trace[:, 0], trace[:, 1], trace[:, 2], color='black')
            ax.scatter(x, y, z, s=50, marker='o', color='blue')
            ax.plot(x, y, z, linestyle='dashed', color='blue', linewidth=1)
        else:
            raise ValueError('Bezier.space is neither 2 nor 3!')

        plt.title(f'a Bezier curve with degree of {self.p}')
        plt.show()

        return 0
