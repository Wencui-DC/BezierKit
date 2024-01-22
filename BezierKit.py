import math
import numpy as np
import matplotlib.pyplot as plt

# global variables for Legendre-Gauss quadrature
abscissae = [0,	-0.201194093997435,	0.201194093997435,
             -0.394151347077563, 0.394151347077563,	-0.570972172608539,
             0.570972172608539,	-0.724417731360170,	0.724417731360170,
             -0.848206583410427, 0.848206583410427, -0.937273392400706,
             0.937273392400706,-0.987992518020485, 0.987992518020485]

weights_LG = [0.202578241925561, 0.198431485327112, 0.198431485327112,
           0.186161000015562, 0.186161000015562, 0.166269205816994,
           0.166269205816994, 0.139570677926154, 0.139570677926154,
           0.107159220467172, 0.107159220467172, 0.0703660474881081,
           0.0703660474881081, 0.0307532419961173, 0.0307532419961173]


# bernstein class
class bernstein:
    @staticmethod
    # calculate the bernstein polynomial 
    def bernsteinPoly(i, p, u):
        if 0 <= i <= p:
            return math.comb(p, i) * math.pow(u, i) * math.pow(1 - u, p - i)
        else:
            return 0


    @staticmethod
    # calculate the derivative of bernstein polynomial
    def derivative(i, p, u, order):
        derB = 0.0
        newP = p - order
        for j in range(order + 1):
            newI = i - order + j
            if (newI >= 0) and (newI <= newP):
                derB += math.pow(-1, j) * math.comb(order, j) * bernstein.bernsteinPoly(newI, newP, u)

        if derB != 0 and order > 0:
            derB *= math.perm(p, order)

        return derB


    @staticmethod
    def pStarDer(weights, ctrlpts, p, u, order):
        derivative = 0.0
        for i in range(p+1):
            derivative += weights[i]* bernstein.derivative(i, p, u, order) * ctrlpts[i,:]

        return derivative


    @staticmethod
    def wDer(weights, p, u, order):
        derivative = 0
        for i in range(p+1):
            derivative += weights[i] * bernstein.derivative(i, p, u, order)

        return derivative

class visualization:
    @staticmethod
    def plot2d(trace, ctrlpts, p):
        plt.figure()
        x = np.array(ctrlpts[:, 0])
        y = np.array(ctrlpts[:, 1])
        plt.plot(trace[:, 0], trace[:, 1], color='blue')
        plt.scatter(x, y, s=50, marker='o', color='black')
        plt.plot(x, y, linestyle='dashed', color='black', linewidth=1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.title(f'a bezier curve with degree of {p}')
        plt.show()

        return 0

    @staticmethod
    def plot3d(trace, ctrlpts, p):
        fig = plt.figure()
        x = np.array(ctrlpts[:, 0])
        y = np.array(ctrlpts[:, 1])
        z = np.array(ctrlpts[:, 2])
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(trace[:, 0], trace[:, 1], trace[:, 2], color='blue')
        ax.scatter(x, y, z, s=50, marker='o', color='black')
        ax.plot(x, y, z, linestyle='dashed', color='black', linewidth=1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_box_aspect([1.0, 1.0, 1.0])
        plt.title(f'a bezier curve with degree of {p}')
        plt.show()

        return 0



# bezier class, including the rational bezier feature
class bezier(bernstein, visualization):
    def __init__(self, ctrlpts):
        super().__init__()
        self.ctrlpts = np.matrix(ctrlpts) # the control points
        self.p = len(ctrlpts)-1 # bezier curve's degree
        self.weights = np.ones(self.p+1) #default it's all one
        self.dimension = len(ctrlpts[0]) # It is either 3 or 2.
        self.sampleSize = 50 # default number of interpolation steps


    # evaluate a bezier curve a single u
    def evaluate(self, u):
        curvePt = 0
        sumWeights = 0
        for i in range(self.p + 1):
            weightTimesBernstein = self.weights[i] * bernstein.bernsteinPoly(i, self.p, u)
            sumWeights += weightTimesBernstein
            curvePt += weightTimesBernstein  * self.ctrlpts[i, :]

        curvePt = curvePt / sumWeights

        return curvePt


    # calculate the order-th derivative of a bezier curve at designated u
    def derivative(self, u, order):
        if order < 0:
            raise ValueError('derivative order must be >= 0')

        bezierDer = 0
        if sum(self.weights) == self.p+1:
            for i in range(self.p + 1):
                bezierDer += bernstein.derivative(i, self.p, u, order) * self.ctrlpts[i, :]
        else:
            for j in range(order):
                i = j+1
                bezierDer += -math.comb(order,i) * self.derivative(u, order-i) * bernstein.wDer(self.weights, self.p, u, i)

            bezierDer += bernstein.pStarDer(self.weights, self.ctrlpts, self.p, u, order)
            bezierDer /= bernstein.wDer(self.weights, self.p, u, 0)

        return bezierDer


    # calculate the trace of a bezier curve
    def trace(self):
        step = 1 / (self.sampleSize - 1)
        U = np.arange(0, 1 + step, step)
        trace = np.zeros((self.sampleSize, self.dimension))
        for i in range(self.sampleSize):
            trace[i, :] = self.evaluate(U[i])

        return trace


    # @brief: calculate the arc length of the bezier curve, using Legendre-Gauss quadrature
    # @a: the lower bound of parameter u \in [0,1]
    # @b: the higher bound of parameter u \in [0,1]
    # @return: the arc length
    def length(self, a=0, b=1):
        if a<0 or b>1:
            raise ValueError('Interval of U is not within [0, 1]')

        coef_1 = (b - a) / 2
        coef_2 = (b + a) / 2
        Len = 0
        abscissaeLen = 15
        for i in range(abscissaeLen):
            u = coef_1 * abscissae[i] + coef_2
            firstDer = self.derivative(u, 1)
            normSquare = 0
            for j in range(self.dimension):
                normSquare += firstDer[0, j] ** 2

            Len += weights_LG[i] * math.sqrt(normSquare)

        return Len * coef_1


    # calculate the curvature at u
    def curvature(self, u):
        firstDer = self.derivative(u, 1)
        secondDer = self.derivative(u, 2)
        k = np.linalg.norm(np.cross(firstDer, secondDer)) / np.linalg.norm(firstDer)**3

        return k


    # plot a bezier
    def vis(self):
        if self.dimension == 2:
            visualization.plot2d(self.trace(), self.ctrlpts, self.p)
        elif self.dimension == 3:
            visualization.plot3d(self.trace(), self.ctrlpts, self.p)
        else:
            raise ValueError('bezier.dimension is neither 2 nor 3!')

        return 0
