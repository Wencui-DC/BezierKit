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
            derivative += weights[i] * bernstein.derivative(i, p, u, order) * ctrlpts[i,:]

        return derivative


    @staticmethod
    def wDer(weights, p, u, order):
        derivative = 0.0
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
        [rown, coln] = np.shape(ctrlpts)
        self.p = rown - 1 # bezier curve's degree
        self.weights = np.ones(self.p+1) #default it's all one
        self.dimension = coln # It is either 3 or 2.
        self.sampleSize = 50 # default number of interpolation steps


    # evaluate a bezier curve a single u
    def evaluate(self, u):
        curvePt = 0.0
        sumWeights = 0.0
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

        m = self.p+1
        n = order+1
        wDerRecord = np.zeros([m, 1])
        pStarDerRecrod = np.zeros([m, self.dimension])
        bezierDer = np.zeros([n, self.dimension])
        for i in range(m):
            wDerRecord[i, 0] = bernstein.wDer(self.weights, self.p, u, i)
            bezierDer[i, :] = bernstein.pStarDer(self.weights, self.ctrlpts, self.p, u, i)

        for k in range(n):
            for i in range(k):
                j = k - i
                if j <= self.p:
                    bezierDer[k, :] -= math.comb(k, j) * bezierDer[i,:] * wDerRecord[j, 0]

            bezierDer[k, :] /= wDerRecord[0, 0]

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
                normSquare += firstDer[-1, j] ** 2

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


class basis:
    #@brief: find the index of u_i for u
    #@n: the number of ctrlpts
    #@p: the curve degree
    #@u: the parameter u
    #@U: the knots vector
    @staticmethod
    def findSpan(n, p, u, U):
        if u == U[n + 1]:
            return n

        low = p
        high = n + 1
        index = (low + high) // 2
        while u < U[index] or u >= U[index + 1]:
            if u < U[index]:
                high = index
            else:
                low = index

            index = (low + high) // 2

        return index


    #@brief: calculate the basis function
    #@i: the index of u_i for u
    #@p: the curve degree
    #@u: the parameter u
    #@U: the knots vector
    @staticmethod
    def evaluate(i, p, u, U):
        N = np.ones([p+1, 1])
        left = np.zeros([p, 1])
        right = np.zeros([p, 1])
        for j in range(1, p+1):
            left[j-1, 0] = u - U[i+1-j]
            right[j-1, 0] = U[i+j] - u
            saved = 0.0;
            for r in range(j):
                temp = N[r, 0] / (right[r, 0] + left[j-r-1, 0])
                N[r, 0] = saved + right[r, 0] * temp
                saved = left[j-r-1, 0] * temp

            N[j, 0] = saved

        return N

    @staticmethod
    def derivatives(i, p, u, U):
        m = p+1
        ndu = np.zeros([m, m])
        ndu[0,0] = 1.0
        a = np.zeros([2, m])
        left = np.zeros([m,1])
        right = np.zeros([m,1])
        basisDers = np.zeros([m, m])
        for j in range(1, m):
            left[j,0] = u - U[i+1-j]
            right[j,0] = U[i+j] - u
            saved = 0.0
            for r in range(j):
                ndu[j, r] = right[r+1,0] + left[j-r,0]
                temp = ndu[r, j-1] / ndu[j, r]
                ndu[r, j] = saved + right[r+1,0]*temp
                saved = left[j-r,0]*temp

            ndu[j,j] = saved
        for j in range(m):
            basisDers[0,j] = ndu[j,p]

        #this section computes the derivatives
        for r in range(m):
            s1 = 0
            s2 = 1
            a[0,0] = 1.0
            for k in range(1, m):
                d = 0.0
                rk = r-k
                pk = p-k
                if (r >= k):
                    a[s2,0] = a[s1,0] / ndu[pk+1, rk]
                    d = a[s2,0] * ndu[rk, pk]

                j1 = 1 if rk >= -1 else -rk
                j2 = k-1 if r-1 <= pk else p-r
                for j in range(j1, j2+1):
                    a[s2,j] = (a[s1,j] - a[s1, j-1]) / ndu[pk+1, rk+j]
                    d += a[s2,j] * ndu[rk+j, pk]

                if (r <= pk):
                    a[s2,k] = -a[s1,k-1] / ndu[pk+1, r]
                    d += a[s2,k] * ndu[r,pk]

                basisDers[k,r] = d
                j = s1
                s1 = s2
                s2 = j

        r = p
        for k in range(1, m):
            for j in range(m):
                basisDers[k, j] *= r

            r *= (p-k)

        return basisDers

    @staticmethod
    def wFunDers(span, p, weights, basisDers):
        m = p+1
        wDers = np.zeros([m, 1])
        for i in range(m):
            for j in range(m):
                wDers[i,0] += basisDers[i, j] * weights[span-p+j]

        return wDers

    @staticmethod
    def pStarFunDers(span, p, dimension, weights, ctrlpts, basisDers):
        m = p+1
        pStarDers = np.zeros([m, dimension])
        pStarDers = np.matrix(pStarDers)
        for i in range(m):
            for j in range(m):
                pStarDers[i, :] += basisDers[i,j] * weights[span-p+j] * ctrlpts[span-p+j,:]

        return pStarDers

    @staticmethod
    def calcRationalBSplineAndBezierDers(ders, wDers, p, order):
        for k in range(order+1):
            for m in range(1,k+1):
                if (m <= p):
                    ders[k,:] -= math.comb(k, m) * ders[k-m,:] * wDers[m]

            ders[k,:] /= wDers[0]

        return ders


class nurbs(basis):
    def __init__(self, ctrlpts, knots, weights, degree):
        super().__init__()
        self.ctrlpts = np.matrix(ctrlpts)
        self.U = knots
        self.p = degree
        [rown, coln] = np.shape(ctrlpts)
        self.n = rown - 1
        self.weights = np.matrix(weights).T
        self.dimension = coln
        self.ctrlptsW = np.hstack((np.multiply(self.ctrlpts, self.weights), self.weights))

    def evaluate(self, u):
        i = basis.findSpan(self.n, self.p, u, self.U)
        N = basis.evaluate(i, self.p, u, self.U)
        tempPt = 0
        for j in range(self.p+1):
            tempPt += N[j, 0] * self.ctrlptsW[i-self.p+j, :]

        pt = tempPt[0, 0:-1] / tempPt[0, -1]

        return pt


    def derivative(self, u, order):
        if order < 0:
            raise ValueError('derivative order must be >= 0')

        nurbsDers = np.zeros([order+1, self.dimension])
        span = basis.findSpan(self.n, self.p, u, self.U)
        basisDers = basis.derivatives(span, self.p, u, self.U)
        wDers = basis.wFunDers(span, self.p, self.weights, basisDers)
        nurbsDers[0:self.p+1,:] = \
            basis.pStarFunDers(span, self.p, self.dimension, self.weights, self.ctrlpts, basisDers)
        nurbsDers = basis.calcRationalBSplineAndBezierDers(nurbsDers, wDers, self.p, order)

        return nurbsDers
