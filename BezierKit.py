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


class bernstein:
    '''bernstein class'''
    @staticmethod
    def bernsteinPoly(i, p, u):
        '''calculate the bernstein polynomial

        @:parameter
        i: the index of u_i for u
        p: the curve degree
        u: the parameter u

        @:return
        bernstein function's value
        '''
        bernsteinFun = 0.
        if 0 <= i <= p:
            bernsteinFun = math.comb(p, i) * math.pow(u, i) * math.pow(1 - u, p - i)

        return bernsteinFun

    @staticmethod
    def derivative(i, p, u, order):
        '''calculate the derivative of bernstein polynomial'''
        derB = 0.0
        newP = p - order
        for j in range(order + 1):
            newI = i - order + j
            if (newI >= 0) and (newI <= newP):
                derB += math.pow(-1, j) * math.comb(order, j) * bernstein.bernsteinPoly(newI, newP, u)

        if derB != 0 and order > 0:
            derB *= math.perm(p, order)

        return derB


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


class bezier(bernstein, visualization):
    '''bezier class, including the rational bezier feature'''
    def __init__(self, ctrlpts):
        super().__init__()
        self.ctrlpts = np.array(ctrlpts) # the control points
        [rown, coln] = np.shape(ctrlpts)
        self.p = rown - 1 # bezier curve's degree
        self._weights = np.ones(rown) #default it's all one
        self.dimension = coln # It is either 3 or 2.
        self.sampleSize = 50 # default number of interpolation steps
        self.__update_ctrlptsW()
    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = np.array(value)
        self.__update_ctrlptsW()

    def __update_ctrlptsW(self):
        rown = self.p + 1
        self.ctrlptsW = np.column_stack((self.ctrlpts * self.weights.reshape(rown, 1), self.weights))

    def __pStarDer(self, u, order):
        pStarDer = 0.0
        for i in range(self.p+1):
            pStarDer += self.weights[i] * bernstein.derivative(i, self.p, u, order) * self.ctrlpts[i,:]

        return pStarDer

    def __wDer(self, u, order):
        wDer = 0.0
        for i in range(self.p+1):
            wDer += self.weights[i] * bernstein.derivative(i, self.p, u, order)

        return wDer

    def evaluate(self, u):
        '''evaluate a bezier curve a single u'''
        curvePt = 0.0
        for i in range(self.p + 1):
            curvePt += self.ctrlptsW[i,:] * bernstein.bernsteinPoly(i, self.p, u)

        pt = curvePt[0:-1] / curvePt[-1]

        return pt

    def derivative(self, u, order):
        '''calculate the order-th derivative of a bezier curve at designated u'''
        if order < 0:
            raise ValueError('derivative order must be >= 0')

        n = order + 1
        m = n if order <= self.p else (self.p + 1)
        wDers = np.zeros([m, 1])
        bezierDers = np.zeros([n, self.dimension])
        for i in range(m):
            wDers[i] = self.__wDer(u, i)
            bezierDers[i, :] = self.__pStarDer(u, i)

        for k in range(n):
            for i in range(k):
                j = k - i
                if j <= self.p:
                    bezierDers[k, :] -= math.comb(k, j) * bezierDers[i,:] * wDers[j]

            bezierDers[k, :] /= wDers[0]

        return bezierDers

    def trace(self):
        '''calculate the trace of a bezier curve

        @:return
        trace: the trace (interpolated points) of the curve
        '''
        step = 1 / (self.sampleSize - 1)
        U = np.arange(0, 1 + step, step)
        trace = np.zeros((self.sampleSize, self.dimension))
        for i in range(self.sampleSize):
            trace[i, :] = self.evaluate(U[i])

        return trace

    def length(self, a=0, b=1):
        '''calculate the arc length of the bezier curve,
        using Legendre-Gauss quadrature

        @:parameter
        a: the lower bound of parameter u \in [0,1]
        b: the higher bound of parameter u \in [0,1]

        @:return
        arcLen: the arc length
        '''
        if a<0 or b>1:
            raise ValueError('Interval of U is not within [0, 1]')

        coef_1 = (b - a) / 2
        coef_2 = (b + a) / 2
        arcLen = 0
        abscissaeLen = 15
        for i in range(abscissaeLen):
            u = coef_1 * abscissae[i] + coef_2
            ders = self.derivative(u, 1)
            firstDer = ders[1, :]
            normSquare = np.sum(firstDer**2)
            arcLen += weights_LG[i] * math.sqrt(normSquare)

        arcLen *= coef_1

        return arcLen

    def curvature(self, u):
        '''calculate the curvature at u'''
        ders = self.derivative(u, 2)
        firstDer = ders[1, :]
        secondDer = ders[2, :]
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
    @staticmethod
    def findSpan(n, p, u, U):
        '''find the index of u_i for u
        @:parameter
        n: the number of ctrlpts
        p: the curve degree
        u: the parameter u
        U: the knots vector

        @:return
        index: the index of u_i for u
        '''
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

    @staticmethod
    def evaluate(i, p, u, U):
        '''evaluate basis function at u
        @:parameter
        i: the index of u_i for u
        p: the curve degree
        u: the parameter u
        U: the knots vector
        @:return
        N: basis function's value at u
        '''
        N = np.ones([p+1, 1])
        left = np.zeros([p, 1])
        right = np.zeros([p, 1])
        for j in range(1, p+1):
            left[j-1] = u - U[i+1-j]
            right[j-1] = U[i+j] - u
            saved = 0.0;
            for r in range(j):
                temp = N[r] / (right[r] + left[j-r-1])
                N[r] = saved + right[r] * temp
                saved = left[j-r-1] * temp

            N[j] = saved

        return N

    @staticmethod
    def derivatives(i, p, u, U, order):
        m = p + 1
        n = (order+1) if order <= p else m
        ndu = np.zeros([m, m])
        ndu[0,0] = 1.0
        a = np.zeros([2, n])
        left = np.zeros([m,1])
        right = np.zeros([m,1])
        basisDers = np.zeros([n, m])
        for j in range(1, m):
            left[j] = u - U[i+1-j]
            right[j] = U[i+j] - u
            saved = 0.0
            for r in range(j):
                ndu[j, r] = right[r+1] + left[j-r]
                temp = ndu[r, j-1] / ndu[j, r]
                ndu[r, j] = saved + right[r+1]*temp
                saved = left[j-r]*temp

            ndu[j,j] = saved

        for j in range(m):
            basisDers[0,j] = ndu[j,p]

        #this section computes the derivatives
        for r in range(m):
            s1 = 0
            s2 = 1
            a[0,0] = 1.0
            for k in range(1, n):
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
        for k in range(1, n):
            for j in range(m):
                basisDers[k, j] *= r
            r *= (p-k)

        return basisDers


class nurbs(basis):
    def __init__(self, ctrlpts, knots, weights, degree):
        super().__init__()
        self.ctrlpts = np.array(ctrlpts)
        self.U = np.array(knots)
        self.p = degree
        rown, coln = np.shape(ctrlpts)
        self.n = rown - 1
        self._weights = np.array(weights)
        self.dimension = coln
        self.sampleSize = 50  # default number of interpolation steps
        self.__update_ctrlptsW()

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = np.array(value)
        self.__update_ctrlptsW()

    def __update_ctrlptsW(self):
        rown = self.n + 1
        self.ctrlptsW = np.column_stack((self.ctrlpts * self.weights.reshape(rown, 1), self.weights))

    def __wDers(self, span, basisDers, k):
        wDers = np.zeros([k, 1])
        for i in range(k):
            for j in range(self.p + 1):
                wDers[i] += basisDers[i, j] * self.weights[span - self.p + j]

        return wDers

    def __pStarDers(self, span, basisDers, k):
        pStarDers = np.zeros([k, self.dimension])
        for i in range(k):
            for j in range(self.p + 1):
                pStarDers[i, :] += basisDers[i, j] * self.weights[span - self.p + j] * self.ctrlpts[span - self.p + j, :]

        return pStarDers

    def __calcRationalBSplineAndBezierDers(self, ders, wDers, order):
        for i in range(order + 1):
            for j in range(1, i + 1):
                if (j <= self.p):
                    ders[i, :] -= math.comb(i, j) * ders[i - j, :] * wDers[j]
            ders[i, :] /= wDers[0]

        return ders

    def evaluate(self, u):
        i = basis.findSpan(self.n, self.p, u, self.U)
        N = basis.evaluate(i, self.p, u, self.U)
        tempPt = 0.
        for j in range(self.p+1):
            tempPt += N[j] * self.ctrlptsW[i-self.p+j, :]

        pt = tempPt[0:-1] / tempPt[-1]

        return pt

    def derivative(self, u, order):
        '''
        @:parameter
        u: the parameter u
        order: the order of the derivative

        @:return: derivatives
        '''
        if order < 0:
            raise ValueError('derivative order must be >= 0')

        k = (order+1) if order <= self.p else (self.p+1)
        nurbsDers = np.zeros([order+1, self.dimension])
        span = basis.findSpan(self.n, self.p, u, self.U)
        basisDers = basis.derivatives(span, self.p, u, self.U, order)
        wDers = self.__wDers(span, basisDers, k)
        nurbsDers[0:k,:] = self.__pStarDers(span, basisDers, k)
        nurbsDers = self.__calcRationalBSplineAndBezierDers(nurbsDers, wDers, order)

        return nurbsDers

    def length(self, a=0, b=1):
        '''calculate the arc length of nurbs curves,
        using Legendre-Gauss quadrature

        @:parameter
        a: the lower bound of parameter u \in [0,1]
        b: the higher bound of parameter u \in [0,1]

        @:return
        arcLen: the arc length
        '''
        if a<0 or b>1:
            raise ValueError('Interval of U is not within [0, 1]')

        coef_1 = (b - a) / 2
        coef_2 = (b + a) / 2
        arcLen = 0
        abscissaeLen = 15
        for i in range(abscissaeLen):
            u = coef_1 * abscissae[i] + coef_2
            ders = self.derivative(u, 1)
            firstDer = ders[1, :]
            normSquare = np.sum(firstDer**2)
            arcLen += weights_LG[i] * math.sqrt(normSquare)

        arcLen *= coef_1

        return arcLen

    def trace(self):
        '''calculate the trace of a nurbs curve

        @:return
        trace: the trace (interpolated points) of the curve
        '''
        step = 1 / (self.sampleSize - 1)
        U = np.arange(0, 1 + step, step)
        trace = np.zeros((self.sampleSize, self.dimension))
        for i in range(self.sampleSize):
            trace[i, :] = self.evaluate(U[i])

        return trace

    def vis(self):
        if self.dimension == 2:
            visualization.plot2d(self.trace(), self.ctrlpts, self.p)
        elif self.dimension == 3:
            visualization.plot3d(self.trace(), self.ctrlpts, self.p)
        else:
            raise ValueError('The nurbs.dimension is neither 2 nor 3!')

        return 0