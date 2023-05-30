import math
import numpy as np
from ..Containers import States
from ..Utilities import MatrixMath
from ..Utilities import Rotations
from ..Constants import VehiclePhysicalConstants as VPC


# Disclaimer: The functions being used in derivative, Rexp, forwardEuler and Integratestate are coming from the
# attitude cheatsheet. I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.


class VehicleDynamicsModel():
    def __init__(self, dT=VPC.dT):
        self.dT = dT
        self.state = States.vehicleState()
        self.dot = States.vehicleState()
        return

    def getVehicleState(self):
        return self.state

    def reset(self):
        self.state = States.vehicleState()
        self.dot = States.vehicleState()
        return

    def resetVehicleState(self):
        self.state = States.vehicleState()
        return self.state

    def setVehicleState(self, state):
        self.state = state
        return

    def derivative(self, state, forcesMoments):

        # for position derivative
        DCM = Rotations.euler2DCM(state.yaw, state.pitch, state.roll)
        gSpeed = [[state.u], [state.v], [state.w]]
        position = MatrixMath.matrixMultiply(MatrixMath.matrixTranspose(DCM), gSpeed)
        pndot = position[0][0]
        pedot = position[1][0]
        pddot = position[2][0]

        # for groundspeed derivative
        gSpeed1 = [[(state.r * state.v) - (state.q * state.w)], [(state.p * state.w) - (state.r * state.u)],
                   [(state.q * state.u) - (state.p * state.v)]]
        forceMat = [[forcesMoments.Fx], [forcesMoments.Fy], [forcesMoments.Fz]]
        forceMat2 = MatrixMath.matrixScalarMultiply(1 / VPC.mass, forceMat)
        groundSpeed = MatrixMath.matrixAdd(gSpeed1, forceMat2)
        udot = groundSpeed[0][0]
        vdot = groundSpeed[1][0]
        wdot = groundSpeed[2][0]

        # for Euler Angles
        rpymat = np.zeros((3, 3))
        rpymat[0][0] = 1
        rpymat[0][1] = math.sin(state.roll) * math.tan(state.pitch)
        rpymat[0][2] = math.cos(state.roll) * math.tan(state.pitch)
        rpymat[1][1] = math.cos(state.roll)
        rpymat[1][2] = -math.sin(state.roll)
        rpymat[2][1] = math.sin(state.roll) / math.cos(state.pitch)
        rpymat[2][2] = math.cos(state.roll) / math.cos(state.pitch)
        bodyRateMat = [[state.p], [state.q], [state.r]]
        rypmatDot = MatrixMath.matrixMultiply(rpymat, bodyRateMat)
        rolldot = rypmatDot[0][0]
        pitchdot = rypmatDot[1][0]
        yawdot = rypmatDot[2][0]

        # for Body Rates
        p = state.p
        q = state.q
        r = state.r

        Gamma1 = VPC.Gamma1
        Gamma2 = VPC.Gamma2
        Gamma7 = VPC.Gamma7
        Gamma3 = VPC.Jzz / VPC.Jdet
        Gamma4 = VPC.Jxz / VPC.Jdet
        Gamma5 = (VPC.Jzz - VPC.Jxx) / VPC.Jyy
        Gamma6 = VPC.Jxz / VPC.Jyy
        Gamma8 = VPC.Jxx / VPC.Jdet

        l = forcesMoments.Mx
        m = forcesMoments.My
        n = forcesMoments.Mz

        BodyMat1 = [[Gamma1 * p * q - Gamma2 * q * r], [Gamma5 * p * r - Gamma6 * (p ** 2 - r ** 2)],
                    [(Gamma7 * p * q) - Gamma1 * q * r]]
        BodyMat2 = [[Gamma3 * l + Gamma4 * n], [(1 / VPC.Jyy) * m], [Gamma4 * l + Gamma8 * n]]
        BodyMatDot = MatrixMath.matrixAdd(BodyMat1, BodyMat2)

        pdot = BodyMatDot[0][0]
        qdot = BodyMatDot[1][0]
        rdot = BodyMatDot[2][0]

        Rskew = MatrixMath.matrixScalarMultiply(-1, MatrixMath.matrixMultiply(MatrixMath.matrixSkew(p, q, r), DCM))
        #dot.R = Rskew
        # Since we should not change the internal variable in derivative
        dot = States.vehicleState(pndot, pedot, pddot, udot, vdot, wdot, yawdot, pitchdot, rolldot, pdot, qdot, rdot)
        dot.R = Rskew

        return dot

    def Rexp(self, dT, state, dot):
        I = np.zeros((3, 3))
        I[0][0] = 1
        I[1][1] = 1
        I[2][2] = 1

        pkD = dot.p
        qkD = dot.q
        rkD = dot.r

        pk = state.p
        qk = state.q
        rk = state.r

        bodyAngularD = [[pkD], [qkD], [rkD]]
        bodyAngular = [[pk], [qk], [rk]]

        w = MatrixMath.matrixAdd(bodyAngular, MatrixMath.matrixScalarMultiply(dT / 2, bodyAngularD))

        p = w[0][0]
        q = w[1][0]
        r = w[2][0]

        magw = math.hypot(p, q, r)

        if magw < 0.2:
            num1 = dT - ((dT ** 3) * (magw ** 2)) / 6 + ((dT ** 5) * (magw ** 4)) / 120
            num2 = (dT ** 2) / 2 - ((dT ** 4) * (magw ** 2)) / 24 + ((dT ** 6) * (magw ** 4)) / 720
        else:
            num1 = math.sin(magw * dT) / magw
            num2 = (1 - math.cos(magw * dT)) / (magw ** 2)

        skew = MatrixMath.matrixSkew(p, q, r)
        skewsquare = MatrixMath.matrixMultiply(skew, skew)
        mat1 = MatrixMath.matrixScalarMultiply(num1, skew)
        mat2 = MatrixMath.matrixScalarMultiply(num2, skewsquare)
        mat3 = MatrixMath.matrixSubtract(I, mat1)
        matExp = MatrixMath.matrixAdd(mat3, mat2)

        Rskew = MatrixMath.matrixScalarMultiply(-1, MatrixMath.matrixMultiply(MatrixMath.matrixSkew(p, q, r), state.R))
        self.dot.R = Rskew

        return matExp

    def IntegrateState(self, dT, state, dot):
        pn = state.pn + (dot.pn * dT)
        pe = state.pe + (dot.pe * dT)
        pd = state.pd + (dot.pd * dT)

        u = state.u + (dot.u * dT)
        v = state.v + (dot.v * dT)
        w = state.w + (dot.w * dT)

        p = state.p + (dot.p * dT)
        q = state.q + (dot.q * dT)
        r = state.r + (dot.r * dT)

        matExp = self.Rexp(dT, state, dot)
        newR = MatrixMath.matrixMultiply(matExp, state.R)
        # # this is to update the DCM
        # R = newR

        euler = Rotations.dcm2Euler(newR)
        yaw = euler[0]
        pitch = euler[1]
        roll = euler[2]
        updatedstate = States.vehicleState(pn, pe, pd, u, v, w, yaw, pitch, roll, p, q, r, newR)
        updatedstate.Va = state.Va
        updatedstate.alpha = state.alpha
        updatedstate.beta = state.beta
        updatedstate.chi = state.chi
        return updatedstate

    def ForwardEuler(self, forcesMoments):
        dot = self.derivative(self.state, forcesMoments)
        state = self.getVehicleState()
        state2 = self.IntegrateState(self.dT, state, dot)
        return state2

    def Update(self, forcesMoments):
        self.state = self.ForwardEuler(forcesMoments)
        # self.setVehicleState(state)
