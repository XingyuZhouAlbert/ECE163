import math
import random
from ..Containers import States
from ..Containers import Inputs
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC
# Disclaimer: The functions being used in are either from the lecture, the book or the prop cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.


class WindModel():
    # def __init__(self, Wn = 0.0, We = 0.0, Wd = 0.0, drydenParameters = Inputs.drydenParameters()):
    #     self.wind = States.windState(Wn, We, Wd, 0, 0, 0)
    #     self.drydenParameters = drydenParameters

    def __init__(self, dT=VPC.dT, Va=VPC.InitialSpeed, drydenParamters=VPC.DrydenNoWind):
        self.Wind = States.windState()
        self.dT = dT
        self.Va = Va
        self.drydenParameters = drydenParamters

        self.x_u = [[0]]
        self.x_v = [[0], [0]]
        self.x_w = [[0], [0]]

        self.Phi_u = [[1]]
        self.Gamma_u = [[0]]
        self.H_u = [[1]]

        self.Phi_v = [[1,0] , [0,1]]
        self.Gamma_v = [[0], [0]]
        self.H_v = [[1,1]]

        self.Phi_w = [[1,0] , [0,1]]
        self.Gamma_w = [[0], [0]]
        self.H_w = [[1,1]]

        self.CreateDrydenTransferFns(dT, Va, drydenParamters)

        return

    def CreateDrydenTransferFns(self, dT, Va, drydenParamters):
        if drydenParamters == VPC.DrydenNoWind:
            self.x_u = [[0]]
            self.x_v = [[0], [0]]
            self.x_w = [[0], [0]]

            self.Phi_u = [[1]]
            self.Gamma_u = [[0]]
            self.H_u = [[1]]

            self.Phi_v = [[1, 0], [0, 1]]
            self.Gamma_v = [[0], [0]]
            self.H_v = [[1, 1]]

            self.Phi_w = [[1, 0], [0, 1]]
            self.Gamma_w = [[0], [0]]
            self.H_w = [[1, 1]]
            return
        if math.isclose(Va, 0.0):
            raise ArithmeticError("Attention! Va is near to 0")

        # constant declaration to make the code look cleaner
        Lu = drydenParamters.Lu
        Lv = drydenParamters.Lv
        Lw = drydenParamters.Lw
        sigmau = drydenParamters.sigmau
        sigmav = drydenParamters.sigmav
        sigmaw = drydenParamters.sigmaw

        # for u, on drydenWind cheatsheet page 5
        self.Phi_u[0][0] = math.exp(-(Va/Lu) * dT)
        self.Gamma_u[0][0] = (Lu/Va) * (1 - (math.exp(-(Va/Lu) * dT)))
        self.H_u[0][0] = sigmau * math.sqrt((2*Va)/(math.pi * Lu))

        # for v, on dryden cheatsheet page 5
        self.Phi_v[0][0] = 1 - ((Va/Lv) * dT)
        self.Phi_v[0][1] = -((Va/Lv) ** 2) * dT
        self.Phi_v[1][0] = dT
        self.Phi_v[1][1] = 1 + ((Va/Lv) * dT)
        self.Phi_v = MatrixMath.matrixScalarMultiply(math.exp(-(Va/Lv) * dT), self.Phi_v)

        self.Gamma_v[0][0] = dT
        self.Gamma_v[1][0] = (((Lv/Va) ** 2) * ((math.exp((Va / Lv) * dT)) - 1)) - (Lv / Va) * dT
        self.Gamma_v = MatrixMath.matrixScalarMultiply(math.exp(-(Va/Lv) * dT), self.Gamma_v)

        self.H_v[0][0] = 1
        self.H_v[0][1] = Va / (math.sqrt(3) * Lv)
        self.H_v = MatrixMath.matrixScalarMultiply(sigmav * math.sqrt((3 * Va)/(Lv * math.pi)), self.H_v)

        # for w, on dryden cheatsheet page 5
        self.Phi_w[0][0] = 1 - ((Va/Lw) * dT)
        self.Phi_w[0][1] = -((Va/Lw) ** 2) * dT
        self.Phi_w[1][0] = dT
        self.Phi_w[1][1] = 1 + ((Va / Lw) * dT)
        self.Phi_w = MatrixMath.matrixScalarMultiply(math.exp(-(Va/Lw) * dT), self.Phi_w)

        self.Gamma_w[0][0] = dT
        self.Gamma_w[1][0] = (((Lw/Va) ** 2) * ((math.exp((Va / Lw) * dT)) - 1)) - (Lw / Va) * dT
        self.Gamma_w = MatrixMath.matrixScalarMultiply(math.exp(-(Va/Lw) * dT), self.Gamma_w)

        self.H_w[0][0] = 1
        self.H_w[0][1] = Va / (math.sqrt(3) * Lw)
        self.H_w = MatrixMath.matrixScalarMultiply(sigmaw * math.sqrt((3 * Va)/(Lw * math.pi)), self.H_w)

        return

    def getDrydenTransferFns(self):
        return self.Phi_u, self.Gamma_u, self.H_u, self.Phi_v, self.Gamma_v, self.H_v, self.Phi_w, self.Gamma_w, self.H_w

    def getWind(self):
        return self.Wind

    def reset(self):
        self.dT = VPC.dT
        self.Va = VPC.InitialSpeed
        self.Wind = States.windState()
        self.drydenParameters = VPC.DrydenNoWind
        self.x_u = [[0]]
        self.x_v = [[0], [0]]
        self.x_w = [[0], [0]]
        return

    def setWind(self,windState):
        self.Wind = windState
        return

    def Update(self, uu=None, uv=None, uw=None):
        # generating noise for windmodels
        if uu == None:
            uu = random.gauss(0,1)
        if uv == None:
            uv = random.gauss(0,1)
        if uw == None:
            uw = random.gauss(0,1)

        # updating u
        umat1 = MatrixMath.matrixMultiply(self.Phi_u, self.x_u)
        umat2 = MatrixMath.matrixScalarMultiply(uu, self.Gamma_u)
        self.x_u = MatrixMath.matrixAdd(umat1, umat2)

        # updating v
        vmat1 = MatrixMath.matrixMultiply(self.Phi_v, self.x_v)
        vmat2 = MatrixMath.matrixScalarMultiply(uv, self.Gamma_v)
        self.x_v = MatrixMath.matrixAdd(vmat1, vmat2)

        #updating w
        wmat1 = MatrixMath.matrixMultiply(self.Phi_w, self.x_w)
        wmat2 = MatrixMath.matrixScalarMultiply(uw, self.Gamma_w)
        self.x_w = MatrixMath.matrixAdd(wmat1, wmat2)

        # MatrixMath.matrixPrint(self.x_u)
        # print("")
        # MatrixMath.matrixPrint(self.x_v)
        # print("")
        # MatrixMath.matrixPrint(self.x_w)
        # print("")

        # Generate gusts from state
        # Wu Wv and Ww are still in a matrix, doing [0][0] will make them scalars, which
        # fixed the typeerror generated by running the Aero test.
        self.Wind.Wu = MatrixMath.matrixMultiply(self.H_u, self.x_u)[0][0]
        self.Wind.Wv = MatrixMath.matrixMultiply(self.H_v, self.x_v)[0][0]
        self.Wind.Ww = MatrixMath.matrixMultiply(self.H_w, self.x_w)[0][0]
        return


