import math
from ..Containers import States
from ..Containers import Inputs
from ..Modeling import VehicleDynamicsModel
from ..Modeling import WindModel
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC


# Disclaimer: The functions being used in are either from the lecture, the book or the prop cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.


class VehicleAerodynamicsModel():
    def __init__(self, iS = VPC.InitialSpeed, iH = VPC.InitialDownPosition):
        self.vehicleDynamics = VehicleDynamicsModel.VehicleDynamicsModel()
        self.vehicleDynamics.state.pn = VPC.InitialNorthPosition
        self.vehicleDynamics.state.pe = VPC.InitialEastPosition
        self.vehicleDynamics.state.pd = iH
        self.vehicleDynamics.state.u = iS
        self.iS = iS
        self.iH = iH
        self.windModel = WindModel.WindModel()
        return

    def setVehicleState(self,state):
        self.vehicleDynamics.state = state
        return

    def setWindModel(self, Wn=0, We=0, Wd=0, drydenParameters = Inputs.drydenParameters()):
        self.windModel.Wind.Wn = Wn
        self.windModel.Wind.We = We
        self.windModel.Wind.Wd = Wd
        self.drydenParameters = drydenParameters
        return

    def getVehicleState(self):
        return self.vehicleDynamics.state

    def getWindState(self):
        return self.windModel.Wind

    def getVehicleDerivative(self):
        return self.vehicleDynamics.dot

    def setVehicleDerivative(self, dot):
        self.vehicleDynamics.dot = dot
        return

    def gravityForces(self, state):
        Fg = Inputs.forcesMoments()
        Fg.Mx = 0
        Fg.My = 0
        Fg.Mz = 0

        Mg = [[0], [0], [VPC.mass * VPC.g0]]
        g = MatrixMath.matrixMultiply(state.R, Mg)

        Fg.Fx = g[0][0]
        Fg.Fy = g[1][0]
        Fg.Fz = g[2][0]
        return Fg

    def reset(self):
        self.vehicleDynamics = VehicleDynamicsModel.VehicleDynamicsModel()
        self.vehicleDynamics.state.pn = VPC.InitialNorthPosition
        self.vehicleDynamics.state.pe = VPC.InitialEastPosition
        self.vehicleDynamics.state.pd = VPC.InitialDownPosition
        self.vehicleDynamics.state.u = VPC.InitialSpeed
        self.iS = VPC.InitialSpeed
        self.iH = VPC.InitialDownPosition
        self.windModel = WindModel.WindModel()
        return

    def aeroForces(self, state):
        Va = state.Va
        alpha = state.alpha
        beta = state.beta
        forceMo = Inputs.forcesMoments()
        CL_alpha, CD_alpha, CM_alpha = self.CalculateCoeff_alpha(alpha)

        # Need to return empty forceMoment when Va = 0
        if Va == 0:
            return Inputs.forcesMoments()

        #Equation 4.19, constant declaration
        CXa = -CD_alpha * math.cos(alpha) + CL_alpha * math.sin(alpha)
        CXqa = -VPC.CDq * math.cos(alpha) + VPC.CLq * math.sin(alpha)
        CXde = -VPC.CDdeltaE * math.cos(alpha) + VPC.CLdeltaE * math.sin(alpha)
        CZa = -CD_alpha * math.sin(alpha) - CL_alpha * math.cos(alpha)
        CZq = -VPC.CDq * math.sin(alpha) - VPC.CLq * math.cos(alpha)
        CZde = -VPC.CDdeltaE * math.sin(alpha) - VPC.CLdeltaE * math.cos(alpha)

        #Equation 4.18, FORCE
        mat1 = [[CXa + CXqa * (VPC.c/(2*Va)) * state.q] ,
                [VPC.CY0 + VPC.CYbeta * beta + VPC.CYp * (VPC.b/(2*Va)) * state.p + VPC.CYr * (VPC.b/2*Va) * state.r] ,
                [CZa + CZq * (VPC.c/(2*Va)) * state.q]]

        constant1 = 1/2 * VPC.rho * (Va ** 2) * VPC.S
        liftDrag = MatrixMath.matrixScalarMultiply(constant1,mat1)

        #Eq 4.20, MOMENT
        moments1 = [[VPC.b * (VPC.Cl0 + (VPC.Clbeta * beta + VPC.Clp * (VPC.b/(2*Va)) * state.p) + (VPC.Clr * (VPC.b/(2*Va)) * state.r))],
                    [VPC.c * (VPC.CM0 + VPC.CMalpha * alpha + VPC.CMq * (VPC.c/(2*Va)) * state.q)],
                    [VPC.b * (VPC.Cn0 + VPC.Cnbeta * beta + VPC.Cnp * (VPC.b/(2*Va)) * state.p + VPC.Cnr * (VPC.b/(2*Va)) * state.r)]]
        Moments = MatrixMath.matrixScalarMultiply(constant1, moments1)

        forceMo.Fx = liftDrag[0][0]
        forceMo.Fy = liftDrag[1][0]
        forceMo.Fz = liftDrag[2][0]

        forceMo.Mx = Moments[0][0]
        forceMo.My = Moments[1][0]
        forceMo.Mz = Moments[2][0]

        return forceMo

    def CalculateCoeff_alpha(self,alpha):
        #Based on the equations on slide p40
        sigma = (1 + math.exp(-VPC.M * (alpha - VPC.alpha0)) + math.exp(VPC.M * (alpha + VPC.alpha0))) / ((1 + math.exp(-VPC.M * (alpha - VPC.alpha0))) * (1 + math.exp(VPC.M * (alpha + VPC.alpha0))))
        CL_alpha = ((1 - sigma) * (VPC.CL0 + (VPC.CLalpha * alpha))) + (sigma * 2 * math.sin(alpha) * math.cos(alpha))
        CD_alpha = (1 - sigma) * (VPC.CDp + ((CL_alpha * alpha) ** 2) / (math.pi * VPC.AR * VPC.e)) + sigma * (2 * (math.sin(alpha) ** 2))
        CM_alpha = VPC.CM0 + (VPC.CMalpha * alpha)
        return CL_alpha, CD_alpha, CM_alpha

    def CalculatePropForces(self, Va, Throttle):
        # Variable declaration
        Vin = VPC.V_max * Throttle
        KE = 60 / (2 * math.pi * VPC.KV)
        KT = KE

        a = (VPC.rho * (VPC.D_prop ** 5) * VPC.C_Q0) / (4 * (math.pi ** 2))
        b = ((VPC.rho * (VPC.D_prop ** 4) * Va * VPC.C_Q1) / (2*math.pi)) + (KE * KT) / VPC.R_motor
        c = VPC.rho * (VPC.D_prop ** 3) * (Va ** 2) * VPC.C_Q2 - KT * (Vin / VPC.R_motor) + KT * VPC.i0


        # check for imaginary omega
        try:
            omega = (-b + math.sqrt((b ** 2) - 4*a*c)) / (2 * a)
        except:
            omega = 100


        J = (2 * math.pi * Va) / (omega * VPC.D_prop)

        CT = VPC.C_T0 + VPC.C_T1 * J + VPC.C_T2 * (J ** 2)
        CQ = VPC.C_Q0 + VPC.C_Q1 * J + VPC.C_Q2 * (J ** 2)

        # Eq 1 and 2 on prop cheatsheet
        Fx_prop = (VPC.rho * (omega ** 2) * (VPC.D_prop ** 4) * CT) / (4 * (math.pi ** 2))
        Mx_prop = -(VPC.rho * (omega ** 2) * (VPC.D_prop ** 5) * CQ) / (4 * (math.pi ** 2))


        return Fx_prop, Mx_prop

    def controlForces(self, state, controls):
        forceMo = Inputs.forcesMoments()
        Va = state.Va
        alpha = state.alpha

        Fx_prop, Mx_prop = self.CalculatePropForces(Va,controls.Throttle)

        # Eq 4.19 from textbook
        CXde = -VPC.CDdeltaE * math.cos(alpha) + VPC.CLdeltaE * math.sin(alpha)
        CZde = -VPC.CDdeltaE * math.sin(alpha) - VPC.CLdeltaE * math.cos(alpha)

        # Eq 4.18 liftdrag force
        mat1 = [[CXde * controls.Elevator],
                [VPC.CYdeltaA * controls.Aileron + VPC.CYdeltaR * controls.Rudder],
                [CZde * controls.Elevator]]

        constant1 = 1/2 * VPC.rho * (Va ** 2) * VPC.S
        liftDrag = MatrixMath.matrixScalarMultiply(constant1,mat1)

        # Eq 4.20 Moment
        mat2 = [[VPC.b * ((VPC.CldeltaA * controls.Aileron) + (VPC.CldeltaR * controls.Rudder))],
                [VPC.c * (VPC.CMdeltaE * controls.Elevator)],
                [VPC.b * (VPC.CndeltaA * controls.Aileron + VPC.CndeltaR* controls.Rudder)]]

        Moments = MatrixMath.matrixScalarMultiply(constant1, mat2)
         #adding prop force inside control

        forceMo.Fx = liftDrag[0][0] + Fx_prop
        forceMo.Fy = liftDrag[1][0]
        forceMo.Fz = liftDrag[2][0]

        forceMo.Mx = Moments[0][0] + Mx_prop
        forceMo.My = Moments[1][0]
        forceMo.Mz = Moments[2][0]

        return forceMo

    def CalculateAirspeed(self, state, wind):

        windMat1 = [[wind.Wn], [wind.We], [wind.Wd]]
        windMat2 = [[wind.Wu], [wind.Wv], [wind.Ww]]

        # We need azimuth and elevation, and check to see if they r zero or not
        if (wind.Wn == 0 and wind.We == 0 and wind.Wd == 0):
            Wazimuth = 0
            Welevation = 0
        else:
            Wazimuth = math.atan2(wind.We, wind.Wn)
            Welevation = -math.asin(wind.Wd / (math.sqrt((wind.Wn ** 2) + (wind.We ** 2) + (wind.Wd ** 2))))

        # And then the Azimuth-Elevation rotation matrix, on cheatsheet page 2, eq3
        Rae=[[math.cos(Wazimuth) * math.cos(Welevation), math.sin(Wazimuth) * math.cos(Welevation), -math.sin(Welevation)],
            [-math.sin(Wazimuth), math.cos(Wazimuth), 0],
            [math.cos(Wazimuth) * math.sin(Welevation), math.sin(Wazimuth) * math.sin(Welevation), math.cos(Welevation)]]
        # Here we are getting the inertial total wind, also from page 2 Dryden cheatsheet
        RaeT = MatrixMath.matrixTranspose(Rae)
        windMat3 = MatrixMath.matrixMultiply(RaeT, windMat2)
        windMat4 = MatrixMath.matrixAdd(windMat1, windMat3)
        totalIwind = MatrixMath.matrixMultiply(state.R, windMat4)

        ur = state.u - totalIwind[0][0]
        vr = state.v - totalIwind[1][0]
        wr = state.w - totalIwind[2][0]


        Va = math.sqrt((ur ** 2) + (vr ** 2) + (wr ** 2))
        alpha = math.atan2(wr, ur)

        #Added a special ccase
        if math.isclose(Va, 0.0):
            beta = 0.0
        else:
            beta = math.asin(vr / Va)

        return Va, alpha, beta

    def updateForces(self, state, wind, controls):
        # get Va, alpha, beta from calculate airspeed
        Va, alpha, beta = self.CalculateAirspeed(state, wind)
        state.Va = Va
        state.alpha = alpha
        state.beta = beta
        #get all forces
        forceMoCon = self.controlForces(state,controls)
        forceMoAero = self.aeroForces(state)
        forceMoGrav = self.gravityForces(state)

        forceMo = Inputs.forcesMoments()

        forceMo.Fx = forceMoCon.Fx + forceMoAero.Fx + forceMoGrav.Fx
        forceMo.Fy = forceMoCon.Fy + forceMoAero.Fy + forceMoGrav.Fy
        forceMo.Fz = forceMoCon.Fz + forceMoAero.Fz + forceMoGrav.Fz

        forceMo.Mx = forceMoCon.Mx + forceMoAero.Mx + forceMoGrav.Mx
        forceMo.My = forceMoCon.My + forceMoAero.My + forceMoGrav.My
        forceMo.Mz = forceMoCon.Mz + forceMoAero.Mz + forceMoGrav.Mz

        return forceMo

    def getVehicleDynamicsModel(self):
        return self.vehicleDynamics


    def Update(self,controls):
        self.windModel.Update()
        wind = self.windModel.Wind
        state = self.vehicleDynamics.state
        forceMo = self.updateForces(state, wind, controls)
        self.vehicleDynamics.Update(forceMo)

        return















