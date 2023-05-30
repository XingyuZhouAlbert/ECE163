import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations
from ece163.Controls import VehicleTrim

# Disclaimer: The functions being used in are either from the lecture, the book or the prop cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.

def dThrust_dThrottle(Va, Throttle, epsilon=0.01):
    VAD = VehicleAerodynamicsModel.VehicleAerodynamicsModel()
    Fx,Mx = VAD.CalculatePropForces(Va, Throttle)
    Fx_hat,Mx_hat = VAD.CalculatePropForces(Va, (Throttle+epsilon))
    deltaFx = (Fx_hat - Fx) / epsilon
    return deltaFx

def dThrust_dVa(Va, Throttle, epsilon=0.5):
    VAD = VehicleAerodynamicsModel.VehicleAerodynamicsModel()
    Fx,Mx = VAD.CalculatePropForces(Va, Throttle)
    Fx_hat,Mx_hat = VAD.CalculatePropForces((Va + epsilon), Throttle)
    deltaFx = (Fx_hat - Fx) / epsilon
    return deltaFx

def CreateTransferFunction(trimState, trimInputs):
    TF = Linearized.transferFunctions()
    Va = trimState.Va

    if Va != 0:
        trimState.beta = math.asin(trimState.v / trimState.Va)
    else:
        trimState.beta = math.copysign(math.pi / 2, trimState.u)
    TF.Va_trim = trimState.Va
    TF.alpha_trim = trimState.alpha
    TF.beta_trim = trimState.beta
    TF.gamma_trim = trimState.pitch - trimState.alpha
    TF.theta_trim = trimState.pitch
    TF.phi_trim = trimState.roll

    # equations are obtained from page 69 of the book. EQ 5.23 and
    TF.a_phi1 = -(1/2) * VPC.rho * (Va ** 2) * VPC.S * VPC.b * VPC.Cpp * (VPC.b/ (2 * Va))
    TF.a_phi2 = (1/2) * VPC.rho * (Va ** 2) * VPC.S * VPC.b * VPC.CpdeltaA

    # a_beta 1 and 2 are from page 71, above EQ5.28
    TF.a_beta1 = -((VPC.rho * Va * VPC.S) / (2 * VPC.mass)) * VPC.CYbeta
    TF.a_beta2 = ((VPC.rho * Va * VPC.S) / (2 * VPC.mass)) * VPC.CYdeltaR

    # theta 1-3 are on page 73
    TF.a_theta1 = -((VPC.rho * (Va ** 2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMq * (VPC.c / (2*Va))
    TF.a_theta2 = -((VPC.rho * (Va ** 2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMalpha
    TF.a_theta3 = ((VPC.rho * (Va ** 2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMdeltaE

    # aV1-3 are on supplement pdf page 26
    TF.a_V1 = ((VPC.rho * Va * VPC.S) / VPC.mass) * (VPC.CD0 + (VPC.CDalpha * trimState.alpha) + (VPC.CDdeltaE * trimInputs.Elevator)) - ((1/VPC.mass) * dThrust_dVa(Va, trimInputs.Throttle))
    TF.a_V2 = (1/VPC.mass) * dThrust_dThrottle(Va, trimInputs.Throttle)
    TF.a_V3 = VPC.g0 * math.cos(trimState.pitch - trimState.alpha)

    return TF




