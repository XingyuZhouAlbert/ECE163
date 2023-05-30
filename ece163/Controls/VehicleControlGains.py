import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Controls
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations

# Disclaimer: The functions being used in are either from the lecture, the book or the prop cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.

def computeGains(tuningParameters=Controls.controlTuning(), linearizedModel=Linearized.transferFunctions()):
    #set up a new controlGain class
    controlGains = Controls.controlGains()
    # Equation from supplement P.36 under EQ6.2, ROLL
    controlGains.kp_roll = (tuningParameters.Wn_roll ** 2) / linearizedModel.a_phi2
    controlGains.ki_roll = 0.001
    controlGains.kd_roll = ((2 * tuningParameters.Zeta_roll * tuningParameters.Wn_roll) - linearizedModel.a_phi1) / linearizedModel.a_phi2

    # Equation from supplement P.38  EQ6.5 and 6.6
    controlGains.kp_course = 2 * tuningParameters.Zeta_course * tuningParameters.Wn_course * (linearizedModel.Va_trim / VPC.g0)
    controlGains.ki_course = (tuningParameters.Wn_course ** 2) * (linearizedModel.Va_trim / VPC.g0)

    # Equation from book 6.15 and 6.17, Page 105
    controlGains.kp_sideslip = ((2 * tuningParameters.Zeta_sideslip * tuningParameters.Wn_sideslip) - linearizedModel.a_beta1) / linearizedModel.a_beta2
    controlGains.ki_sideslip = (1 / linearizedModel.a_beta2) * (((linearizedModel.a_beta1 + (linearizedModel.a_beta2 * controlGains.kp_sideslip)) / (2 * tuningParameters.Zeta_sideslip)) ** 2)

    # Equation from supplement P44 under EQ 6.12
    controlGains.kp_pitch = ((tuningParameters.Wn_pitch) ** 2 - linearizedModel.a_theta2) / linearizedModel.a_theta3
    controlGains.kd_pitch = ((2 * tuningParameters.Zeta_pitch * tuningParameters.Wn_pitch) - linearizedModel.a_theta1) / linearizedModel.a_theta3

    # Equation from supplement P46 EQ 6.13 and 6.14
    Kpitch_DC = (controlGains.kp_pitch * linearizedModel.a_theta3) / (tuningParameters.Wn_pitch ** 2)
    controlGains.ki_altitude = (tuningParameters.Wn_altitude ** 2) / (Kpitch_DC * linearizedModel.Va_trim)
    controlGains.kp_altitude = (2 * tuningParameters.Zeta_altitude * tuningParameters.Wn_altitude) / (Kpitch_DC * linearizedModel.Va_trim)

    # Equation from supplement P 47 EQ 6.15 and 6.16
    controlGains.ki_SpeedfromThrottle = (tuningParameters.Wn_SpeedfromThrottle ** 2) / linearizedModel.a_V2
    controlGains.kp_SpeedfromThrottle = ((2 * tuningParameters.Zeta_SpeedfromThrottle * tuningParameters.Wn_SpeedfromThrottle) - linearizedModel.a_V1) / linearizedModel.a_V2

    # Equation From UAV book P111 EQ 6.27 and 6.28
    controlGains.kp_SpeedfromElevator = (linearizedModel.a_V1 - (2 * tuningParameters.Zeta_SpeedfromElevator * tuningParameters.Wn_SpeedfromElevator)) / (Kpitch_DC * VPC.g0)
    controlGains.ki_SpeedfromElevator = -(tuningParameters.Wn_SpeedfromElevator ** 2) / (Kpitch_DC * VPC.g0)

    return controlGains

def computeTuningParameters(controlGains=Controls.controlGains(), linearizedModel=Linearized.transferFunctions()):
    controlTuning = Controls.controlTuning()
    try:
        # From Page 100 of the UAV book, EQ6.5 and 6.6
        controlTuning.Wn_roll = math.sqrt(controlGains.kp_roll * linearizedModel.a_phi2)
        controlTuning.Zeta_roll = (linearizedModel.a_phi1 + (linearizedModel.a_phi2 * controlGains.kd_roll)) / (2 * controlTuning.Wn_roll)

        # From Page 103 of the UAV book, above EQ 6.12 Course
        controlTuning.Wn_course = math.sqrt((VPC.g0 / linearizedModel.Va_trim) * controlGains.ki_course)
        controlTuning.Zeta_course = ((VPC.g0 / linearizedModel.Va_trim) * controlGains.kp_course) / (2 * controlTuning.Wn_course)

        # From Page 105 of UAV book, EQ 6.14 and 6.15 Sideslip
        controlTuning.Wn_sideslip = math.sqrt(linearizedModel.a_beta2 * controlGains.ki_sideslip)
        controlTuning.Zeta_sideslip = (linearizedModel.a_beta1 + (linearizedModel.a_beta2 * controlGains.kp_sideslip)) / (2 * controlTuning.Wn_sideslip)

        # From page 107 of UAV book, EQ 6.19 and 6.20 PITCH
        controlTuning.Wn_pitch = math.sqrt(linearizedModel.a_theta2 + (controlGains.kp_pitch * linearizedModel.a_theta3))
        controlTuning.Zeta_pitch = (linearizedModel.a_theta1 + (controlGains.kd_pitch * linearizedModel.a_theta3)) / (2 * controlTuning.Wn_pitch)

        # From page 109 of UAV book, above EQ 6.24
        Kpitch_DC = (controlGains.kp_pitch * linearizedModel.a_theta3) / (controlTuning.Wn_pitch ** 2)
        controlTuning.Wn_altitude = math.sqrt(Kpitch_DC * linearizedModel.Va_trim * controlGains.ki_altitude)
        controlTuning.Zeta_altitude = (Kpitch_DC * linearizedModel.Va_trim * controlGains.kp_altitude) / (2 * controlTuning.Wn_altitude)

        # From page 111 of UAV book, above EQ 6.27
        controlTuning.Wn_SpeedfromElevator = math.sqrt(-Kpitch_DC * VPC.g0 * controlGains.ki_SpeedfromElevator)
        controlTuning.Zeta_SpeedfromElevator = (linearizedModel.a_V1 - (Kpitch_DC * VPC.g0 * controlGains.kp_SpeedfromElevator)) / (2 * controlTuning.Wn_SpeedfromElevator)

        # From page 112 of UAV book, above EQ 6.29
        controlTuning.Wn_SpeedfromThrottle = math.sqrt(linearizedModel.a_V2 * controlGains.ki_SpeedfromThrottle)
        controlTuning.Zeta_SpeedfromThrottle = (linearizedModel.a_V1 + (linearizedModel.a_V2 * controlGains.kp_SpeedfromThrottle)) / (2 * controlTuning.Wn_SpeedfromThrottle)
    except ValueError:
        return Controls.controlTuning()

    return controlTuning

