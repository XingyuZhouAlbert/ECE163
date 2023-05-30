import math
from . import MatrixMath
import numpy as np

def dcm2Euler(DCM):
    if DCM[0][2] < -1:
        DCM[0][2] = -1
    elif DCM[0][2] > 1:
        DCM[0][2] = 1

    pitch = -math.asin(DCM[0][2])
    roll = math.atan2(DCM[1][2], DCM[2][2])
    yaw = math.atan2(DCM[0][1], DCM[0][0])

    return yaw, pitch, roll

def euler2DCM(yaw, pitch, roll):
    Ryaw = np.zeros((3,3))
    Rpitch = np.zeros((3,3))
    Rroll = np.zeros((3,3))

    Ryaw[0][0] = math.cos(yaw)
    Ryaw[0][1] = math.sin(yaw)
    Ryaw[1][0] = -math.sin(yaw)
    Ryaw[1][1] = math.cos(yaw)
    Ryaw[2][2] = 1

    Rpitch[0][0] = math.cos(pitch)
    Rpitch[0][2] = -math.sin(pitch)
    Rpitch[1][1] = 1
    Rpitch[2][0] = math.sin(pitch)
    Rpitch[2][2] = math.cos(pitch)

    Rroll[0][0] = 1
    Rroll[1][1] = math.cos(roll)
    Rroll[1][2] = math.sin(roll)
    Rroll[2][1] = -math.sin(roll)
    Rroll[2][2] = math.cos(roll)

    DCM1 = MatrixMath.matrixMultiply(Rroll,Rpitch)
    DCM2 = MatrixMath.matrixMultiply(DCM1,Ryaw)



    # DCM[2][0] is not giving me the correct result

    # DCM[0][0] = math.cos(yaw) * math.cos(pitch)
    # DCM[0][1] = math.sin(yaw) * math.cos(pitch)
    # DCM[0][2] = -math.sin(pitch)
    # DCM[1][0] = math.cos(yaw) * math.sin(pitch) * math.sin(roll) - math.sin(yaw) * math.cos(roll)
    # DCM[1][1] = math.sin(yaw) * math.sin(pitch) * math.sin(roll) + math.cos(yaw) * math.cos(roll)
    # DCM[1][2] = math.cos(pitch) * math.sin(roll)
    # DCM[2][0] = math.cos(yaw) * math.sin(pitch) * math.cos(roll) - math.sin(yaw) * math.sin(roll)
    # DCM[2][1] = math.sin(yaw) * math.sin(pitch) * math.cos(roll) - math.cos(yaw) * math.sin(roll)
    # DCM[2][2] = math.cos(pitch) * math.cos(roll)



    return DCM2


def ned2enu(points):
    transform = np.zeros((3, 3))
    transform[0][1] = 1
    transform[1][0] = 1
    transform[2][2] = -1
    enu = MatrixMath.matrixMultiply(points, transform)
    return enu



















