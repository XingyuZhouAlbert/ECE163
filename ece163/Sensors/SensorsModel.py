import math
import random
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Utilities import MatrixMath
from ..Containers import Sensors
from ..Constants import VehiclePhysicalConstants as VPC
from ..Constants import VehicleSensorConstants as VSC
from ..Modeling import VehicleAerodynamicsModel

# Disclaimer: The functions being used in are either from the lecture, the book or the sensor cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.


class GaussMarkov():
    def __init__(self, dT=VPC.dT, tau=1e6, eta=0.0):
        self.dT = dT
        self.tau = tau
        self.eta = eta
        self.v = 0.0
        return

    def update(self, vnoise=None):
        if (vnoise == None):
            vnoise = random.gauss(0, self.eta)

        self.v = (math.exp((-self.dT / self.tau)) * self.v) + vnoise
        return self.v

    def reset(self):
        self.v = 0.0
        return


class GaussMarkovXYZ():
    def __init__(self, dT=VPC.dT, tauX=1e6, etaX=0.0, tauY=1e6, etaY=0.0, tauZ=1e6, etaZ=0.0):
        self.dT = dT
        self.tauX = tauX
        self.etaX = etaX
        self.tauY = tauY
        self.etaY = etaY
        self.tauZ = tauZ
        self.etaZ = etaZ

        self.gmX = GaussMarkov(dT, tauX, etaX)
        self.gmY = GaussMarkov(dT, tauY, etaY)
        self.gmZ = GaussMarkov(dT, tauZ, etaZ)
        self.vX = 0.0
        self.vY = 0.0
        self.vZ = 0.0
        return

    def update(self, vXnoise=None, vYnoise=None, vZnoise=None):
        if vXnoise == None:
            vXnoise = random.gauss(0, self.gmX.eta)
        if vYnoise == None:
            vYnoise = random.gauss(0, self.gmY.eta)
        if vZnoise == None:
            vZnoise = random.gauss(0, self.gmZ.eta)

        self.vX = self.gmX.update(vXnoise)
        self.vY = self.gmY.update(vYnoise)
        self.vZ = self.gmZ.update(vZnoise)

        return self.vX, self.vY, self.vZ

    def reset(self):
        self.gmX.reset()
        self.gmY.reset()
        self.gmZ.reset()
        self.vX = 0.0
        self.vY = 0.0
        self.vZ = 0.0
        return


class SensorsModel():
    def __init__(self, aeroModel=VehicleAerodynamicsModel.VehicleAerodynamicsModel(), taugyro=VSC.gyro_tau, etagyro=VSC.gyro_eta, tauGPS=VSC.GPS_tau, etaGPSHorizontal=VSC.GPS_etaHorizontal, etaGPSVertical=VSC.GPS_etaVertical, gpsUpdateHz=VSC.GPS_rate):
        self.aeroModel = aeroModel
        self.taugyro = taugyro
        self.etagyro = etagyro
        self.tauGPS = tauGPS
        self.etaGPSHorizontal = etaGPSHorizontal
        self.etaGPSVertical = etaGPSVertical
        self.gpsUpdateHz = gpsUpdateHz
        self.sensorTrue = Sensors.vehicleSensors()
        self.sensorNoisy = Sensors.vehicleSensors()
        self.sensorBias = self.initializeBiases()
        self.sensorSigma = self.initializeSigmas()
        vehicleDynamics = aeroModel.getVehicleDynamicsModel()
        self.dT = vehicleDynamics.dT
        self.gyroGM = GaussMarkovXYZ(self.dT, taugyro, etagyro)
        self.gpsGM = GaussMarkovXYZ(1/gpsUpdateHz, tauGPS, etaGPSHorizontal, tauGPS, etaGPSHorizontal, tauGPS, etaGPSVertical)
        self.updateTicks = 0
        self.gpsTickUpdate = 1 / (self.dT * gpsUpdateHz)
        return

    #Wrapper function to return the noisy sensor values
    def getSensorsNoisy(self):
        return self.sensorNoisy
    # Wrapper function to return the true sensor values
    def getSensorsTrue(self):
        return self.sensorTrue

    # Function to generate the biases for each of the sensors.
    def initializeBiases(self, gyroBias=VSC.gyro_bias, accelBias=VSC.accel_bias, magBias=VSC.mag_bias, baroBias=VSC.baro_bias, pitotBias=VSC.pitot_bias):

        sensorBiases = Sensors.vehicleSensors()

        sensorBiases.gyro_x = gyroBias * random.uniform(-1,1)
        sensorBiases.gyro_y = gyroBias * random.uniform(-1,1)
        sensorBiases.gyro_z = gyroBias * random.uniform(-1, 1)

        sensorBiases.accel_x = accelBias * random.uniform(-1, 1)
        sensorBiases.accel_y = accelBias * random.uniform(-1, 1)
        sensorBiases.accel_z = accelBias * random.uniform(-1, 1)

        sensorBiases.mag_x = magBias * random.uniform(-1, 1)
        sensorBiases.mag_y = magBias * random.uniform(-1, 1)
        sensorBiases.mag_z = magBias * random.uniform(-1, 1)

        sensorBiases.baro= baroBias * random.uniform(-1, 1)

        sensorBiases.pitot = pitotBias * random.uniform(-1, 1)

        sensorBiases.gps_e = 0.0
        sensorBiases.gps_cog = 0.0
        sensorBiases.gps_n = 0.0
        sensorBiases.gps_alt = 0.0
        sensorBiases.gps_sog = 0.0

        return sensorBiases

    # Function to gather all of the white noise standard deviations into a single vehicleSensor class object. These will be used as the input to generating the white noise added to each sensor when generating the noisy sensor data.
    def initializeSigmas(self, gyroSigma=VSC.gyro_sigma, accelSigma=VSC.accel_sigma, magSigma=VSC.mag_sigma, baroSigma=VSC.baro_sigma, pitotSigma=VSC.pitot_sigma, gpsSigmaHorizontal=VSC.GPS_sigmaHorizontal, gpsSigmaVertical=VSC.GPS_sigmaVertical, gpsSigmaSOG=VSC.GPS_sigmaSOG, gpsSigmaCOG=VSC.GPS_sigmaCOG):

        sensorSigma = Sensors.vehicleSensors()

        sensorSigma.gyro_x = gyroSigma
        sensorSigma.gyro_y = gyroSigma
        sensorSigma.gyro_z = gyroSigma

        sensorSigma.accel_x = accelSigma
        sensorSigma.accel_y = accelSigma
        sensorSigma.accel_z = accelSigma

        sensorSigma.mag_x = magSigma
        sensorSigma.mag_y = magSigma
        sensorSigma.mag_z = magSigma

        sensorSigma.baro= baroSigma

        sensorSigma.pitot = pitotSigma

        sensorSigma.gps_e = gpsSigmaHorizontal
        sensorSigma.gps_cog = gpsSigmaCOG
        sensorSigma.gps_n = gpsSigmaHorizontal
        sensorSigma.gps_alt = gpsSigmaVertical
        sensorSigma.gps_sog = gpsSigmaSOG

        return sensorSigma

    # Function to update the accelerometer sensor. Will be called within the updateSensors functions.
    def updateGPSTrue(self, state, dot):
        gps_e = state.pe
        gps_n = state.pn
        gps_alt = -state.pd
        gps_sog = math.hypot(state.u, state.v, state.w)
        gps_cog = math.atan2(dot.pe, dot.pn)
        return gps_n, gps_e, gps_alt, gps_sog, gps_cog
    # Function to update the magnetometer sensor. Will be called within the updateSensors functions.
    def updateMagsTrue(self, state):
        magTrue = MatrixMath.matrixMultiply(state.R, VSC.magfield)
        return magTrue[0][0], magTrue[1][0], magTrue[2][0]

    # Function to update the rate gyro sensor. Will be called within the updateSensors functions.
    def updateGyrosTrue(self, state):
        return state.p, state.q, state.r

    # Function to update the accelerometer sensor. Will be called within the updateSensors functions.
    def updateAccelsTrue(self, state, dot):
        # Equation coming from the UAV handbook page 122
        accel_x = dot.u + (state.q * state.w) - (state.r * state.v) + (VPC.g0 * math.sin(state.pitch))
        accel_y = dot.v + (state.r * state.u) - (state.p * state.w) - (VPC.g0 * math.cos(state.pitch) * math.sin(state.roll))
        accel_z = dot.w + (state.p * state.v) - (state.q * state.u) - (VPC.g0 * math.cos(state.pitch) * math.cos(state.roll))
        return accel_x, accel_y, accel_z

    # Function to update the pressure sensors onboard the aircraft. Will be called within the updateSensors functions.
    def updatePressureSensorsTrue(self, state):
        Pbaro = VPC.rho * VPC.g0 * (state.pd) + VSC.Pground
        Ppitot = VPC.rho * ((state.Va ** 2) / 2)
        return Pbaro, Ppitot

    # Function to generate the true sensors given the current state and state derivative.
    def updateSensorsTrue(self, prevTrueSensors, state, dot):
        sensorTrue = Sensors.vehicleSensors()
        sensorTrue.gyro_x, sensorTrue.gyro_y, sensorTrue.gyro_z = self.updateGyrosTrue(state)
        sensorTrue.accel_x, sensorTrue.accel_y, sensorTrue.accel_z = self.updateAccelsTrue(state,dot)
        sensorTrue.mag_x, sensorTrue.mag_y, sensorTrue.mag_z = self.updateMagsTrue(state)
        sensorTrue.baro, sensorTrue.pitot = self.updatePressureSensorsTrue(state)

        # Check to see if the GPS needs updating
        if (self.updateTicks % self.gpsTickUpdate):
            sensorTrue.gps_n, sensorTrue.gps_e, sensorTrue.gps_alt, sensorTrue.gps_sog, sensorTrue.gps_cog = self.updateGPSTrue(state,dot)
        else:
            sensorTrue.gps_n = prevTrueSensors.gps_n
            sensorTrue.gps_e = prevTrueSensors.gps_e
            sensorTrue.gps_alt = prevTrueSensors.gps_alt
            sensorTrue.gps_sog = prevTrueSensors.gps_sog
            sensorTrue.gps_cog = prevTrueSensors.gps_cog

        return sensorTrue
    # Function to generate the noisy sensor data given the true sensor readings, the biases, and the sigmas for the white noise on each sensor.
    def updateSensorsNoisy(self, trueSensors=Sensors.vehicleSensors(), noisySensors=Sensors.vehicleSensors(), sensorBiases=Sensors.vehicleSensors(), sensorSigmas=Sensors.vehicleSensors()):
        sensorNoise = Sensors.vehicleSensors()
        # Getting the drift for gyro
        gx,gy,gz = self.gyroGM.update()

        sensorNoise.gyro_x = trueSensors.gyro_x + sensorBiases.gyro_x + gx + random.gauss(0, sensorSigmas.gyro_x)
        sensorNoise.gyro_y = trueSensors.gyro_y + sensorBiases.gyro_y + gy + random.gauss(0, sensorSigmas.gyro_y)
        sensorNoise.gyro_z = trueSensors.gyro_z + sensorBiases.gyro_z + gz + random.gauss(0, sensorSigmas.gyro_z)

        # For acceleration theres no drift
        sensorNoise.accel_x = trueSensors.accel_x + sensorBiases.accel_x + random.gauss(0, sensorSigmas.accel_x)
        sensorNoise.accel_y = trueSensors.accel_y + sensorBiases.accel_y + random.gauss(0, sensorSigmas.accel_y)
        sensorNoise.accel_z = trueSensors.accel_z + sensorBiases.accel_z + random.gauss(0, sensorSigmas.accel_z)

        # Same for mag, no drift as well
        sensorNoise.mag_x = trueSensors.mag_x + sensorBiases.mag_x + random.gauss(0, sensorSigmas.mag_x)
        sensorNoise.mag_y = trueSensors.mag_y + sensorBiases.mag_y + random.gauss(0, sensorSigmas.mag_y)
        sensorNoise.mag_z = trueSensors.mag_z + sensorBiases.mag_z + random.gauss(0, sensorSigmas.mag_z)

        # Update Pressure Sensor
        sensorNoise.baro = trueSensors.baro + sensorBiases.baro
        sensorNoise.pitot = trueSensors.pitot + sensorBiases.pitot

        # Check to see if the GPS needs updating
        if (self.updateTicks % self.gpsTickUpdate):
            nNoise, eNoise, altNoise = self.gyroGM.update()
            sensorNoise.gps_n = trueSensors.gps_n + nNoise + random.gauss(0, sensorSigmas.gps_n)
            sensorNoise.gps_e = trueSensors.gps_e + eNoise + random.gauss(0, sensorSigmas.gps_e)
            sensorNoise.gps_alt = trueSensors.gps_alt + altNoise + random.gauss(0, sensorSigmas.gps_alt)
            sensorNoise.gps_sog = trueSensors.gps_sog + random.gauss(0, sensorSigmas.gps_sog)

            if math.isclose(sensorNoise.gps_sog, 0.0):
                sensorNoise.gps_cog = trueSensors.gps_cog + random.gauss(0, sensorSigmas.gps_cog * 100)
            else:
                sensorNoise.gps_cog = trueSensors.gps_cog + random.gauss(0, (sensorSigmas.gps_cog * VPC.InitialSpeed) / trueSensors.gps_sog)
            sensorNoise.gps_cog = math.fmod(sensorNoise.gps_cog, math.pi)

        else:
            sensorNoise.gps_n = noisySensors.gps_n
            sensorNoise.gps_e = noisySensors.gps_e
            sensorNoise.gps_alt = noisySensors.gps_alt
            sensorNoise.gps_cog = noisySensors.gps_cog
            sensorNoise.gps_sog = noisySensors.gps_sog

        return sensorNoise

    # Wrapper function to update the Sensors (both true and noisy) using the state and dot held within the self.AeroModel.
    def update(self):
        self.sensorTrue = self.updateSensorsTrue(self.sensorTrue, self.aeroModel.getVehicleState(), self.aeroModel.getVehicleDerivative())
        self.sensorNoisy = self.updateSensorsNoisy(self.sensorTrue, self.sensorNoisy, self.sensorBias, self.sensorSigma)
        self.updateTicks += 1
        return

    # Function to reset the module to run again. Should reset the Gauss-Markov models, re-initialize the sensor biases, and reset the sensors true and noisy to pristine conditions
    def reset(self):
        self.sensorTrue = Sensors.vehicleSensors()
        self.sensorNoisy = Sensors.vehicleSensors()
        self.sensorBias = self.initializeBiases()
        self.sensorSigma = self.initializeSigmas()

        self.gyroGM.reset()
        self.gpsGM.reset()

        self.updateTicks = 0
        return



