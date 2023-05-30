import math
import sys
import pickle
import enum
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Controls as Controls
import ece163.Constants.VehiclePhysicalConstants as VPC
import ece163.Modeling.VehicleAerodynamicsModel as VehicleAerodynamicsModule

# Disclaimer: The functions being used in are either from the lecture, the book or the prop cheatsheet.
# I also discuss the general coding concept with Jimmy Chen, Leonid Shuster as well as
# Christian Sabile. The coding are done independently.

class PDControl():
    def __init__(self, kp=0.0, kd=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return


    def setPDGains(self, kp=0.0, kd=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        # Find out the error of command and current measurement
        error = command - current
        # compute the output for PD controller
        u = (self.kp * error) - (self.kd * derivative) + self.trim
        if u < self.lowLimit:
            u = self.lowLimit
        elif u > self.highLimit:
            u = self.highLimit
        return u


class PIControl():
    def __init__(self, dT=VPC.dT, kp=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        self.accumulator = 0.0
        self.prevError = 0.0
        return

    def setPIGains(self, dT=VPC.dT, kp=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def resetIntegrator(self):
        self.accumulator = 0.0
        self.prevError = 0.0
        return

    def Update(self, command=0.0, current=0.0):
        # calculate error
        error = command - current
        self.accumulator += (1/2) * self.dT * (error + self.prevError)
        # compute the output for PI controller
        u = self.trim + (self.kp * error) + (self.ki * self.accumulator)
        # if the output is out of bound, the error should stop accumulating
        if u > self.highLimit:
            u = self.highLimit
            self.accumulator -= (1/2) * self.dT * (error + self.prevError)
        elif u < self.lowLimit:
            u = self.lowLimit
            self.accumulator -= (1/2) * self.dT * (error + self.prevError)

        # Update previous error
        self.prevError = error
        return u

class PIDControl():
    def __init__(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        self.accumulator = 0.0
        self.prevError = 0.0
        return

    def setPIDGains(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def resetIntegrator(self):
        self.accumulator = 0.0
        self.prevError = 0.0
        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        error = command - current
        self.accumulator += (1/2) * self.dT * (error + self.prevError)
        #compute the output for PID controller
        u = self.trim + (self.kp * error) - (self.kd * derivative) + (self.ki * self.accumulator)
        # if the output is out of bound, the error should stop accumulating
        if u > self.highLimit:
            u = self.highLimit
            self.accumulator -= (1/2) * self.dT * (error + self.prevError)
        elif u < self.lowLimit:
            u = self.lowLimit
            self.accumulator -= (1 / 2) * self.dT * (error + self.prevError)
        # Update previous error
        self.prevError = error

        return u

class VehicleClosedLoopControl():
    def __init__(self, dT = VPC.dT):
        self.controlGain = Controls.controlGains()
        self.trimInputs = Inputs.controlInputs()
        self.consurOutputs = Inputs.controlInputs()
        self.VAM = VehicleAerodynamicsModule.VehicleAerodynamicsModel()
        self.VAM.vehicleDynamics.dT = dT
        self.dT = self.VAM.vehicleDynamics.dT
        self.aileronFromRoll = PIDControl()
        self.rollFromCourse = PIControl()
        self.rudderFromSideslip = PIControl()
        self.elevatorFromPitch = PDControl()
        self.throttleFromAirspeed = PIControl()
        self.pitchFromAltitude = PIControl()
        self.pitchFromAirspeed = PIControl()
        self.climbMode = Controls.AltitudeStates.HOLDING
        return

    # Wrapper function to extract control gains from the class.
    def getControlGains(self):
        return self.controlGain
    # Wrapper function to extract the internal VehicleAerodynamicsModel in order to access the various function that are associated with the Aero model
    def getVehicleAerodynamicsModel(self):
        return self.VAM
    # Wrapper function to extract control outputs (Throttle, Aileron, Elevator, Rudder) from the class.
    def getVehicleControlSurfaces(self):
        return self.consurOutputs
    # Wrapper function to extract vehicle state from the class.
    def getVehicleState(self):
        return self.VAM.vehicleDynamics.state

    def getTrimInputs(self):
        return self.trimInputs

    # Function to set all of the gains from the controlGains previously computed to the correct places within the various control loops to do the successive loop closure
    def setControlGains(self, controlGains = Controls.controlGains()):
        self.controlGain = controlGains
        self.aileronFromRoll.setPIDGains(self.dT, self.controlGain.kp_roll, self.controlGain.kd_roll, self.controlGain.ki_roll, self.trimInputs.Aileron, VPC.minControls.Aileron, VPC.maxControls.Aileron)
        # roll has something to do with bank angle
        self.rollFromCourse.setPIGains(self.dT, self.controlGain.kp_course, self.controlGain.ki_course, 0.0, -math.radians(VPC.bankAngleLimit), math.radians(VPC.bankAngleLimit))
        self.rudderFromSideslip.setPIGains(self.dT, self.controlGain.kp_sideslip, self.controlGain.ki_sideslip, self.trimInputs.Rudder, VPC.minControls.Rudder, VPC.maxControls.Rudder)
        self.elevatorFromPitch.setPDGains(self.controlGain.kp_pitch, self.controlGain.kd_pitch, self.trimInputs.Elevator, VPC.minControls.Elevator, VPC.maxControls.Elevator)
        self.throttleFromAirspeed.setPIGains(self.dT, self.controlGain.kp_SpeedfromThrottle, self.controlGain.ki_SpeedfromThrottle, self.trimInputs.Throttle, VPC.minControls.Throttle, VPC.maxControls.Throttle)
        # pitch has something to do with pitchingAngle
        self.pitchFromAltitude.setPIGains(self.dT, self.controlGain.kp_altitude, self.controlGain.ki_altitude, 0.0, -math.radians(VPC.pitchAngleLimit), math.radians(VPC.pitchAngleLimit))
        self.pitchFromAirspeed.setPIGains(self.dT, self.controlGain.kp_SpeedfromElevator, self.controlGain.ki_SpeedfromElevator, 0.0, -math.radians(VPC.pitchAngleLimit), math.radians(VPC.pitchAngleLimit))
        return

    # Wrapper function to inject vehicle state into the class.
    def setVehicleState(self,state):
        self.VAM.vehicleDynamics.state = state
        return

    # Wrapper function to inject the trim inputs into the class.
    def setTrimInputs(self, trimInputs = Inputs.controlInputs()):
        self.trimInputs = trimInputs
        return


    # Resets the module to run again. Does not overwrite control gains, but does reset the integral states of all of the PI control loops.
    def reset(self):
        self.VAM.reset()
        self.aileronFromRoll.resetIntegrator()
        self.rollFromCourse.resetIntegrator()
        self.rudderFromSideslip.resetIntegrator()
        self.throttleFromAirspeed.resetIntegrator()
        self.pitchFromAltitude.resetIntegrator()
        self.pitchFromAirspeed.resetIntegrator()
        return

    # Function that implements the full closed loop controls using the commanded inputs of airspeed, altitude, and course (chi).
    def Update(self, referenceCommands=Controls.referenceCommands()):
        chi = self.VAM.vehicleDynamics.state.chi
        roll = self.VAM.vehicleDynamics.state.roll
        p = self.VAM.vehicleDynamics.state.p
        beta = self.VAM.vehicleDynamics.state.beta
        pd = self.VAM.vehicleDynamics.state.pd
        pitch = self.VAM.vehicleDynamics.state.pitch
        q = self.VAM.vehicleDynamics.state.q
        altitude = -pd


        # define the course error
        courseError = referenceCommands.commandedCourse - chi
        # When the error is outside of +- pi, subtract or add 2pi from my chi
        if courseError >= math.pi:
            chi += (2 * math.pi)
        elif courseError <= -math.pi:
            chi += -(2 * math.pi)

        # Updating control commands, such as rudder and aileron
        rollCmd = self.rollFromCourse.Update(referenceCommands.commandedCourse, chi)
        aileronCmd = self.aileronFromRoll.Update(rollCmd, roll, p)
        rudderCmd = self.rudderFromSideslip.Update(0.0, beta)

        # Minimise controller effort so that it doesnt react to every small changes in altitude
        # Longitudinal controller
        if altitude > (referenceCommands.commandedAltitude + VPC.altitudeHoldZone):
            # Checking for climbMode changes
            if self.climbMode is not Controls.AltitudeStates.DESCENDING:
                self.climbMode = Controls.AltitudeStates.DESCENDING
                # Get rid of the potentially large error
                self.pitchFromAirspeed.resetIntegrator()
            # Reset error
            throttleCmd = VPC.minControls.Throttle
            pitchCmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)

        elif altitude < (referenceCommands.commandedAltitude - VPC.altitudeHoldZone):
            if self.climbMode is not Controls.AltitudeStates.CLIMBING:
                self.climbMode = Controls.AltitudeStates.CLIMBING
                self.pitchFromAirspeed.resetIntegrator()
            throttleCmd = VPC.maxControls.Throttle
            pitchCmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)
        # dealing with the holding zone
        else:
            if self.climbMode is not Controls.AltitudeStates.HOLDING:
                self.climbMode = Controls.AltitudeStates.HOLDING
                self.pitchFromAirspeed.resetIntegrator()
            throttleCmd = self.throttleFromAirspeed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)
            pitchCmd = self.pitchFromAltitude.Update(referenceCommands.commandedAltitude, altitude)

        # Dealing with the elevator
        elevatorCmd = self.elevatorFromPitch.Update(pitchCmd, pitch, q)

        # Feeding commands back to controlsurface outputs
        self.consurOutputs.Throttle = throttleCmd
        self.consurOutputs.Rudder = rudderCmd
        self.consurOutputs.Elevator = elevatorCmd
        self.consurOutputs.Aileron = aileronCmd
        referenceCommands.commandedRoll = rollCmd
        referenceCommands.commandedPitch = pitchCmd
        self.VAM.Update(self.consurOutputs)
        return
















