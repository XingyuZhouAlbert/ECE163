"""
File contains the classes for controls primitives, gains, and tuning knobs for tuning controls models.
"""

import math
import enum

from ..Constants import VehiclePhysicalConstants as VPC

testingAbs_tol = 1e-4

class AltitudeStates(enum.Enum):
	"""
	class AltitudeStates(enum.Enum):
	Enumeration class for the altitude hold state machine. Defines three states that we will be using to reset the PI integrators
	when switching between different altitudes.
	"""
	CLIMBING = enum.auto()
	HOLDING = enum.auto()
	DESCENDING = enum.auto()

class referenceCommands():
	def __init__(self, courseCommand=VPC.InitialYawAngle, altitudeCommand=-VPC.InitialDownPosition, airspeedCommand=VPC.InitialSpeed):
		"""
		def __init__(self, courseCommand=VPC.InitialYawAngle, altitudeCommand=-VPC.InitialDownPosition, airspeedCommand=VPC.InitialSpeed):
		Class to hold the commanded inputs for closed loop control.

		:param courseCommand: commanded course over ground [rad], measured + from Inertial North
		:param altitudeCommand: commanded height [m]
		:param airspeedCommand: commanded airspeed [m/s]
		:return: none
		"""
		self.commandedCourse = courseCommand
		self.commandedAltitude = altitudeCommand
		self.commandedAirspeed = airspeedCommand
		self.commandedRoll = 0.0	# These will be set by the control loops internally to the successive loop closure
		self.commandedPitch = 0.0	# These will be set by the control loops internally to the successive loop closure
		return

	def __str__(self):
		return "referenceCommands:(course={}, altitude={}, " \
			   "airspeed={}, roll={}, pitch={}".format(math.degrees(self.commandedCourse), self.commandedAltitude,
													   self.commandedAirspeed, math.degrees(self.commandedRoll),
													   math.degrees(self.commandedPitch))

	def __eq__(self, other):
		if isinstance(other, type(self)):
			if not all([math.isclose(getattr(self, member), getattr(other, member), abs_tol=testingAbs_tol) for member in ['commandedCourse',
																								   'commandedAltitude',
																								   'commandedAirspeed',
																								   'commandedRoll',
																								   'commandedPitch']]):
				return False
			else:
				return True
		else:
			return NotImplemented


class controlGains():
	def __init__(self):
		"""
		Class to hold the control gains for both lateral and longitudinal autopilots in the successive loop closure method
		described in Beard Chapter 6.

		:return: none
		"""
		# Lateral Gains
		self.kp_roll = 0.0
		self.kd_roll = 0.0
		self.ki_roll = 0.0
		self.kp_sideslip = 0.0
		self.ki_sideslip = 0.0
		self.kp_course = 0.0
		self.ki_course = 0.0
		# Longitudinal Gains
		self.kp_pitch = 0.0
		self.kd_pitch = 0.0
		self.kp_altitude = 0.0
		self.ki_altitude = 0.0
		self.kp_SpeedfromThrottle = 0.0
		self.ki_SpeedfromThrottle = 0.0
		self.kp_SpeedfromElevator = 0.0
		self.ki_SpeedfromElevator = 0.0
		return

	def __eq__(self, other):
		if isinstance(other, type(self)):
			if not all([math.isclose(getattr(self, member), getattr(other, member), abs_tol=testingAbs_tol)
						for member in ['kp_roll', 'kd_roll', 'ki_roll', 'kp_sideslip', 'ki_sideslip', 'kp_course',
									   'ki_course', 'kp_pitch', 'kd_pitch', 'kp_altitude', 'ki_altitude',
									   'kp_SpeedfromThrottle', 'ki_SpeedfromThrottle', 'kp_SpeedfromElevator',
									   'ki_SpeedfromElevator']]):
				return False
			else:
				return True
		else:
			return NotImplemented

	def __str__(self):
		return "{0.__name__}(kp_roll={1.kp_roll}, kd_roll={1.kd_roll}, ki_roll={1.ki_roll}, kp_sideslip={1.kp_sideslip}, " \
			   "ki_sideslip={1.ki_sideslip}, kp_course={1.kp_course}, ki_course={1.ki_course}, kp_pitch={1.kp_pitch}, " \
			   "kd_pitch={1.kd_pitch}, kp_altitude={1.kp_altitude}, ki_altitude={1.ki_altitude}, kp_SpeedfromThrottle={1.kp_SpeedfromThrottle}, " \
			   "ki_SpeedfromThrottle={1.ki_SpeedfromThrottle}, kp_SpeedfromElevator={1.kp_SpeedfromElevator}, " \
			   "ki_SpeedfromElevator={1.ki_SpeedfromElevator})".format(type(self), self)


class controlTuning():
	def __init__(self):
		"""
		Class to hold the tuning knobs for both lateral and longitudinal autopilots in the successive loop closure method
		described in Beard Chapter 6. Note that in the successive loop closure methodology, the gains are determined for
		the inner-most loop and the bandwidth separation is used to ensure that the outer loops do not interfere with the
		inner loops. Typical 5-10 ratios at a minimum between natural frequency of the loops.

		:return: none
		"""
		# tuning knobs for lateral control (ignoring Ki_phi)
		self.Wn_roll = 0.0
		self.Zeta_roll = 0.0
		self.Wn_course = 0.0	# Wn_roll should be 5-10x larger
		self.Zeta_course = 0.0
		self.Wn_sideslip = 0.0
		self.Zeta_sideslip = 0.0
		#tuning knows for longitudinal control
		self.Wn_pitch = 0.0
		self.Zeta_pitch = 0.0
		self.Wn_altitude = 0.0	# Wn_pitch should be 5-10x larger
		self.Zeta_altitude = 0.0
		self.Wn_SpeedfromThrottle = 0.0
		self.Zeta_SpeedfromThrottle = 0.0
		self.Wn_SpeedfromElevator = 0.0
		self.Zeta_SpeedfromElevator = 0.0
		return

	def __str__(self):
		return "{0.__name__}(Wn_roll={1.Wn_roll}, Zeta_roll={1.Zeta_roll}, Wn_course={1.Wn_course}, Zeta_course={1.Zeta_course}, " \
			   "Wn_sideslip={1.Wn_sideslip}, Zeta_sideslip={1.Zeta_sideslip}, Wn_pitch={1.Wn_pitch}, Zeta_pitch={1.Zeta_pitch}, " \
			   "Wn_altitude={1.Wn_altitude}, Zeta_altitude={1.Zeta_altitude}, Wn_SpeedfromThrottle={1.Wn_SpeedfromThrottle}, " \
			   "Zeta_SpeedfromThrottle={1.Zeta_SpeedfromThrottle}, Wn_SpeedfromElevator={1.Wn_SpeedfromElevator}, " \
			   "Zeta_SpeedfromElevator={1.Zeta_SpeedfromElevator})".format(type(self), self)



	def __eq__(self, other):
		if isinstance(other, type(self)):
			if not all([math.isclose(getattr(self, member), getattr(other, member), abs_tol=testingAbs_tol)
						for member in ['Wn_roll', 'Zeta_roll', 'Wn_course', 'Zeta_course', 'Wn_sideslip', 'Zeta_sideslip',
									   'Wn_pitch', 'Zeta_pitch', 'Wn_altitude', 'Zeta_altitude', 'Wn_SpeedfromThrottle',
									   'Zeta_SpeedfromThrottle', 'Wn_SpeedfromElevator', 'Zeta_SpeedfromElevator']]):
				return False
			else:
				return True
		else:
			return NotImplemented