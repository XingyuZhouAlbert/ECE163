"""
File contains the class primitive for the collection of sensors available on the UAV.
"""

import math

testingAbs_tol = 1e-6

class vehicleSensors:
	def __init__(self):
		"""
		Defines the typical sensor suite on the UAV. This includes a 3-axis accelerometer, a 3-axis gyro, and a 3-axis
		magnetometer. A simple barometer (absolute pressure sensor) and a pitot tube (differential pressure sensor) are
		also included. Lastly, GPS is included, providing position and both speed over ground (SOG) and course over
		ground (COG).
		"""
		# gyros
		self.gyro_x = 0.0
		self.gyro_y = 0.0
		self.gyro_z = 0.0
		# accelerometers
		self.accel_x = 0.0
		self.accel_y = 0.0
		self.accel_z = 0.0
		# magnetometers
		self.mag_x = 0.0
		self.mag_y = 0.0
		self.mag_z = 0.0
		# pressure sensors
		self.baro = 0.0
		self.pitot = 0.0
		# gps
		self.gps_n = 0.0
		self.gps_e = 0.0
		self.gps_alt = 0.0
		self.gps_sog = 0.0
		self.gps_cog = 0.0
		return

	def __repr__(self):
		return "{0.__name__}(gyro_x={1.gyro_x}, gyro_y={1.gyro_y}, gyro_z={1.gyro_z}, accel_x={1.accel_x}, accel_y={1.accel_y}, " \
			   "accel_z={1.accel_z}, mag_x={1.mag_x}, mag_y={1.mag_y}, mag_z={1.mag_z}, baro={1.baro}, pitot={1.pitot}, " \
			   "gps_n={1.gps_n}, gps_e={1.gps_e}, gps_alt={1.gps_alt}, gps_sog={1.gps_sog}, gps_cog={1.gps_cog})".format(type(self), self)

	def __eq__(self, other):
		if isinstance(other, type(self)):
			if not all(
					[math.isclose(getattr(self, member), getattr(other, member), abs_tol=testingAbs_tol) for member in
					 ['gyro_x', 'gyro_y', 'gyro_z', 'accel_x', 'accel_y', 'accel_z', 'mag_x', 'mag_y', 'mag_z', 'baro',
					  'pitot', 'gps_n', 'gps_e', 'gps_alt', 'gps_sog', 'gps_cog']]):
				return False
			else:
				return True
		else:
			return NotImplemented

