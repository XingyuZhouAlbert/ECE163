"""
.. modle:: VehicleAeroDynamicsModel_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from VehicleDynamicsModelGeneration.py with 
	students'  function implementations
.. moduleauthor:: Pavlo Vlastos <pvlastos@ucsc.edu>
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Modeling.VehicleAerodynamicsModel as VAM
import ece163.Containers.States as States
import ece163.Containers.Inputs as Inputs
import ece163.Modeling.WindModel as WindModel
import argparse
import random
import pickle
import traceback
import copy

###############################################################################
inFailureMode = False
inContinueMode = False 

###############################################################################
# Private methods
def printTestBlockResult(function, testsPassed, testCount):
	if testsPassed != testCount:
		addendum = " (TESTS FAILED)"
	else:
		addendum = ""
	print("\t{}/{} tests passed for {}{}()".format(testsPassed, testCount, function.__name__, addendum))

def printTestFailure(function, inputs, outputs, expectedoutputs):
	print("Test Failed for {}. Please find repr version of the inputs below for testing".format(function.__name__))
	print("Inputs: {}".format(repr(inputs)))
	print("Outputs: {}".format(repr(outputs)))
	print("Expected Outputs: {}".format(repr(expectedoutputs)))

def matrixCompare(A, B):
	"""
	Compare the elements of two matrices

	:param A: matrix (list of lists) of [m x n]
	:param B: matrix (list of lists) of [n x r]
	:return: [True, expTot, resTot]: True or False if the matrices match, the
	number of expected matching elements, and the number of resulting matching
	elements between A and B
	"""
	[m, r] = MatrixMath.matrixSize(A)
	[m_c, r_c] = MatrixMath.matrixSize(B)

	expTot = (m_c * r_c)
	resTot = 0

	if (m == m_c) and (r == r_c):
		for row in range(m):
			for col in range(r):
				if math.isclose(A[row][col], B[row][col]) is not True:
					if inFailureMode and (not inContinueMode):
						print("Element [{0},{1}] is incorrect".format(row, col))
					return [False, expTot, resTot]
				else:
					resTot += 1
		if expTot != resTot:
			if inFailureMode and (not inContinueMode):
				print("\r\nResulting matrix dimensions match the expected matrix dimensions")
			return [True, expTot, resTot]
	else:
		if inFailureMode and (not inContinueMode):
			print("Error: Resulting matrix dimensions do not match the expected matrix dimensions")
		return [False, expTot, resTot]

	return [True, expTot, resTot]

def stateCompare(resultingState, expectedState):
	if resultingState != expectedState:
		print("\tState differences:")
		print("\t{}".format(repr(resultingState)))

		# Positions
		resultingElement = resultingState.pn
		expectedElement = expectedState.pn
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Position [m]", "pn", resultingElement, expectedElement))

		resultingElement = resultingState.pe
		expectedElement = expectedState.pe
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Position [m]", "pe", resultingElement, expectedElement))

		resultingElement = resultingState.pd
		expectedElement = expectedState.pd
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Position [m]", "pd", resultingElement, expectedElement))

		# Velocities
		resultingElement = resultingState.u
		expectedElement = expectedState.u
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Velocities [m/s]", "u", resultingElement, expectedElement))

		resultingElement = resultingState.v
		expectedElement = expectedState.v
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Velocities [m/s]", "v", resultingElement, expectedElement))

		resultingElement = resultingState.w
		expectedElement = expectedState.w
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Velocities [m/s]", "w", resultingElement, expectedElement))

		# Euler Angles
		resultingElement = resultingState.yaw
		expectedElement = expectedState.yaw
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Euler Angle [rad]", "yaw", resultingElement, expectedElement))

		resultingElement = resultingState.pitch
		expectedElement = expectedState.pitch
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Euler Angle [rad]", "pitch", resultingElement, expectedElement))

		resultingElement = resultingState.roll
		expectedElement = expectedState.roll
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Euler Angle [rad]", "roll", resultingElement, expectedElement))

		resultingElement = resultingState.R
		expectedElement = expectedState.R
		[matsEqual, expTot, resTot] = matrixCompare(resultingElement, expectedElement)
		if matsEqual is False:
			print("\tDirection Cosine Matrix (DCM) [3x3]:")
			print("\tResulting DCM:")
			MatrixMath.matrixPrint(resultingState.R)
			print("\tExpected DCM:")
			MatrixMath.matrixPrint(expectedState.R)

		resultingElement = resultingState.p
		expectedElement = expectedState.p
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Body Rates [rad/s]", "p", resultingElement, expectedElement))

		resultingElement = resultingState.q
		expectedElement = expectedState.q
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Body Rates [rad/s]", "q", resultingElement, expectedElement))

		resultingElement = resultingState.r
		expectedElement = expectedState.r
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Body Rates [rad/s]", "r", resultingElement, expectedElement))

		resultingElement = resultingState.Va
		expectedElement = expectedState.Va
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Airspeed [m/s]", "Va", resultingElement, expectedElement))

		resultingElement = resultingState.alpha
		expectedElement = expectedState.alpha
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Angle of attack [rad]", "alpha", resultingElement, expectedElement))

		resultingElement = resultingState.beta
		expectedElement = expectedState.beta
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Sideslip Angle", "beta", resultingElement, expectedElement))

		resultingElement = resultingState.chi
		expectedElement = expectedState.chi
		if (resultingElement != expectedElement):
			print("\t{}: {}, result: {}, expected: {}".format("Angle", "chi", resultingElement, expectedElement))

	return

###############################################################################
# Instantiate a vehicle dynamics model object to reference module method names
vehicle = VAM.VehicleAerodynamicsModel()

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all tests regardless of failures')
parser.add_argument('picklePath', nargs='?', default='VehicleAerodynamicsModel_TestData.pickle', help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()

picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Vehicle Aerodynamics Model using file {}".format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed

testBlockIterator = iter(allTests)  # we hard code the tests as well so we need an iterator

###############################################################################
# getVehicleState() and setVehicleState()
print("Comparing outputs for {}() and {}()".format(vehicle.getVehicleState.__name__, vehicle.setVehicleState.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, expectedState in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState = States.vehicleState()

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.getVehicleState, state, resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.getVehicleState, testsPassed, len(curTestBlock))
printTestBlockResult(vehicle.setVehicleState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# gravityForces()
print("Comparing outputs for {}()".format(vehicle.gravityForces.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, expectedOutput in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.gravityForces(state)

		if inFailureMode:
			resultingOutput.Fx += epsilon

		if resultingOutput == expectedOutput:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.gravityForces, state, resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.gravityForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CalculateCoeff_alpha()
print("Comparing outputs for {}()".format(vehicle.CalculateCoeff_alpha.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, expectedOutput in curTestBlock:  # and now we can iterate through them
	expectedOutput = list(expectedOutput)
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.CalculateCoeff_alpha(state.alpha)
		resultingOutput = list(resultingOutput)

		if inFailureMode:
			resultingOutput[0] += epsilon

		lenExpected = len(expectedOutput)
		countAlphaMatch = 0

		for i in range(lenExpected):
			if math.isclose(resultingOutput[i], expectedOutput[i]):
				countAlphaMatch += 1

		if countAlphaMatch == lenExpected:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.CalculateCoeff_alpha, state, resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.CalculateCoeff_alpha, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# aeroForces()
print("Comparing outputs for {}()".format(vehicle.aeroForces.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, expectedOutput in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.aeroForces(state)

		if inFailureMode:
			resultingOutput.Fx += epsilon

		if resultingOutput == expectedOutput:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.aeroForces, state, resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.aeroForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# controlForces()
print("Comparing outputs for {}()".format(vehicle.controlForces.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, control, expectedOutput in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.controlForces(state, control)

		if inFailureMode:
			resultingOutput.Fx += epsilon

		if resultingOutput == expectedOutput:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.controlForces, (state, control), resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.controlForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CalculatePropForces()
print("Comparing outputs for {}()".format(vehicle.CalculatePropForces.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, throttle, expectedOutput in curTestBlock:  # and now we can iterate through them
	expectedOutput = list(expectedOutput)
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.CalculatePropForces(state.Va, throttle)

		if inFailureMode:
			resultingOutput = list(resultingOutput)
			resultingOutput[0] += epsilon

		if math.isclose(resultingOutput[0], expectedOutput[0]) and math.isclose(resultingOutput[1], expectedOutput[1]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.CalculatePropForces, (state, throttle), resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.CalculatePropForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# setWindModel() and getWindState()
print("Comparing outputs for {}() and {}()".format(vehicle.setWindModel.__name__, vehicle.getWindState.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, Wn, We, Wd, expectedOutput in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		vehicle.setWindModel(Wn, We, Wd)
		resultingOutput = vehicle.getWindState()

		if inFailureMode:
			resultingOutput.Wn += epsilon

		if resultingOutput == expectedOutput:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.setWindModel, (state, Wn, We, Wd), resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.setWindModel, testsPassed, len(curTestBlock))
printTestBlockResult(vehicle.getWindState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CalculateAirspeed()
print("Comparing outputs for {}()".format(vehicle.CalculateAirspeed.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, expectedOutput in curTestBlock:  # and now we can iterate through them
	expectedOutput = list(expectedOutput)
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingOutput = vehicle.CalculateAirspeed(state, wind)
		resultingOutput = list(resultingOutput)

		if inFailureMode:
			resultingOutput[0] += epsilon

		lenExpected = len(expectedOutput)
		countMatch = 0

		for i in range(lenExpected):
			if math.isclose(resultingOutput[i], expectedOutput[i]):
				countMatch += 1

		if countMatch == lenExpected:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.CalculateAirspeed, (state, wind), resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.CalculateAirspeed, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces()
print("Comparing outputs for {}()".format(vehicle.updateForces.__name__))
updateForcesTestBlockCount = 20
updateForcesBlockCount = 0

###############################################################################
# updateForces() Block 0 - ALl
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 1 - Wn Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 2 - We Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 3 - Wd Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 4 - Wu Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 5 - Wv Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 6 - Ww Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 7 - Throttle Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 8 - Aileron Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 9 - Elevator Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# updateForces() Block 10 - Rudder Only
print("\t{} Test block {}".format(vehicle.updateForces.__name__, updateForcesBlockCount))
updateForcesBlockCount += 1

curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, wind, control, expectedForce in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(state)
		resultingForce = vehicle.updateForces(copy.deepcopy(state), wind, control)

		if inFailureMode:
			resultingForce.Fx += epsilon

		if resultingForce == expectedForce:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.updateForces, (state, wind, control), resultingForce, expectedForce)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.updateForces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update()
print("Comparing outputs for {}()".format(vehicle.Update.__name__))
UpdateTestBlockCount = 0

###############################################################################
## Update() Block 0
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 1 - Test Throttle Only
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 2 - Test Aileron Only
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 3 - Test Elevator Only
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 4 - Test Rudder Only
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 5 - Test Zero Control Input
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1


###############################################################################
# Update() Block 6 - Walk out
print("\t{} Test block {}".format(vehicle.Update.__name__, UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for UpdateTestCount, previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)

		for i in range(UpdateTestCount):
			vehicle.Update(control)
			
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.Update, (UpdateTestCount, previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# reset()
print("Comparing outputs for {}()".format(vehicle.reset.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for previousState, control, expectedState in curTestBlock:  # and now we can iterate through them
	try:

		# Instantiate a vehicle dynamics model object
		vehicle = VAM.VehicleAerodynamicsModel()

		vehicle.setVehicleState(previousState)
		vehicle.Update(control)
		vehicle.reset()
		resultingState = vehicle.getVehicleState()

		if inFailureMode:
			resultingState.pn += epsilon

		stateCompare(resultingState, expectedState)

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehicle.reset, (previousState, control), resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehicle.reset, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Vehicle Aerodynamics Model")
else:
	print("{}/{} tests blocks passed for Vehicle Aerodynamics Model".format(testBlocksPassed, len(allTests)))

