"""
.. modle:: VehicleDynamicsModel_TestHarness.py
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
import ece163.Modeling.VehicleDynamicsModel as VDM
import ece163.Containers.States as States
import argparse
import random
import pickle
import traceback

inFailureMode = False
inContinueMode = False 

#  these two functions allow for more standardized output, they should be copied to each test harness and customized
def printTestBlockResult(function, testsPassed, testCount):
	if testsPassed != testCount:
		addendum = " (TESTS FAILED)"
	else:
		addendum = ""
	print("{}/{} tests passed for {}{}()".format(testsPassed, testCount, function.__name__, addendum))

def printTestFailure(function, inputs, outputs, expectedoutputs):
	print("Test Failed for {}. Please find repr version of the inputs below for testing".format(function.__name__))
	print("Inputs: {}".format(repr(inputs)))
	print("Outputs: {}".format(repr(outputs)))
	print("Expected Outputs: {}".format(repr(expectedoutputs)))

# Private methods
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
###############################################################################
# Instantiate a vehicle dynamics model object to reference module method names
vehDyMo = VDM.VehicleDynamicsModel()

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all tests regardless of failures')
parser.add_argument('picklePath', nargs='?', default='VehicleDynamicsModel_TestData.pickle', help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()

picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Vehilce Dynamics Model using file {}".format(picklePath))
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
print("Comparing outputs for {}() and {}()".format(vehDyMo.getVehicleState.__name__, vehDyMo.setVehicleState.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for state, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		vehDyMo.setVehicleState(state)
		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.getVehicleState, state, result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.getVehicleState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# resetVehicleState()
print("Comparing outputs for {}()".format(vehDyMo.resetVehicleState.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for previousState, expectedState in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		vehDyMo.setVehicleState(previousState)
		vehDyMo.resetVehicleState()
		resultState = vehDyMo.getVehicleState()
		
		if inFailureMode:
			resultState.pn += epsilon

		if resultState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.resetVehicleState, previousState, resultState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.resetVehicleState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# derivative()
print("Comparing outputs for {}()".format(vehDyMo.derivative.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for previousState, fmInputs, expectedStateDot in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		resultStateDot = vehDyMo.derivative(previousState, fmInputs)

		if inFailureMode:
			resultStateDot.pn += epsilon

		if resultStateDot == expectedStateDot:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.derivative, (previousState, fmInputs), resultStateDot, expectedStateDot)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.derivative, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# reset()
print("Comparing outputs for {}()".format(vehDyMo.reset.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for previousState, expectedState in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		vehDyMo.setVehicleState(previousState)
		vehDyMo.reset()
		resultState = vehDyMo.getVehicleState()
		
		if inFailureMode:
			resultState.pn += epsilon

		if resultState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.resetVehicleState, previousState, resultState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.resetVehicleState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Rexp()
print("Comparing outputs for {}()".format(vehDyMo.Rexp.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for dT, previousState, dot, expectedResult in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		result = vehDyMo.Rexp(dT, previousState, dot)

		if inFailureMode:
			result[0][0] += epsilon

		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Rexp, (dT, previousState, dot), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Rexp, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# IntegrateState()
print("Comparing outputs for {}()".format(vehDyMo.IntegrateState.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for dT, previousState, dot, expectedResult in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		result = vehDyMo.IntegrateState(dT, previousState, dot)

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.IntegrateState, (dT, previousState, dot), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.IntegrateState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# ForwardEuler()
print("Comparing outputs for {}()".format(vehDyMo.ForwardEuler.__name__))
curTestBlock = next(testBlockIterator)  
testsPassed = 0
for previousState, fmInputs, expectedResult in curTestBlock: 
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel()

		vehDyMo.setVehicleState(previousState)
		result = vehDyMo.ForwardEuler(fmInputs)

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.ForwardEuler, (previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.ForwardEuler, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update()
print("Comparing outputs for {}()".format(vehDyMo.Update.__name__))
UpdateTestBlockCount = 0

###############################################################################
# Update() Block 0
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 1
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 2
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 3
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 4
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 5
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update() Block 6
print("Update Test Block {}".format(UpdateTestBlockCount))
UpdateTestBlockCount += 1
curTestBlock = next(testBlockIterator)
testsPassed = 0
for dt, previousState, fmInputs, expectedResult in curTestBlock:
	try:
		# Instantiate a vehicle dynamics model object
		vehDyMo = VDM.VehicleDynamicsModel(dt)
		vehDyMo.setVehicleState(previousState)

		vehDyMo.Update(fmInputs)

		result = vehDyMo.getVehicleState()

		if inFailureMode:
			result.pn += epsilon

		if result == expectedResult:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(vehDyMo.Update, (dt, previousState, fmInputs), result, expectedResult)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(vehDyMo.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Vehicle Dynamics Model")
else:
	print("{}/{} tests blocks passed for Vehicle Dynamics Model".format(testBlocksPassed, len(allTests)))

