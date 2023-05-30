"""
.. module:: ClosedLoopControl_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from ClosedLoopControl_Generation.py with 
	students' function implementations
.. moduleauthor:: Pavlo Vlastos <pvlastos@ucsc.edu>
"""

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Modeling.VehicleAerodynamicsModel as VAM
import ece163.Modeling.WindModel as WM
import ece163.Containers.States as States
import ece163.Containers.Inputs as Inputs
import ece163.Containers.Controls as Controls
import ece163.Controls.VehicleControlGains as CG
import ece163.Controls.VehicleClosedLoopControl as CLC
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
	print("\t{}/{} tests passed for {}{}()".format(testsPassed, testCount, \
		function.__name__, addendum))

def printTestFailure(function, inputs, outputs, expectedoutputs):
	print("Test Failed for {}. Please find repr version of the inputs below for\
	 testing".format(function.__name__))
	print("Inputs: {}".format(repr(inputs)))
	print("Outputs: {}".format(repr(outputs)))
	print("Expected Outputs: {}".format(repr(expectedoutputs)))

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all\
 tests regardless of failures')
parser.add_argument('picklePath', nargs='?', \
	default='ClosedLoopControl_TestData.pickle', \
	help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail\
 all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()

picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Closed Loop Control using file {}".\
	format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed
testBlockIterator = iter(allTests)  # we hard code the tests as well so we need 
# an iterator

###############################################################################
# CLC.PDControl.setPDGains()
print("Comparing outputs for {}()".format(CLC.PDControl.setPDGains.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for kp, kd, trim, lowLimit, highLimit, expectedPD in curTestBlock:
	try:
		PD = CLC.PDControl()
		PD.setPDGains(kp, kd, trim, lowLimit, highLimit)

		if inFailureMode:
			PD.kp += epsilon

		eqFlagStr = list()
		eqFlagStr.append(math.isclose(PD.kp, expectedPD.kp))
		eqFlagStr.append(math.isclose(PD.kd, expectedPD.kd))
		eqFlagStr.append(math.isclose(PD.trim, expectedPD.trim))
		eqFlagStr.append(math.isclose(PD.lowLimit, expectedPD.lowLimit))
		eqFlagStr.append(math.isclose(PD.highLimit, expectedPD.highLimit))

		if all(elems == True for elems in eqFlagStr):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PDControl.setPDGains, (kp, kd, trim, \
					lowLimit, highLimit), PD, expectedPD)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PDControl.setPDGains, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PDControl.Update()
print("Comparing outputs for {}()".format(CLC.PDControl.Update.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for kp, kd, trim, lowLimit, highLimit, command, current, derivative, \
expectedOutput in curTestBlock:
	try:
		PD = CLC.PDControl()
		PD.setPDGains(kp, kd, trim, lowLimit, highLimit)
		result = PD.Update(command, current, derivative)

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expectedOutput):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PDControl.Update, (kp, kd, trim, \
					lowLimit, highLimit, command, current, derivative), result,\
				expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PDControl.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIControl.setPIGains()
print("Comparing outputs for {}()".format(CLC.PIControl.setPIGains.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dt, kp, ki, trim, lowLimit, highLimit, expectedPI in curTestBlock:
	try:
		PI = CLC.PIControl()
		PI.setPIGains(dt, kp, ki, trim, lowLimit, highLimit)

		if inFailureMode:
			result += epsilon

		if inFailureMode:
			PI.kp += epsilon

		eqFlagStr = list()
		eqFlagStr.append(math.isclose(PI.dT, expectedPI.dT))
		eqFlagStr.append(math.isclose(PI.kp, expectedPI.kp))
		eqFlagStr.append(math.isclose(PI.ki, expectedPI.ki))
		eqFlagStr.append(math.isclose(PI.accumulator, expectedPI.accumulator))
		eqFlagStr.append(math.isclose(PI.prevError, expectedPI.prevError))
		eqFlagStr.append(math.isclose(PI.trim, expectedPI.trim))
		eqFlagStr.append(math.isclose(PI.lowLimit, expectedPI.lowLimit))
		eqFlagStr.append(math.isclose(PI.highLimit, expectedPI.highLimit))

		if all(elems == True for elems in eqFlagStr):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIControl.setPIGains, (dt, kp, ki, trim, \
					lowLimit, highLimit), PI, expectedPI)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIControl.setPIGains, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIControl.Update()
print("Comparing outputs for {}()".format(CLC.PIControl.Update.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dt, kp, ki, trim, lowLimit, highLimit, command, current, expectedOutput in \
curTestBlock:
	try:
		PI = CLC.PIControl()
		PI.setPIGains(dt, kp, ki, trim, lowLimit, highLimit)
		result = PI.Update(command, current)

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expectedOutput):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIControl.Update, (dt, kp, ki, trim, \
					lowLimit, highLimit, command, current), result,\
				expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIControl.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIControl.resetIntegrator()
print("Comparing outputs for {}()".format(CLC.PIControl.resetIntegrator.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for curTestNumber, dt, kp, ki, trim, lowLimit, highLimit, command, current, expectedOutput in \
curTestBlock:
	try:
		PI = CLC.PIControl()
		PI.setPIGains(dt, kp, ki, trim, lowLimit, highLimit)

		for i in range(curTestNumber):
			PI.Update(command, current)

		PI.resetIntegrator()

		if inFailureMode:
			PI.accumulator += epsilon
			PI.prevError += epsilon

		if math.isclose(PI.accumulator, expectedOutput.accumulator) and \
		math.isclose(PI.prevError, expectedOutput.prevError):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIControl.resetIntegrator,(curTestNumber,\
					dt, kp, ki, trim, lowLimit, highLimit, command, current),\
				result, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIControl.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIDControl.setPIDGains()
print("Comparing outputs for {}()".format(CLC.PIDControl.setPIDGains.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dt, kp, kd, ki, trim, lowLimit, highLimit, expectedPID in curTestBlock:
	try:
		PID = CLC.PIDControl()
		PID.setPIDGains(dt, kp, kd, ki, trim, lowLimit, highLimit)

		if inFailureMode:
			PID.kp += epsilon

		eqFlagStr = list()
		eqFlagStr.append(math.isclose(PID.dT, expectedPID.dT))
		eqFlagStr.append(math.isclose(PID.kp, expectedPID.kp))
		eqFlagStr.append(math.isclose(PID.ki, expectedPID.ki))
		eqFlagStr.append(math.isclose(PID.kd, expectedPID.kd))
		eqFlagStr.append(math.isclose(PID.trim, expectedPID.trim))
		eqFlagStr.append(math.isclose(PID.lowLimit, expectedPID.lowLimit))
		eqFlagStr.append(math.isclose(PID.highLimit, expectedPID.highLimit))

		if all(elems == True for elems in eqFlagStr):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIDControl.setPIDGains, (dt, kp, kd, ki, trim, \
					lowLimit, highLimit), PID, expectedPID)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIDControl.setPIDGains, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIDControl.Update()
print("Comparing outputs for {}()".format(CLC.PIDControl.Update.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dt, kp, kd, ki, trim, lowLimit, highLimit, command, current, derivative, \
expectedOutput in curTestBlock:
	try:
		PID = CLC.PIDControl()
		PID.setPIDGains(dt, kp, kd, ki, trim, lowLimit, highLimit)
		result = PID.Update(command, current, derivative)

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expectedOutput):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIDControl.Update, (dt, kp, ki, trim, \
					lowLimit, highLimit, command, current), result,\
				expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIDControl.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.PIDControl.resetIntegrator()
print("Comparing outputs for {}()".format(CLC.PIDControl.resetIntegrator.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for curTestNumber, dt, kp, kd, ki, trim, lowLimit, highLimit, command, current, expectedOutput in \
curTestBlock:
	try:
		PID = CLC.PIDControl()
		PID.setPIDGains(dt, kp, kd, ki, trim, lowLimit, highLimit)

		for i in range(curTestNumber):
			PID.Update(command, current)

		PID.resetIntegrator()

		if inFailureMode:
			PID.accumulator += epsilon
			PID.prevError += epsilon

		if math.isclose(PID.accumulator, expectedOutput.accumulator) and \
		math.isclose(PID.prevError, expectedOutput.prevError):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.PIDControl.resetIntegrator,(curTestNumber,\
					dt, kp, ki, trim, lowLimit, highLimit, command, current),\
				result, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.PIDControl.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.VehicleClosedLoopControl.setControlGains()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.setControlGains.__name__))
# CLC.VehicleClosedLoopControl.getControlGains()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.getControlGains.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for cg, expectedCG in curTestBlock:
	try:
		VCLC = CLC.VehicleClosedLoopControl()
		VCLC.setControlGains(cg)
		resultCG = VCLC.getControlGains()

		if inFailureMode:
			resultCG.kp_roll += epsilon

		if math.isclose(resultCG.kp_roll, expectedCG.kp_roll) and \
		math.isclose(resultCG.kd_roll, expectedCG.kd_roll) and \
		math.isclose(resultCG.ki_roll, expectedCG.ki_roll) and \
		math.isclose(resultCG.kp_sideslip, expectedCG.kp_sideslip) and \
		math.isclose(resultCG.ki_sideslip, expectedCG.ki_sideslip) and \
		math.isclose(resultCG.kp_course, expectedCG.kp_course) and \
		math.isclose(resultCG.ki_course, expectedCG.ki_course) and \
		math.isclose(resultCG.kp_pitch, expectedCG.kp_pitch) and \
		math.isclose(resultCG.kd_pitch, expectedCG.kd_pitch) and \
		math.isclose(resultCG.kp_altitude, expectedCG.kp_altitude) and \
		math.isclose(resultCG.ki_altitude, expectedCG.ki_altitude) and \
		math.isclose(resultCG.kp_SpeedfromThrottle, expectedCG.kp_SpeedfromThrottle) and \
		math.isclose(resultCG.ki_SpeedfromThrottle, expectedCG.ki_SpeedfromThrottle) and \
		math.isclose(resultCG.kp_SpeedfromElevator, expectedCG.kp_SpeedfromElevator) and \
		math.isclose(resultCG.ki_SpeedfromElevator, expectedCG.ki_SpeedfromElevator):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.VehicleClosedLoopControl.setControlGains,(cg),\
				resultCG, expectedCG)
				printTestFailure(CLC.VehicleClosedLoopControl.getControlGains,(cg),\
				resultCG, expectedCG)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.VehicleClosedLoopControl.setControlGains, testsPassed, len(curTestBlock))
printTestBlockResult(CLC.VehicleClosedLoopControl.getControlGains, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.VehicleClosedLoopControl.setTrimInputs()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.setTrimInputs.__name__))
# CLC.VehicleClosedLoopControl.getControlGains()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.getTrimInputs.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for ct, expectedOutput in curTestBlock:
	try:
		VCLC = CLC.VehicleClosedLoopControl()
		VCLC.setTrimInputs(ct)
		result = VCLC.getTrimInputs()
		
		if inFailureMode:
			result = Inputs.controlInputs()

		if result == expectedOutput:
				testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.VehicleClosedLoopControl.setTrimInputs,(ct),\
				result, expectedOutput)
				printTestFailure(CLC.VehicleClosedLoopControl.getTrimInputs,(ct),\
				result, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.VehicleClosedLoopControl.setTrimInputs, testsPassed, len(curTestBlock))
printTestBlockResult(CLC.VehicleClosedLoopControl.getTrimInputs, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.VehicleClosedLoopControl.setVehicleState()
# CLC.VehicleClosedLoopControl.getVehicleState()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.setVehicleState.__name__))
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.getVehicleState.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, expectedState in curTestBlock:
	try:
		VCLC = CLC.VehicleClosedLoopControl()
		VCLC.setVehicleState(state)
		resultingState = VCLC.getVehicleState()
		
		if inFailureMode:
			resultingState.pn += epsilon

		if resultingState == expectedState:
				testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.VehicleClosedLoopControl.setVehicleState,\
					state, resultingState, expectedState)
				printTestFailure(CLC.VehicleClosedLoopControl.getVehicleState,\
					state, resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.VehicleClosedLoopControl.setVehicleState, testsPassed, len(curTestBlock))
printTestBlockResult(CLC.VehicleClosedLoopControl.getVehicleState, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.VehicleClosedLoopControl.Update()
# CLC.VehicleClosedLoopControl.getVehicleControlSurfaces()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.Update.__name__))
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.getVehicleControlSurfaces.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for curSubtestNum, cg, ct, ref, expSurf, expState, expVAM in curTestBlock:
	try:
		VCLC = CLC.VehicleClosedLoopControl()
		VCLC.setControlGains(cg)
		VCLC.setTrimInputs(ct)

		for i in range(curSubtestNum):
			VCLC.Update(ref)

		resSurf = VCLC.getVehicleControlSurfaces()
		resState = VCLC.getVehicleState()
		resVAM = VCLC.getVehicleAerodynamicsModel()
		
		if inFailureMode:
			resSurf.Throttle += epsilon

		if (resSurf == expSurf) and (resState == expState) and \
		(resVAM.getVehicleState() == expVAM.getVehicleState()):
				testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.VehicleClosedLoopControl.Update,\
					(curSubtestNum, cg, ct, ref), (resSurf, resState), (expSurf, expState))
				printTestFailure(CLC.VehicleClosedLoopControl.getVehicleControlSurfaces,\
					(curSubtestNum, cg, ct, ref), (resSurf, resState), (expSurf, expState))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.VehicleClosedLoopControl.Update, testsPassed, len(curTestBlock))
printTestBlockResult(CLC.VehicleClosedLoopControl.getVehicleControlSurfaces, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CLC.VehicleClosedLoopControl.reset()
print("Comparing outputs for {}()".\
	format(CLC.VehicleClosedLoopControl.Update.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for curSubtestNum, cg, ct, ref, expSurf, expState, expVAM in curTestBlock:
	try:
		VCLC = CLC.VehicleClosedLoopControl()
		VCLC.setControlGains(cg)
		VCLC.setTrimInputs(ct)

		for i in range(curSubtestNum):
			VCLC.Update(ref)

		VCLC.reset()

		resSurf = VCLC.getVehicleControlSurfaces()
		resState = VCLC.getVehicleState()
		resVAM = VCLC.getVehicleAerodynamicsModel()
		
		if inFailureMode:
			resSurf.Throttle += epsilon

		if (resSurf == expSurf) and (resState == expState) and \
		(resVAM.getVehicleState() == expVAM.getVehicleState()):
				testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CLC.VehicleClosedLoopControl.reset,\
					(curSubtestNum, cg, ct, ref), (resSurf, resState), (expSurf, expState))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CLC.VehicleClosedLoopControl.reset, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1
###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Closed Loop Control")
else:
	print("{}/{} tests blocks passed for Closed Loop Control".\
		format(testBlocksPassed, len(allTests)))

