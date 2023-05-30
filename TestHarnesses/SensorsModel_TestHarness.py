"""
.. module:: SensorsModel_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from SensorsModel_Generation.py with 
	students' function implementations
.. moduleauthor:: Pavlo Vlastos <pvlastos@ucsc.edu>
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Sensors.SensorsModel as SM
import ece163.Containers.Sensors as Sensors
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
def printTestBlockResult(currClass, function, testsPassed, testCount):
	if testsPassed != testCount:
		addendum = " (TESTS FAILED)"
	else:
		addendum = ""
	print("\t{}/{} tests passed for {}.{}{}()".format(testsPassed, testCount, \
		currClass.__name__, function.__name__, addendum))

def printTestFailure(currClass, function, inputs, outputs, expectedoutputs):
	print("Test Failed for {}.{}. Please find repr version of the inputs below for\
	 testing".format(currClass.__name__, function.__name__))
	print("Inputs          : {}".format(repr(inputs)))
	print("Outputs         : {}".format(repr(outputs)))
	print("Expected Outputs: {}".format(repr(expectedoutputs)))

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all\
 tests regardless of failures')
parser.add_argument('picklePath', nargs='?', \
	default='SensorsModel_TestData.pickle', \
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

print("Beginning Test Harness for Sensors Model using file {}".\
	format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed
testBlockIterator = iter(allTests)  # we hard code the tests as well so we need 

###############################################################################
# SM.GaussMarkov.update()
currentFuncDef = SM.GaussMarkov.update
currentClass = SM.GaussMarkov
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, dt, tau, eta, vNoiseList, expected in curTestBlock:
	try:
		gm = SM.GaussMarkov(dt, tau, eta)

		for j in range(subTestCount):
			result = gm.update(vNoiseList[j])

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expected):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, dt, tau, eta, \
					vNoiseList), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.GaussMarkov.reset()
currentFuncDef = SM.GaussMarkov.reset
currentClass = SM.GaussMarkov
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dt, tau, eta, expected in curTestBlock:
	try:
		gm = SM.GaussMarkov(dt, tau, eta)
		gm.update()
		gm.reset()
		result = gm.update(0) #this is to get the value, without adding noise, because there is no getter function for self.v
		
		if inFailureMode:
			result += epsilon

		if math.isclose(result, expected):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (dt, tau, eta), result, \
					expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.GaussMarkovXYZ.update()
currentFuncDef = SM.GaussMarkovXYZ.update
currentClass = SM.GaussMarkovXYZ
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, dt, tauX, tauY, tauZ, etaX, etaY, etaZ, vNoiseList, expected \
in curTestBlock:
	try:
		gmXYZ = SM.GaussMarkovXYZ(dt, tauX, tauY, tauZ, etaX, etaY, etaZ)
		
		for j in range(subTestCount):
			vNoiseX, vNoiseY, vNoiseZ = vNoiseList[j]
			result = gmXYZ.update(vNoiseX, vNoiseY, vNoiseZ)

		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, dt, tauX, tauY,\
					tauZ, etaX, etaY, etaZ, vNoiseList), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.GaussMarkovXYZ.reset()
currentFuncDef = SM.GaussMarkovXYZ.reset
currentClass = SM.GaussMarkovXYZ
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, dt, tauX, tauY, tauZ, etaX, etaY, etaZ, vNoiseList, expected \
in curTestBlock:
	try:
		gmXYZ = SM.GaussMarkovXYZ(dt, tauX, tauY, tauZ, etaX, etaY, etaZ)
		
		for j in range(subTestCount):
			vNoiseX, vNoiseY, vNoiseZ = vNoiseList[j]
			gmXYZ.update(vNoiseX, vNoiseY, vNoiseZ)
		gmXYZ.reset()
		result = gmXYZ.update(0, 0, 0)

		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, dt, tauX, tauY,\
					tauZ, etaX, etaY, etaZ, vNoiseList), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.initializeBiases()
currentFuncDef = SM.SensorsModel.initializeBiases
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, gyroBias, accelBias, magBias, baroBias, pitotBias, \
	expectedMaxBiases, expectedMinBiases in curTestBlock:
	try:
		gyro_x_biasList = list()
		gyro_y_biasList = list()
		gyro_z_biasList = list()

		accel_x_biasList = list()
		accel_y_biasList = list()
		accel_z_biasList = list()

		mag_x_biasList = list()
		mag_y_biasList = list()
		mag_z_biasList = list()

		baro_biasList = list()
		pitot_biasList = list()

		for j in range(subTestCount):
			sm = SM.SensorsModel()

			biases = sm.initializeBiases(gyroBias, accelBias, magBias, baroBias, \
				pitotBias)

			gyro_x_biasList.append(biases.gyro_x)
			gyro_y_biasList.append(biases.gyro_y)
			gyro_z_biasList.append(biases.gyro_z)
			
			accel_x_biasList.append(biases.accel_x)
			accel_y_biasList.append(biases.accel_y)
			accel_z_biasList.append(biases.accel_z)

			mag_x_biasList.append(biases.mag_x)
			mag_y_biasList.append(biases.mag_y)
			mag_z_biasList.append(biases.mag_z)

			baro_biasList.append(biases.baro)
			pitot_biasList.append(biases.pitot)

		# We simply check that for a large number of calls to intializeBiases() 
		# that the minimum and maximum bounds are close enough. If a value is 
		# outside the bounds by large amount, then there is definitely 
		# something wrong with intializeBiases()
		resultingMaxBiases = Sensors.vehicleSensors()
		resultingMinBiases = Sensors.vehicleSensors()

		resultingMaxBiases.gyro_x = round(max(gyro_x_biasList), 2)
		resultingMinBiases.gyro_x = round(min(gyro_x_biasList), 2)

		resultingMaxBiases.gyro_y = round(max(gyro_y_biasList), 2)
		resultingMinBiases.gyro_y = round(min(gyro_y_biasList), 2)

		resultingMaxBiases.gyro_z = round(max(gyro_z_biasList), 2)
		resultingMinBiases.gyro_z = round(min(gyro_z_biasList), 2)

		resultingMaxBiases.accel_x = round(max(accel_x_biasList), 2)
		resultingMinBiases.accel_x = round(min(accel_x_biasList), 2)

		resultingMaxBiases.accel_y = round(max(accel_y_biasList), 2)
		resultingMinBiases.accel_y = round(min(accel_y_biasList), 2)

		resultingMaxBiases.accel_z = round(max(accel_z_biasList), 2)
		resultingMinBiases.accel_z = round(min(accel_z_biasList), 2)

		resultingMaxBiases.mag_x = round(max(mag_x_biasList), -1)
		resultingMinBiases.mag_x = round(min(mag_x_biasList), -1)

		resultingMaxBiases.mag_y = round(max(mag_y_biasList), -1)
		resultingMinBiases.mag_y = round(min(mag_y_biasList), -1)

		resultingMaxBiases.mag_z = round(max(mag_z_biasList), -1)
		resultingMinBiases.mag_z = round(min(mag_z_biasList), -1)

		resultingMaxBiases.baro = round(max(baro_biasList), 0)
		resultingMinBiases.baro= round(min(baro_biasList), 0)

		resultingMaxBiases.pitot = round(max(pitot_biasList), 0)
		resultingMinBiases.pitot= round(min(pitot_biasList), 0)

		if inFailureMode:
			resultingMaxBiases.gyro_x += epsilon

		if (resultingMaxBiases == expectedMaxBiases) and \
			(resultingMinBiases == expectedMinBiases):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, \
					gyroBias, accelBias, magBias, baroBias, pitotBias), \
					(resultingMaxBiases, resultingMinBiases), \
					(expectedMaxBiases, expectedMinBiases))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.initializeSigmas()
currentFuncDef = SM.SensorsModel.initializeSigmas
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for gyroSigma, accelSigma, magSigma, baroSigma, pitotSigma, gpsSigmaHorizontal,\
		gpsSigmaVertical, gpsSigmaSOG, gpsSigmaCOG, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.initializeSigmas(gyroSigma, accelSigma, magSigma, \
			baroSigma,  pitotSigma, gpsSigmaHorizontal, gpsSigmaVertical, \
			gpsSigmaSOG, gpsSigmaCOG)
		
		if inFailureMode:
			result.gyro_x += epsilon

		if result == expected:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (gyroSigma, \
					accelSigma, magSigma, baroSigma, pitotSigma, \
					gpsSigmaHorizontal, gpsSigmaVertical, gpsSigmaSOG, \
					gpsSigmaCOG), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updateGPSTrue()
currentFuncDef = SM.SensorsModel.updateGPSTrue
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, dot, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.updateGPSTrue(state, dot)
		
		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]) and \
		math.isclose(result[3], expected[3]) and \
		math.isclose(result[4], expected[4]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (state, expected,\
					dot), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updateAccelsTrue()
currentFuncDef = SM.SensorsModel.updateAccelsTrue
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, dot, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.updateAccelsTrue(state, dot)
		
		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (state, expected,\
					dot), result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updateMagsTrue()
currentFuncDef = SM.SensorsModel.updateMagsTrue
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.updateMagsTrue(state)
		
		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (state, expected),\
					result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updateGyrosTrue()
currentFuncDef = SM.SensorsModel.updateGyrosTrue
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.updateGyrosTrue(state)
		
		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]) and \
		math.isclose(result[2], expected[2]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (state, expected),\
					result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updatePressureSensorsTrue()
currentFuncDef = SM.SensorsModel.updatePressureSensorsTrue
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for state, expected in curTestBlock:
	try:
		sm = SM.SensorsModel()
		result = sm.updatePressureSensorsTrue(state)
		
		if inFailureMode:
			result = list(result)
			result[0] += epsilon
			result = tuple(result)

		if math.isclose(result[0], expected[0]) and \
		math.isclose(result[1], expected[1]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (state, expected),\
					result, expected)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.updateSensorsNoisy()
currentFuncDef = SM.SensorsModel.updateSensorsNoisy
currentClass = SM.SensorsModel
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, trueSensors, noisySensors, sensorBiases, sensorSigmas, \
	expectedMean in curTestBlock:
	try:
		sm = SM.SensorsModel()

		gyro_x_biasList = list()
		gyro_y_biasList = list()
		gyro_z_biasList = list()

		accel_x_biasList = list()
		accel_y_biasList = list()
		accel_z_biasList = list()

		mag_x_biasList = list()
		mag_y_biasList = list()
		mag_z_biasList = list()

		baro_biasList = list()
		pitot_biasList = list()

		if inFailureMode:
			sensorBiases.gyro_x += epsilon

		for j in range(subTestCount):
			result = sm.updateSensorsNoisy(trueSensors, noisySensors, \
					sensorBiases, sensorSigmas)

			gyro_x_biasList.append(result.gyro_x)
			gyro_y_biasList.append(result.gyro_y)
			gyro_z_biasList.append(result.gyro_z)
			
			accel_x_biasList.append(result.accel_x)
			accel_y_biasList.append(result.accel_y)
			accel_z_biasList.append(result.accel_z)

			mag_x_biasList.append(result.mag_x)
			mag_y_biasList.append(result.mag_y)
			mag_z_biasList.append(result.mag_z)

			baro_biasList.append(result.baro)
			pitot_biasList.append(result.pitot)

		resultingMean = Sensors.vehicleSensors()
		resultingMean.gyro_x = round(sum(gyro_x_biasList)/float(len(gyro_x_biasList)), 0)
		resultingMean.gyro_y = round(sum(gyro_y_biasList)/float(len(gyro_y_biasList)), 0)
		resultingMean.gyro_z = round(sum(gyro_z_biasList)/float(len(gyro_z_biasList)), 0)
		resultingMean.accel_x = round(sum(accel_x_biasList)/float(len(accel_x_biasList)), -1)
		resultingMean.accel_y = round(sum(accel_y_biasList)/float(len(accel_y_biasList)), -1)
		resultingMean.accel_z = round(sum(accel_z_biasList)/float(len(accel_z_biasList)), -1)
		resultingMean.mag_x = round(sum(mag_x_biasList)/float(len(mag_x_biasList)), -1)
		resultingMean.mag_y = round(sum(mag_y_biasList)/float(len(mag_y_biasList)), -1)
		resultingMean.mag_z = round(sum(mag_z_biasList)/float(len(mag_z_biasList)), -1)
		resultingMean.baro = round(sum(baro_biasList)/float(len(baro_biasList)), 0)
		resultingMean.pitot = round(sum(pitot_biasList)/float(len(pitot_biasList)), 0)

		if resultingMean == expectedMean:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, \
					trueSensors, noisySensors, sensorBiases, sensorSigmas),\
					resultingMean, expectedMean)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.update()
currentFuncDef = SM.SensorsModel.update
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
currentFuncDef = SM.SensorsModel.getSensorsTrue
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
currentFuncDef = SM.SensorsModel.getSensorsNoisy
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, aeroModel, expectedTrue, expectedNoisy in curTestBlock:
	try:
		sm = SM.SensorsModel(aeroModel)

		for i in range(subTestCount):
			sm.update()

		resultTrue = sm.getSensorsTrue()
		resultNoisy = sm.getSensorsNoisy()

		resultNoisy.gyro_x = round(resultNoisy.gyro_x, 0)
		resultNoisy.gyro_y = round(resultNoisy.gyro_y, 0)
		resultNoisy.gyro_z = round(resultNoisy.gyro_z, 0)
		resultNoisy.accel_x = round(resultNoisy.accel_x, -2)
		resultNoisy.accel_y = round(resultNoisy.accel_y, -2)
		resultNoisy.accel_z = round(resultNoisy.accel_z, -2)
		resultNoisy.mag_x = round(resultNoisy.mag_x, -5)
		resultNoisy.mag_y = round(resultNoisy.mag_y, -5)
		resultNoisy.mag_z = round(resultNoisy.mag_z, -5)
		resultNoisy.baro = round(resultNoisy.baro, -3)
		resultNoisy.pitot = round(resultNoisy.pitot, -3)
		resultNoisy.gps_n = round(resultNoisy.gps_n, -1)
		resultNoisy.gps_e = round(resultNoisy.gps_e, -1)
		resultNoisy.gps_alt = round(resultNoisy.gps_alt, -2)
		resultNoisy.gps_sog = round(resultNoisy.gps_sog, -2)
		resultNoisy.gps_cog = round(resultNoisy.gps_cog, -2)
		
		if inFailureMode:
			resultNoisy.gyro_x += epsilon

		if (resultTrue == expectedTrue) and (resultNoisy == expectedNoisy):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, \
					aeroModel), (resultTrue, resultNoisy), (expectedTrue, expectedNoisy))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# SM.SensorsModel.reset()
currentClass = SM.SensorsModel
currentFuncDef = SM.SensorsModel.reset
print("Comparing outputs for {}.{}()".format(currentClass.__name__, currentFuncDef.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for subTestCount, aeroModel, expectedTrue, expectedNoisy in curTestBlock:
	try:
		sm = SM.SensorsModel(aeroModel)

		for i in range(subTestCount):
			sm.update()

		sm.reset()

		resultTrue = sm.getSensorsTrue()
		resultNoisy = sm.getSensorsNoisy()
		
		if inFailureMode:
			resultNoisy.gyro_x += epsilon

		if (resultTrue == expectedTrue) and (resultNoisy == expectedNoisy):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(currentClass, currentFuncDef, (subTestCount, \
					aeroModel), (resultTrue, resultNoisy), (expectedTrue, expectedNoisy))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(currentClass, currentFuncDef, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Sensors Model")
else:
	print("{}/{} tests blocks passed for Sensors Model".\
		format(testBlocksPassed, len(allTests)))