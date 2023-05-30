"""
.. modle:: ControlGains_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from ControlGains_Generation.py with students'  
	function implementations
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
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Controls as Controls
import ece163.Containers.Linearized as Linearized
import ece163.Constants.VehiclePhysicalConstants as VPC
import ece163.Controls.VehicleControlGains as CG
import ece163.Controls.VehiclePerturbationModels as VPM
import ece163.Controls.VehicleTrim as VT
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

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all tests regardless of failures')
parser.add_argument('picklePath', nargs='?', default='ControlGains_TestData.pickle', help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()


picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Control Gains using file {}".format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed
testBlockIterator = iter(allTests)  # we hard code the tests as well so we need an iterator

###############################################################################
# computeGains()
print("Comparing outputs for {}()".format(CG.computeGains.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for Va, Kappa, Gamma, ct, expectedOutput in curTestBlock:
	try:
		vTrim = VT.VehicleTrim()
		validInput = vTrim.computeTrim(Va, Kappa, Gamma)

		if validInput:
			trimState = vTrim.getTrimState()
			trimInputs = vTrim.getTrimControls()
			tf = VPM.CreateTransferFunction(trimState, trimInputs)

			result = CG.computeGains(ct, tf)

		if inFailureMode:
			result.kp_roll += epsilon

		if result == expectedOutput:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CG.computeGains, (Va, Kappa, Gamma, str(ct)),\
					str(result), str(expectedOutput))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CG.computeGains, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# computeTuningParameters()
print("Comparing outputs for {}()".format(CG.computeTuningParameters.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for Va, Kappa, Gamma, gains, expectedOutput in curTestBlock:
	try:
		vTrim = VT.VehicleTrim()
		validInput = vTrim.computeTrim(Va, Kappa, Gamma)

		if validInput:
			trimState = vTrim.getTrimState()
			trimInputs = vTrim.getTrimControls()
			tf = VPM.CreateTransferFunction(trimState, trimInputs)

			result = CG.computeTuningParameters(gains, tf)

		if inFailureMode:
			result.Wn_roll += epsilon

		if result == expectedOutput:
				testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(CG.computeTuningParameters, (Va, Kappa, Gamma, str(ct), str(gains)),\
					str(result), str(expectedOutput))
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(CG.computeTuningParameters, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Control Gains")
else:
	print("{}/{} tests blocks passed for Control Gains".format(testBlocksPassed, len(allTests)))

