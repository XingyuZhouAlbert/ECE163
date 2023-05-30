"""
.. modle:: PerturbationModel_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from PerturbationModel_Generation.py with 
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
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Linearized
import ece163.Modeling.WindModel as WindModel
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
parser.add_argument('picklePath', nargs='?', default='PerturbationModel_TestData.pickle', help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()

picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Perturbation Model using file {}".format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed

testBlockIterator = iter(allTests)  # we hard code the tests as well so we need an iterator

###############################################################################
# dThrust_dVa()
print("Comparing outputs for {}()".format(VPM.dThrust_dVa.__name__))
curTestBlock = next(testBlockIterator)
testsPassed = 0

for Va, Throttle, ep, expectedOutput in curTestBlock: 
	try:
		if (ep != None):
			result = VPM.dThrust_dVa(Va, Throttle, ep)
		else:
			result = VPM.dThrust_dVa(Va, Throttle)

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expectedOutput):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(VPM.dThrust_dVa, (Va, Throttle), result, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(VPM.dThrust_dVa, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# dThrust_dThrottle()
print("Comparing outputs for {}()".format(VPM.dThrust_dThrottle.__name__))
curTestBlock = next(testBlockIterator)
testsPassed = 0

for Va, Throttle, ep, expectedOutput in curTestBlock: 
	try:
		if (ep != None):
			result = VPM.dThrust_dThrottle(Va, Throttle, ep)
		else:
			result = VPM.dThrust_dThrottle(Va, Throttle)

		if inFailureMode:
			result += epsilon

		if math.isclose(result, expectedOutput):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(VPM.dThrust_dThrottle, (Va, Throttle), result, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(VPM.dThrust_dThrottle, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# CreateTransferFunction()
print("Comparing outputs for {}()".format(VPM.CreateTransferFunction.__name__))
curTestBlock = next(testBlockIterator)
testsPassed = 0

for Va, Kappa, Gamma, expect in curTestBlock: 
	try:
		vTrim = VT.VehicleTrim()
		validInput = vTrim.computeTrim(Va, Kappa, Gamma)

		trimState = vTrim.getTrimState()
		trimInputs = vTrim.getTrimControls()

		if validInput:
			result = VPM.CreateTransferFunction(trimState, trimInputs)
		else:
			result = Linearized.transferFunctions()

		resultRep = ("Va_trim=",result.Va_trim, "alpha_trim=",result.alpha_trim,\
			"beta_trim=",result.beta_trim, "gamma_trim=",result.gamma_trim,\
			"theta_trim=",result.theta_trim, "phi_trim=",result.phi_trim,\
			"a_phi1=",result.a_phi1, "a_phi2=", result.a_phi2, "a_beta1=", \
			result.a_beta1, "a_beta2=", result.a_beta2, "a_theta1=", \
			result.a_theta1, "a_theta2=", result.a_theta2, "a_theta3=", \
			result.a_theta3, "a_V1=", result.a_V1, "a_V2=", result.a_V2,\
			"a_V3=", result.a_V3)

		expectRep = ("Va_trim=",expect.Va_trim, "alpha_trim=",expect.alpha_trim\
			,"beta_trim=",expect.beta_trim, "gamma_trim=",expect.gamma_trim,\
			"theta_trim=",expect.theta_trim, "phi_trim=",expect.phi_trim,\
			"a_phi1=",expect.a_phi1, "a_phi2=", expect.a_phi2, "a_beta1=", \
			expect.a_beta1, "a_beta2=", expect.a_beta2, "a_theta1=", \
			expect.a_theta1, "a_theta2=", expect.a_theta2, "a_theta3=", \
			expect.a_theta3, "a_V1=", expect.a_V1, "a_V2=", expect.a_V2,\
			"a_V3=", expect.a_V3)

		if inFailureMode:
			result.Va_trim += epsilon

		if result == expect: # comparison is valid here, __eq__ is overloaded, see Linearized.py
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(VPM.CreateTransferFunction, (Va, Throttle), resultRep, expectRep)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(VPM.CreateTransferFunction, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Perturbation Model")
else:
	print("{}/{} tests blocks passed for Perturbation Model".format(testBlocksPassed, len(allTests)))

