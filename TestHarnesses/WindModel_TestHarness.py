"""
.. modle:: WindModel_TestHarness.py
	:platform: MacOS, Unix, Windows,
	:synopsis: Compares output from WindModel_Generation.py with students'  
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
import ece163.Containers.States as States
import ece163.Containers.Inputs as Inputs
import ece163.Modeling.WindModel as WM
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
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all tests regardless of failures')
parser.add_argument('picklePath', nargs='?', default='WindModel_TestData.pickle', help='valid path to pickle for input')
parser.add_argument('-f','--failure', action='store_true', help='Purposely fail all tests to ensure failures modes do not have errors')
arguments = parser.parse_args()


picklePath = arguments.picklePath
inContinueMode = arguments.continueMode
inFailureMode = arguments.failure

if inFailureMode:
	epsilon = 1
else:
	epsilon = 0

print("Beginning Test Harness for Wind Model using file {}".format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed
testBlockIterator = iter(allTests)  # we hard code the tests as well so we need an iterator

###############################################################################
# setWind() and getWind()
print("Comparing outputs for {}() and {}()".format(WM.WindModel.setWind.__name__, WM.WindModel.getWind.__name__))
curTestBlock = next(testBlockIterator) # we now have all the tests for this block
testsPassed = 0
for state, expectedState in curTestBlock: 
	try:
		wind = WM.WindModel()

		wind.setWind(state)
		resultingState = wind.getWind()

		if inFailureMode:
			resultingState.Wn += epsilon

		if resultingState == expectedState:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(wind.setWind, state, resultingState, expectedState)
				printTestFailure(wind.getWind, state, resultingState, expectedState)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(wind.setWind, testsPassed, len(curTestBlock))
printTestBlockResult(wind.getWind, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1


###############################################################################
# CreateDrydenTransferFns()
print("Comparing outputs for {}()".format(WM.WindModel.CreateDrydenTransferFns.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for dT, Va, dryden, Phi_u, Gamma_u, H_u, Phi_v, Gamma_v, H_v, Phi_w, Gamma_w, \
H_w in curTestBlock:  # and now we can iterate through them
	try:
		resGustModel = WM.WindModel(dT, Va, dryden)

		res_Phi_u, res_Gamma_u, res_H_u, res_Phi_v, res_Gamma_v, res_H_v, \
		res_Phi_w, res_Gamma_w, res_H_w,  = resGustModel.getDrydenTransferFns()

		if inFailureMode:
			res_Phi_u[0][0] += epsilon

		eqFlag = list(range(9))
		eqFlag[0], et, rt = matrixCompare(res_Phi_u, Phi_u)
		eqFlag[1], et, rt = matrixCompare(res_Gamma_u, Gamma_u)
		eqFlag[2], et, rt = matrixCompare(res_H_u, H_u)

		eqFlag[3], et, rt = matrixCompare(res_Phi_v, Phi_v)
		eqFlag[4], et, rt = matrixCompare(res_Gamma_v, Gamma_v)
		eqFlag[5], et, rt = matrixCompare(res_H_v, H_v)

		eqFlag[6], et, rt = matrixCompare(res_Phi_w, Phi_w)
		eqFlag[7], et, rt = matrixCompare(res_Gamma_w, Gamma_w)
		eqFlag[8], et, rt = matrixCompare(res_H_w, H_w)

		resultingOutput = ("Phi_u=", res_Phi_u, "Gamma_u=", \
			res_Gamma_u, "H_u=", res_H_u, "Phi_v=", \
			res_Phi_v, "Gamma_v=", res_Gamma_v, "H_v=", \
			res_H_v, "Phi_w=", res_Phi_w, "Gamma_w=", \
			res_Gamma_w, "H_w=", res_H_w)

		expectedOutput = ("Phi_u=", Phi_u, "Gamma_u=", Gamma_u, "H_u=", H_u, \
			"Phi_v=", Phi_v, "Gamma_v=", Gamma_v, "H_v=", H_v, "Phi_w=", Phi_w,\
			"Gamma_w=", Gamma_w, "H_w=", H_w)

		if all(elems == True for elems in eqFlag):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(wind.CreateDrydenTransferFns, (dT, Va, dryden),\
					resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(wind.CreateDrydenTransferFns, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# Update()
print("Comparing outputs for {}()".format(WM.WindModel.Update.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for windState, dT, Va, dryden, uu, uv, uw, UpdateTestCount, x_u, x_v, x_w, expectedWind in curTestBlock:
	try:
		resGustModel = WM.WindModel(dT, Va, dryden)
		resGustModel.setWind(windState)

		for i in range(UpdateTestCount):
			resGustModel.Update(uu, uv, uw)

		resultingWind = resGustModel.getWind()
		if inFailureMode:
			resGustModel.x_u[0][0] += epsilon

		eqFlag = list(range(3))
		eqFlag[0], et, rt = matrixCompare(resGustModel.x_u, x_u)
		eqFlag[1], et, rt = matrixCompare(resGustModel.x_v, x_v)
		eqFlag[2], et, rt = matrixCompare(resGustModel.x_w, x_w)


		resultingOutput = ("x_u=", resGustModel.x_u, "x_v=", resGustModel.x_v, \
			"x_w=", resGustModel.x_w)

		expectedOutput = ("x_u=",x_u, "x_v=",x_v, "x_w=",x_w)

		if all(elems == True for elems in eqFlag) and \
		(resultingWind == expectedWind):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(wind.Update, (dT, Va, dryden),\
					resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(wind.Update, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# reset()
print("Comparing outputs for {}()".format(WM.WindModel.reset.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0

for initialGustModel, expGustModel in curTestBlock:
	try:
		resGustModel = copy.deepcopy(initialGustModel)
		resGustModel.reset()

		res_Phi_u, res_Gamma_u, res_H_u, res_Phi_v, res_Gamma_v, res_H_v, \
		res_Phi_w, res_Gamma_w, res_H_w,  = resGustModel.getDrydenTransferFns()

		exp_Phi_u, exp_Gamma_u, exp_H_u, exp_Phi_v, exp_Gamma_v, exp_H_v, \
		exp_Phi_w, exp_Gamma_w, exp_H_w,  = expGustModel.getDrydenTransferFns()

		if inFailureMode:
			resGustModel.x_u[0][0] += epsilon

		eqFlag = list(range(12))
		eqFlag[0], et, rt = matrixCompare(resGustModel.x_u, expGustModel.x_u)
		eqFlag[1], et, rt = matrixCompare(resGustModel.x_v, expGustModel.x_v)
		eqFlag[2], et, rt = matrixCompare(resGustModel.x_w, expGustModel.x_w)

		eqFlag[3], et, rt = matrixCompare(res_Phi_u, exp_Phi_u)
		eqFlag[4], et, rt = matrixCompare(res_Gamma_u, exp_Gamma_u)
		eqFlag[5], et, rt = matrixCompare(res_H_u, exp_H_u)

		eqFlag[6], et, rt = matrixCompare(res_Phi_v, exp_Phi_v)
		eqFlag[7], et, rt = matrixCompare(res_Gamma_v, exp_Gamma_v)
		eqFlag[8], et, rt = matrixCompare(res_H_v, exp_H_v)

		eqFlag[9], et, rt = matrixCompare(res_Phi_w, exp_Phi_w)
		eqFlag[10], et, rt = matrixCompare(res_Gamma_w, exp_Gamma_w)
		eqFlag[11], et, rt = matrixCompare(res_H_w, exp_H_w)

		resultingOutput = ("x_u=",resGustModel.x_u, "x_v=",resGustModel.x_v, \
			"x_w=",resGustModel.x_w, "Phi_u=", res_Phi_u, "Gamma_u=", \
			res_Gamma_u, "H_u=", res_H_u, "Phi_v=", \
			res_Phi_v, "Gamma_v=", res_Gamma_v, "H_v=", \
			res_H_v, "Phi_w=", res_Phi_w, "Gamma_w=", \
			res_Gamma_w, "H_w=", res_H_w)

		expectedOutput = ("x_u=",expGustModel.x_u, "x_v=",expGustModel.x_v, \
			"x_w=",expGustModel.x_w, "Phi_u=", exp_Phi_u, "Gamma_u=", \
			exp_Gamma_u, "H_u=", exp_H_u, "Phi_v=", \
			exp_Phi_v, "Gamma_v=", exp_Gamma_v, "H_v=", \
			exp_H_v, "Phi_w=", exp_Phi_w, "Gamma_w=", \
			exp_Gamma_w, "H_w=", exp_H_w)

		if all(elems == True for elems in eqFlag) and \
		(resGustModel.getWind() == expGustModel.Wind):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(wind.reset, (initialGustModel),\
					resultingOutput, expectedOutput)
				if not inFailureMode:
					sys.exit(-1)

	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(wind.reset, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Wind Model")
else:
	print("{}/{} tests blocks passed for Wind Model".format(testBlocksPassed, len(allTests)))

