"""
.. modle:: matrixMathTestHarness.py
    :platform: MacOS, Unix, Windows,
    :synopsis: Compares output from matrixMathGeneration.py with students'  
    function implementations
    matrixMathTestHarness.py
.. moduleauthor:: Pavlo Vlastos <pvlastos@ucsc.edu>
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import argparse
import random
import pickle
import traceback

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
					print("Element [{0},{1}] is incorrect".format(row, col))
					return [False, expTot, resTot]
				else:
					resTot += 1
		if expTot != resTot:
			print("\r\nResulting matrix dimensions match the expected matrix dimensions")
			return [True, expTot, resTot]
	else:
		print("Error: Resulting matrix dimensions do not match the expected matrix dimensions")
		return [False, expTot, resTot]

	return [True, expTot, resTot]

#  these two functions allow for more standardized output, they should be copied to each test harness and customized
def printTestBlockResult(function, testsPassed, testCount):
	if testsPassed != testCount:
		addendum = " (TESTS FAILED)"
	else:
		addendum = ""
	print("{}/{} tests passed for {}{}".format(testsPassed, testCount, function.__name__, addendum))

def printTestFailure(function, inputs, outputs, expectedoutputs):
	print("Test Failed for {}. Please find the corresponding inputs and expected output(s) below for testing".format(function.__name__))
	print("Inputs: {}".format(repr(inputs)))
	print("Outputs: {}".format(repr(outputs)))
	print("Expected Outputs: {}".format(repr(expectedoutputs)))


###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c','--continueMode', action='store_true', help='Runs all tests regardless of failures')
parser.add_argument('picklePath', nargs='?', default='MatrixMath_TestData.pickle', help='valid path to pickle for input')
arguments = parser.parse_args()

picklePath = arguments.picklePath
inContinueMode = arguments.continueMode

print("Beginning Test Harness for Matirx Math Library using file {}".format(picklePath))
try:
	with open(picklePath, 'rb') as f:
		allTests = pickle.load(f)
except FileNotFoundError:
	print('Test file not found, exiting')
	sys.exit(-1)

testBlocksPassed = 0  # we keep track of the number of test blocks passed

testBlockIterator = iter(allTests)  # we hard code the tests as well so we need an iterator

###############################################################################
# matrixMultiply()
print("Comparing outputs for {}".format(MatrixMath.matrixMultiply.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixMultiply(A, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixMultiply, (A, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixMultiply, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixTranspose()
print("Comparing outputs for {}".format(MatrixMath.matrixTranspose.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixTranspose(A)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixTranspose, A, result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixTranspose, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixAdd()
print("Comparing outputs for {}".format(MatrixMath.matrixAdd.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixAdd(A, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixAdd, (A, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixAdd, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixSubtract()
print("Comparing outputs for {}".format(MatrixMath.matrixSubtract.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixSubtract(A, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixSubtract, (A, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixSubtract, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixScalarMultiply()
print("Comparing outputs for {}".format(MatrixMath.matrixScalarMultiply.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for a, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixScalarMultiply(a, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixScalarMultiply, (a, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixScalarMultiply, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixDotProduct()
print("Comparing outputs for {}".format(MatrixMath.matrixDotProduct.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for a, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixDotProduct(a, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixDotProduct, (A, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixDotProduct, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixSkew()
print("Comparing outputs for {}".format(MatrixMath.matrixSkew.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for a, b, c, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixSkew(a, b, c)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixSkew, (a, b, c), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixSkew, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixCrossProduct()
print("Comparing outputs for {}".format(MatrixMath.matrixCrossProduct.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, B, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixCrossProduct(A, B)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixCrossProduct, (A, B), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixCrossProduct, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixOffset()
print("Comparing outputs for {}".format(MatrixMath.matrixOffset.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, x, y, z, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixOffset(A, x, y, z)
		[matricesEqual, expTot, resTot] = matrixCompare(result, expectedResult)
		if matricesEqual is True:
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixOffset, (A, x, y, z), result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixOffset, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
# matrixSize()
print("Comparing outputs for {}".format(MatrixMath.matrixSize.__name__))
curTestBlock = next(testBlockIterator)  # we now have all the tests for this block
testsPassed = 0
for A, expectedResult in curTestBlock:  # and now we can iterate through them
	try:
		result = MatrixMath.matrixSize(A)

		if (expectedResult[0] == result[0]) and (expectedResult[1] == result[1]):
			testsPassed += 1
		else:
			if not inContinueMode:
				printTestFailure(MatrixMath.matrixSize, A, result, expectedResult)
				sys.exit(-1)
	except Exception as e:  # overly broad exception clause but we want it to catch everything
		print('Test harness failed with the exception below. It will not continue')
		print(traceback.format_exc())
		sys.exit(-1)

printTestBlockResult(MatrixMath.matrixSize, testsPassed, len(curTestBlock))
if testsPassed == len(curTestBlock):
	testBlocksPassed += 1

###############################################################################
if testBlocksPassed == len(allTests):
	print("All tests Passed for Matrix Math")
else:
	print("{}/{} tests blocks passed for Matrix Math".format(testBlocksPassed, len(allTests)))