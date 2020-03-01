# encoding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>
"""All defined functionality tests for sudodem."""
import unittest,inspect,sys

# add any new test suites to the list here, so that they are picked up by testAll
allTests=['wrapper','core','pbc','clump','cohesive-chain']

# all sudodem modules (ugly...)
import sudodem.export,sudodem.linterpolation,sudodem.pack,sudodem.plot,sudodem.post2d,sudodem.timing,sudodem.utils,sudodem.ymport,sudodem.geom
allModules=(sudodem.export,sudodem.linterpolation,sudodem.pack,sudodem.plot,sudodem.post2d,sudodem.timing,sudodem.utils,sudodem.ymport,sudodem.geom)
try:
	import sudodem.qt
	allModules+=(sudodem.qt,)
except ImportError: pass

# fully qualified module names
allTestsFQ=['sudodem.tests.'+test for test in allTests]

def testModule(module):
	"""Run all tests defined in the module specified, return TestResult object
	(http://docs.python.org/library/unittest.html#unittest.TextTestResult)
	for further processing.

	@param module: fully-qualified module name, e.g. sudodem.tests.wrapper
	"""
	suite=unittest.defaultTestLoader().loadTestsFromName(module)
	return unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)

def testAll():
	"""Run all tests defined in all sudodem.tests.* modules and return
	TestResult object for further examination."""
	suite=unittest.defaultTestLoader.loadTestsFromNames(allTestsFQ)
	import doctest
	for mod in allModules:
		suite.addTest(doctest.DocTestSuite(mod))
	return unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)



