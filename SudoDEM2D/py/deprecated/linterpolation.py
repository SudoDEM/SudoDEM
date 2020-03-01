# encoding: utf-8
#
# © 2009 Václav Šmilauer <eudoxos@arcig.cz>
#

"""
Module for rudimentary support of manipulation with piecewise-linear
functions (which are usually interpolations of higher-order functions,
whence the module name). Interpolation is always given as two lists
of the same length, where the x-list must be increasing.

Periodicity is supported by supposing that the interpolation
can wrap from the last x-value to the first x-value (which
should be 0 for meaningful results).

Non-periodic interpolation can be converted to periodic one
by padding the interpolation with constant head and tail using
the sanitizeInterpolation function.

There is a c++ template function for interpolating on such sequences in
pkg/common/Engine/PartialEngine/LinearInterpolate.hpp (stateful, therefore
fast for sequential reads).

TODO: Interpolating from within python is not (yet) supported.
"""

def revIntegrateLinear(I,x0,y0,x1,y1):
	"""Helper function, returns value of integral variable x for
	linear function f passing through (x0,y0),(x1,y1) such that
	1. x∈[x0,x1]
	2. ∫_x0^x f dx=I
	and raise exception if such number doesn't exist or the solution
	is not unique (possible?)
	"""
	from math import sqrt
	dx,dy=x1-x0,y1-y0
	if dy==0: # special case, degenerates to linear equation
		return (x0*y0+I)/y0
	a=dy/dx; b=2*(y0-x0*dy/dx); c=x0**2*dy/dx-2*x0*y0-2*I
	det=b**2-4*a*c; assert(det>0)
	p,q=(-b+sqrt(det))/(2*a),(-b-sqrt(det))/(2*a)
	pOK,qOK=x0<=p<=x1,x0<=q<=x1
	if pOK and qOK: raise ValueError("Both radices within interval!?")
	if not pOK and not qOK: raise ValueError("No radix in given interval!")
	return p if pOK else q

def integral(x,y):
	"""Return integral of piecewise-linear function given by points
	x0,x1,… and y0,y1,… """
	assert(len(x)==len(y))
	sum=0
	for i in range(1,len(x)): sum+=(x[i]-x[i-1])*.5*(y[i]+y[i-1])
	return sum

def xFractionalFromIntegral(integral,x,y):
	"""Return x within range x0…xn such that ∫_x0^x f dx==integral.
	Raises error if the integral value is not reached within the x-range.
	"""
	sum=0
	for i in range(1,len(x)):
		diff=(x[i]-x[i-1])*.5*(y[i]+y[i-1])
		if sum+diff>integral:
			return revIntegrateLinear(integral-sum,x[i-1],y[i-1],x[i],y[i])
		else: sum+=diff
	raise "Integral not reached within the interpolation range!"


def xFromIntegral(integralValue,x,y):
	"""Return x such that ∫_x0^x f dx==integral.
	x wraps around at xn. For meaningful results, therefore, x0 should == 0 """
	from math import floor
	period=x[-1]-x[0]
	periodIntegral=integral(x,y)
	numPeriods=floor(integralValue/periodIntegral)
	xFrac=xFractionalFromIntegral(integralValue-numPeriods*periodIntegral,x,y)
	#print '### wanted _%g; period=%g; periodIntegral=_%g (numPeriods=%g); rests _%g (xFrac=%g)'%(integralValue,period,periodIntegral,numPeriods,integralValue-numPeriods*periodIntegral,xFrac)
	#print '### returning %g*%g+%g=%g'%(period,numPeriods,xFrac,period*numPeriods+xFrac)
	return period*numPeriods+xFrac

def sanitizeInterpolation(x,y,x0,x1):
	"""Extends piecewise-linear function in such way that it spans at least
	the x0…x1 interval, by adding constant padding at the beginning (using y0)
	and/or at the end (using y1) or not at all."""
	xx,yy=[],[]
	if x0<x[0]:
		xx+=[x0]; yy+=[y[0]]
	xx+=x; yy+=y
	if x1>x[-1]:
		xx+=[x1]; yy+=[y[-1]]
	return xx,yy

if __name__=="main":
	xx,yy=sanitizeInterpolation([1,2,3],[1,1,2],0,4)
	print xx,yy
	print integral(xx,yy) # 5.5
	print revIntegrateLinear(.625,1,1,2,2) # 1.5
	print xFractionalFromIntegral(1.625,xx,yy) # 1.625
	print xFractionalFromIntegral(2.625,xx,yy) # 2.5

