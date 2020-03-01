1# this module is populated at initialization from the c++ part of PythonUI
"""Runtime variables, populated at sudodem startup."""
# default value
hasDisplay=False


# find out about which ipython version we use -- 0.10* and 0.11 are supported, but they have different internals
import IPython
try: # attempt to get numerical version
	ipython_version=int(IPython.__version__.split('.')[0])*100 + int(IPython.__version__.split('.')[1])
except ValueError:
	print 'WARN: unable to extract IPython version from %s, defaulting to 10'%(IPython.__version__)
	ipython_version=10
if (ipython_version < 10): #set version 10 for very old systems
	newipver=10
	ipython_version=newipver
