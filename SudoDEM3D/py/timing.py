# encoding: utf-8
# 2008 © Václav Šmilauer <eudoxos@arcig.cz>
"""Functions for accessing timing information stored in engines and functors.

See :ref:`timing` section of the programmer's manual, `wiki page <http://sudodem-dem.org/index.php/Speed_profiling_using_TimingInfo_and_TimingDeltas_classes>`_ for some examples.

"""

from sudodem.wrapper import *


def _resetEngine(e):
	if e.timingDeltas: e.timingDeltas.reset()
	if isinstance(e,Functor): return
	if isinstance(e,Dispatcher):
		for f in e.functors: _resetEngine(f)
	elif isinstance(e,ParallelEngine):
		for s in e.slaves: _resetEngine(s)
	e.execTime,e.execCount=0,0

def reset():
	"Zero all timing data."
	for e in O.engines: _resetEngine(e)

_statCols={'label':40,'count':20,'time':20,'relTime':20}
_maxLev=3

def _formatLine(label,time,count,totalTime,level):
	sp,negSp=' '*level*2,' '*(_maxLev-level)*2
	raw=[]
	raw.append(label)
	raw.append(str(count) if count>=0 else '')
	raw.append((str(time/1000)+u'us') if time>=0 else '')
	raw.append(('%6.2f%%'%(time*100./totalTime)) if totalTime>0 else '')
	return u' '.join([
		(sp+raw[0]).ljust(_statCols['label']),
		(raw[1]+negSp).rjust(_statCols['count']),
		(raw[2]+negSp).rjust(_statCols['time']),
		(raw[3]+negSp).rjust(_statCols['relTime']),
	])

def _delta_stats(deltas,totalTime,level):
	ret=0
	deltaTime=sum([d[1] for d in deltas.data])
	for d in deltas.data:
		print _formatLine(d[0],d[1],d[2],totalTime,level); ret+=1
	if len(deltas.data)>1:
		print _formatLine('TOTAL',deltaTime,sum(d[2] for d in deltas.data),totalTime,level); ret+=1
	return ret

def _engines_stats(engines,totalTime,level):
	lines=0; hereLines=0
	for e in engines:
		if not isinstance(e,Functor): print _formatLine(u'"'+e.label+'"' if e.label else e.__class__.__name__,e.execTime,e.execCount,totalTime,level); lines+=1; hereLines+=1
		if e.timingDeltas:
			if isinstance(e,Functor):
				print _formatLine(e.__class__.__name__,sum(d[1] for d in e.timingDeltas.data),sum(d[2] for d in e.timingDeltas.data),totalTime,level); lines+=1; hereLines+=1
				execTime=sum([d[1] for d in e.timingDeltas.data])
			else: execTime=e.execTime
			lines+=_delta_stats(e.timingDeltas,execTime,level+1)
		if isinstance(e,Dispatcher): lines+=_engines_stats(e.functors,e.execTime,level+1)
		if isinstance(e,InteractionLoop):
			lines+=_engines_stats(e.geomDispatcher.functors,e.execTime,level+1)
			lines+=_engines_stats(e.physDispatcher.functors,e.execTime,level+1)
			lines+=_engines_stats(e.lawDispatcher.functors,e.execTime,level+1)
		elif isinstance(e,ParallelEngine): lines+=_engines_stats(e.slave,e.execTime,level+1)
	if hereLines>1 and not isinstance(e,Functor):
		print _formatLine('TOTAL',totalTime,-1,totalTime,level); lines+=1
	return lines

def stats():
	"""Print summary table of timing information from engines and functors. Absolute times as well as percentages are given. Sample output:

	.. code-block:: none

		Name                                                    Count                 Time            Rel. time
		-------------------------------------------------------------------------------------------------------
		ForceResetter                                       102               2150us                0.02%
		"collider"                                            5              64200us                0.60%
		InteractionLoop                                     102           10571887us               98.49%
		"combEngine"                                        102               8362us                0.08%
		"newton"                                            102              73166us                0.68%
		"cpmStateUpdater"                                     1               9605us                0.09%
		PyRunner                                              1                136us                0.00%
		"plotDataCollector"                                   1                291us                0.00%
		TOTAL                                                             10733564us              100.00%


	sample output (compiled with -DCMAKE_CXX_FLAGS="-DUSE_TIMING_DELTAS" option):

	.. code-block:: none

		Name                                                    Count                 Time            Rel. time
		-------------------------------------------------------------------------------------------------------
		ForceResetter                                       102               2150us                0.02%
		"collider"                                            5              64200us                0.60%
		InteractionLoop                                     102           10571887us               98.49%
		  Ig2_Sphere_Sphere_ScGeom                        1222186            1723168us               16.30%
		    Ig2_Sphere_Sphere_ScGeom                        1222186            1723168us              100.00%
		  Ig2_Facet_Sphere_ScGeom                             753               1157us                0.01%
		    Ig2_Facet_Sphere_ScGeom                             753               1157us              100.00%
		  Ip2_CpmMat_CpmMat_CpmPhys                         11712              26015us                0.25%
		    end of Ip2_CpmPhys                                11712              26015us              100.00%
		  Ip2_FrictMat_CpmMat_FrictPhys                         0                  0us                0.00%
		  Law2_ScGeom_CpmPhys_Cpm                         3583872            4819289us               45.59%
		    GO A                                            1194624            1423738us               29.54%
		    GO B                                            1194624            1801250us               37.38%
		    rest                                            1194624            1594300us               33.08%
		    TOTAL                                           3583872            4819289us              100.00%
		  Law2_ScGeom_FrictPhys_CundallStrack                   0                  0us                0.00%
		"combEngine"                                        102               8362us                0.08%
		"newton"                                            102              73166us                0.68%
		"cpmStateUpdater"                                     1               9605us                0.09%
		PyRunner                                              1                136us                0.00%
		"plotDataCollector"                                   1                291us                0.00%
		TOTAL                                                             10733564us              100.00%

	"""

	print 'Name'.ljust(_statCols['label'])+' '+'Count'.rjust(_statCols['count'])+' '+'Time'.rjust(_statCols['time'])+' '+'Rel. time'.rjust(_statCols['relTime'])
	print '-'*(sum([_statCols[k] for k in _statCols])+len(_statCols)-1)
	_engines_stats(O.engines,sum([e.execTime for e in O.engines]),0)
	print
