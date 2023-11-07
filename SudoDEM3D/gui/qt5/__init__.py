# encoding: utf-8
import sudodem.runtime
if not sudodem.runtime.hasDisplay:
	msg = "Connecting to DISPLAY at SudoDEM startup failed, unable to activate the qt4 interface."
	import os
	if 'SUDODEM_BATCH' in os.environ:
		msg += "\nDo not import qt when running in batch mode."
	raise ImportError(msg)

from PyQt5.QtGui import *
from PyQt5 import QtCore

from sudodem.qt.ui_controller import Ui_Controller

from sudodem.qt.Inspector import *
from sudodem import *
import sudodem.system

from sudodem.qt._GLViewer import *

maxWebWindows=1
"Number of webkit windows that will be cycled to show help on clickable objects"
webWindows=[]
"holds instances of QtWebKit windows; clicking an url will open it in the window that was the least recently updated"
sphinxOnlineDocPath='https://www.sudodem.org/doc/'
"Base URL for the documentation. Packaged versions should change to the local installation directory."


import os.path
# find if we have docs installed locally from package
sphinxLocalDocPath=os.environ['SUDODEM_PREFIX']+'/share/doc/sudodem-doc/html/'
#sphinxBuildDocPath=sudodem.config.sourceRoot+'/doc/sphinx/_build/html/'
# we prefer the packaged documentation for this version, if installed
if   os.path.exists(sphinxLocalDocPath+'/index.html'): sphinxPrefix='file://'+sphinxLocalDocPath
# otherwise look for documentation generated in the source tree
#elif  os.path.exists(sphinxBuildDocPath+'/index.html'): sphinxPrefix='file://'+sphinxBuildDocPath
# fallback to online docs
else: sphinxPrefix=sphinxOnlineDocPath

sphinxDocWrapperPage=sphinxPrefix+'/sudodem.wrapper.html'




def openUrl(url):
	from PyQt5 import QtWebKit
	global maxWebWindows,webWindows
	reuseLast=False
	# use the last window if the class is the same and only the attribute differs
	try:
		reuseLast=(len(webWindows)>0 and str(webWindows[-1].url()).split('#')[-1].split('.')[2]==url.split('#')[-1].split('.')[2])
		#print str(webWindows[-1].url()).split('#')[-1].split('.')[2],url.split('#')[-1].split('.')[2]
	except: pass
	if not reuseLast:
		if len(webWindows)<maxWebWindows: webWindows.append(QtWebKit.QWebView())
		else: webWindows=webWindows[1:]+[webWindows[0]]
	web=webWindows[-1]
	web.load(QUrl(url)); web.setWindowTitle(url);
	if 0:
		def killSidebar(result):
			frame=web.page().mainFrame()
			frame.evaluateJavaScript("var bv=$('.bodywrapper'); bv.css('margin','0 0 0 0');")
			frame.evaluateJavaScript("var sbw=$('.sphinxsidebarwrapper'); sbw.css('display','none');")
			frame.evaluateJavaScript("var sb=$('.sphinxsidebar'); sb.css('display','none'); ")
			frame.evaluateJavaScript("var sb=$('.sidebar'); sb.css('width','0px'); ")
			web.loadFinished.disconnect(killSidebar)
		web.loadFinished.connect(killSidebar)
	web.show();	web.raise_()


controller=None

class ControllerClass(QWidget,Ui_Controller):
	def __init__(self,parent=None):
		QWidget.__init__(self)
		self.setupUi(self)
		self.generator=None # updated automatically
		self.renderer=Renderer() # only hold one instance, managed by OpenGLManager
		self.addPreprocessors()
		self.addRenderers()
		global controller
		controller=self
		self.refreshTimer=QtCore.QTimer()
		self.refreshTimer.timeout.connect(self.refreshEvent)
		self.refreshTimer.start(200)
		self.iterPerSecTimeout=1 # how often to recompute the number of iterations per second
		self.iterTimes,self.iterValues,self.iterPerSec=[],[],0 # arrays to keep track of the simulation speed
		self.dtEditUpdate=True # to avoid updating while being edited
		# show off with this one as well now
	def addPreprocessors(self):
		for c in sudodem.system.childClasses('FileGenerator'):
			self.generatorCombo.addItem(c)
	def addRenderers(self):
		self.displayCombo.addItem('OpenGLRenderer'); afterSep=1
		for bc in ('GlShapeFunctor','GlStateFunctor','GlBoundFunctor','GlIGeomFunctor','GlIPhysFunctor'):
			if afterSep>0: self.displayCombo.insertSeparator(10000); afterSep=0
			for c in sudodem.system.childClasses(bc) | set([bc]):
				inst=eval(c+'()');
				if len(set(inst.dict().keys())-set(['label']))>0:
					self.displayCombo.addItem(c); afterSep+=1
	def inspectSlot(self):
		self.inspector=SimulationInspector(parent=None)
		self.inspector.show()
	def setTabActive(self,what):
		if what=='simulation': ix=0
		elif what=='display': ix=1
		elif what=='generator': ix=2
		elif what=='python': ix=3
		else: raise ValueErorr("No such tab: "+what)
		self.controllerTabs.setCurrentIndex(ix)
	def generatorComboSlot(self,genStr):
		"update generator parameters when a new one is selected"
		gen=eval(str(genStr)+'()')
		self.generator=gen
		se=SerializableEditor(gen,parent=self.generatorArea,showType=True)
		self.generatorArea.setWidget(se)
	def pythonComboSlot(self,cmd):
		try:
			code=compile(str(cmd),'<UI entry>','exec')
			exec(code, globals())
		except:
			import traceback
			traceback.print_exc()
	def generateSlot(self):
		filename=str(self.generatorFilenameEdit.text())
		if self.generatorMemoryCheck.isChecked():
			filename=':memory:'+filename
			print('BUG: Saving to memory slots freezes SudoDEM (cause unknown). Cross fingers.')
		#print 'Will save to ',filename
		self.generator.generate(filename)
		if self.generatorAutoCheck:
			O.load(filename)
			self.setTabActive('simulation')
			if len(views())==0:
				v=View(); v.center()
	def displayComboSlot(self,dispStr):
		ser=(self.renderer if dispStr=='OpenGLRenderer' else eval(str(dispStr)+'()'))
		path='sudodem.qt.Renderer()' if dispStr=='OpenGLRenderer' else dispStr
		se=SerializableEditor(ser,parent=self.displayArea,ignoredAttrs=set(['label']),showType=True,path=path)
		self.displayArea.setWidget(se)
	def loadSlot(self):
		f=QFileDialog.getOpenFileName(self,'Load simulation','','SudoDEM simulations (*.xml *.xml.bz2 *.xml.gz *.sudodem *.sudodem.gz *.sudodem.bz2);; *.*')
		f=str(f)
		if not f: return # cancelled
		self.deactivateControls()
		O.load(f)
	def saveSlot(self):
		f=QFileDialog.getSaveFileName(self,'Save simulation','','SudoDEM simulations (*.xml *.xml.bz2 *.xml.gz *.sudodem *.sudodem.gz *.sudodem.bz2);; *.*')
		f=str(f)
		if not f: return # cancelled
		O.save(f)
	def reloadSlot(self):
		self.deactivateControls()
		from sudodem import plot
		plot.splitData()
		O.reload()
	def dtFixedSlot(self):
		O.dt=O.dt
		O.dynDt=False
	def dtDynSlot(self):
		O.dt=-O.dt
	def dtEditNoupdateSlot(self):
		self.dtEditUpdate=False
	def dtEditedSlot(self):
		try:
			t=float(self.dtEdit.text())
			O.dt=t
		except ValueError: pass
		self.dtEdit.setText(str(O.dt))
		self.dtEditUpdate=True
	def playSlot(self):	O.run()
	def pauseSlot(self): O.pause()
	def stepSlot(self):  O.step()
	def subStepSlot(self,value): O.subStepping=bool(value)
	def show3dSlot(self, show):
		vv=views()
		assert(len(vv) in (0,1))
		if show:
			if len(vv)==0: View()
		else:
			if len(vv)>0: vv[0].close()
	def setReferenceSlot(self):
		# sets reference periodic cell as well
		utils.setRefSe3()
	def centerSlot(self):
		for v in views(): v.center()
	def setViewAxes(self,dir,up):
		try:
			v=views()[0]
			v.viewDir=dir
			v.upVector=up
			v.center()
		except IndexError: pass
	def xyzSlot(self): self.setViewAxes((0,0,-1),(0,1,0))
	def yzxSlot(self): self.setViewAxes((-1,0,0),(0,0,1))
	def zxySlot(self): self.setViewAxes((0,-1,0),(1,0,0))
	def refreshEvent(self):
		self.refreshValues()
		self.activateControls()
	def deactivateControls(self):
		self.realTimeLabel.setText('')
		self.virtTimeLabel.setText('')
		self.iterLabel.setText('')
		self.fileLabel.setText('<i>[loading]</i>')
		self.playButton.setEnabled(False)
		self.pauseButton.setEnabled(False)
		self.stepButton.setEnabled(False)
		self.subStepCheckbox.setEnabled(False)
		self.reloadButton.setEnabled(False)
		self.dtFixedRadio.setEnabled(False)
		self.dtDynRadio.setEnabled(False)
		self.dtEdit.setEnabled(False)
		self.dtEdit.setText('')
		self.dtEditUpdate=True
	def activateControls(self):
		hasSim=len(O.engines)>0
		running=O.running
		if hasSim:
			self.playButton.setEnabled(not running)
			self.pauseButton.setEnabled(running)
			self.reloadButton.setEnabled(O.filename is not None)
			self.stepButton.setEnabled(not running)
			self.subStepCheckbox.setEnabled(not running)
		else:
			self.playButton.setEnabled(False)
			self.pauseButton.setEnabled(False)
			self.reloadButton.setEnabled(False)
			self.stepButton.setEnabled(False)
			self.subStepCheckbox.setEnabled(False)
		self.dtFixedRadio.setEnabled(True)
		self.dtDynRadio.setEnabled(O.dynDtAvailable)
		dynDt=O.dynDt
		self.dtFixedRadio.setChecked(not dynDt)
		self.dtDynRadio.setChecked(dynDt)
		if dynDt or self.dtEditUpdate:
			self.dtEdit.setText(str(O.dt))
		if dynDt: self.dtEditUpdate=True
		self.dtEdit.setEnabled(not dynDt)
		fn=O.filename
		self.fileLabel.setText(fn if fn else '<i>[no file]</i>')

	def refreshValues(self):
		rt=int(O.realtime); t=O.time; iter=O.iter;
		assert(len(self.iterTimes)==len(self.iterValues))
		if len(self.iterTimes)==0: self.iterTimes.append(rt); self.iterValues.append(iter); self.iterPerSec=0 # update always for the first time
		elif rt-self.iterTimes[-1]>self.iterPerSecTimeout: # update after a timeout
			if len(self.iterTimes)==1: self.iterTimes.append(self.iterTimes[0]); self.iterValues.append(self.iterValues[0]) # 2 values, first one is bogus
			self.iterTimes[0]=self.iterTimes[1]; self.iterValues[0]=self.iterValues[1]
			self.iterTimes[1]=rt; self.iterValues[1]=iter;
			self.iterPerSec=(self.iterValues[-1]-self.iterValues[-2])/(self.iterTimes[-1]-self.iterTimes[-2])
		if not O.running: self.iterPerSec=0
		stopAtIter=O.stopAtIter
		subStepInfo=''
		if O.subStepping:
			subStep=O.subStep
			if subStep==-1: subStepInfo=u'→ <i>prologue</i>'
			elif subStep>=0 and subStep<len(O.engines):
				e=O.engines[subStep]; subStepInfo=u'→ %s'%(e.label if e.label else e.__class__.__name__)
			elif subStep==len(O.engines): subStepInfo=u'→ <i>epilogue</i>'
			else: raise RuntimeError("Invalid O.subStep value %d, should be ∈{-1,…,len(o.engines)}"%subStep)
			subStepInfo="<br><small>sub %d/%d [%s]</small>"%(subStep,len(O.engines),subStepInfo)
		self.subStepCheckbox.setChecked(O.subStepping) # might have been changed async
		if stopAtIter<=iter:
			self.realTimeLabel.setText('%02d:%02d:%02d'%(rt//3600,(rt%3600)//60,rt%60))
			self.iterLabel.setText('#%ld, %.1f/s %s'%(iter,self.iterPerSec,subStepInfo))
		else:
			e=int((stopAtIter-iter)*self.iterPerSec)
			self.realTimeLabel.setText('%02d:%02d:%02d (ETA %02d:%02d:%02d)'%(rt//3600,rt//60,rt%60,e//3600,e//60,e%60))
			self.iterLabel.setText('#%ld / %ld, %.1f/s %s'%(O.iter,stopAtIter,self.iterPerSec,subStepInfo))
		if t!=float('inf'):
			s=int(t); ms=int(t*1000)%1000; us=int(t*1000000)%1000; ns=int(t*1000000000)%1000
			self.virtTimeLabel.setText(u'%03ds%03dm%03dμ%03dn'%(s,ms,us,ns))
		else: self.virtTimeLabel.setText(u'[ ∞ ] ?!')
		self.show3dButton.setChecked(len(views())>0)

def Generator():
	global controller
	if not controller: controller=ControllerClass();
	controller.show(); controller.raise_()
	controller.setTabActive('generator')
def Controller():
	global controller
	if not controller: controller=ControllerClass();
	controller.show(); controller.raise_()
	controller.setTabActive('simulation')
def Inspector():
	global controller
	if not controller: controller=ControllerClass();
	controller.inspectSlot()

#if __name__=='__main__':
#	from PyQt4 import QtGui
#	import sys
#	qapp=QtGui.QApplication(sys.argv)
#	c=Controller().show()
#	qapp.exec_()
