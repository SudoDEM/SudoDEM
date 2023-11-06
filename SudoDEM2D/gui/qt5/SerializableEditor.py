# encoding: utf-8
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *

import re,itertools
import logging
logging.trace=logging.debug
logging.basicConfig(level=logging.INFO)

from sudodem import *
import sudodem.qt

try:
  from minieigen import *
except ImportError:
  from miniEigen import *

seqSerializableShowType=True # show type headings in serializable sequences (takes vertical space, but makes the type hyperlinked)

# BUG: cursor is moved to the beginnign of the input field even if it has focus
#
# checking for focus seems to return True always and cursor is never moved
#
# the 'True or' part effectively disables the condition (so that the cursor is moved always), but it might be fixed in the future somehow
#
# if True or w.hasFocus(): w.home(False)
#
#

def makeWrapperHref(text,className,attr=None,static=False):
	"""Create clickable HTML hyperlink to a SudoDEM class or its attribute.

	:param className: name of the class to link to.
	:param attr: attribute to link to. If given, must exist directly in given *className*; if not given or empty, link to the class itself is created and *attr* is ignored.
	:return: HTML with the hyperref.
	"""
	if not static: return '<a href="%s#sudodem.wrapper.%s%s">%s</a>'%(sudodem.qt.sphinxDocWrapperPage,className,(('.'+attr) if attr else ''),text)
	else:          return '<a href="%s#ystaticattr-%s.%s">%s</a>'%(sudodem.qt.sphinxDocWrapperPage,className,attr,text)

def serializableHref(ser,attr=None,text=None):
	"""Return HTML href to a *ser* optionally to the attribute *attr*.
	The class hierarchy is crawled upwards to find out in which parent class is *attr* defined,
	so that the href target is a valid link. In that case, only single inheritace is assumed and
	the first class from the top defining *attr* is used.

	:param ser: object of class deriving from :yref:`Serializable`, or string; if string, *attr* must be empty.
	:param attr: name of the attribute to link to; if empty, linke to the class itself is created.
	:param text: visible text of the hyperlink; if not given, either class name or attribute name without class name (when *attr* is not given) is used.

	:returns: HTML with the hyperref.
	"""
	# klass is a class name given as string
	if isinstance(ser,str):
		if attr: raise InvalidArgument("When *ser* is a string, *attr* must be empty (only class link can be created)")
		return makeWrapperHref(text if text else ser,ser)
	# klass is a type object
	if attr:
		klass=ser.__class__
		while attr in dir(klass.__bases__[0]): klass=klass.__bases__[0]
		if not text: text=attr
	else:
		klass=ser.__class__
		if not text: text=klass.__name__
	return makeWrapperHref(text,klass.__name__,attr,static=(attr and getattr(klass,attr)==getattr(ser,attr)))

class AttrEditor():
	"""Abstract base class handing some aspects common to all attribute editors.
	Holds exacly one attribute which is updated whenever it changes."""
	def __init__(self,getter=None,setter=None):
		self.getter,self.setter=getter,setter
		self.hot,self.focused=False,False
		self.widget=None
	def refresh(self): pass
	def update(self): pass
	def isHot(self,hot=True):
		"Called when the widget gets focus; mark it hot, change colors etc."
		if hot==self.hot: return
		self.hot=hot
		if hot: self.setStyleSheet('QWidget { background: red }')
		else: self.setStyleSheet('QWidget { background: none }')
	def sizeHint(self): return QSize(150,12)
	def trySetter(self,val):
		try: self.setter(val)
		except AttributeError: self.setEnabled(False)
		self.isHot(False)

class AttrEditor_Bool(AttrEditor,QFrame):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QFrame.__init__(self,parent)
		self.checkBox=QCheckBox(self)
		lay=QVBoxLayout(self); lay.setSpacing(0); lay.setContentsMargins(0, 0, 0, 0); lay.addStretch(1); lay.addWidget(self.checkBox); lay.addStretch(1)
		self.checkBox.clicked.connect(self.update)
	def refresh(self): self.checkBox.setChecked(self.getter())
	def update(self): self.trySetter(self.checkBox.isChecked())

class AttrEditor_Int(AttrEditor,QSpinBox):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QSpinBox.__init__(self,parent)
		self.setRange(int(-1e9),int(1e9)); self.setSingleStep(1);
		self.valueChanged.connect(self.update)
	def refresh(self): self.setValue(self.getter())
	def update(self):  self.trySetter(self.value())

class AttrEditor_Str(AttrEditor,QLineEdit):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QLineEdit.__init__(self,parent)
		self.textEdited.connect(self.isHot)
		self.selectionChanged.connect(self.isHot)
		self.editingFinished.connect(self.update)
	def refresh(self): self.setText(self.getter())
	def update(self):  self.trySetter(str(self.text()))

class AttrEditor_Float(AttrEditor,QLineEdit):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QLineEdit.__init__(self,parent)
		self.textEdited.connect(self.isHot)
		self.selectionChanged.connect(self.isHot)
		self.editingFinished.connect(self.update)
	def refresh(self):
		self.setText(str(self.getter()));
		if True or not self.hasFocus(): self.home(False)
	def update(self):
		try: self.trySetter(float(self.text()))
		except ValueError: self.refresh()

class AttrEditor_Quaternion(AttrEditor,QFrame):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QFrame.__init__(self,parent)
		self.grid=QHBoxLayout(self); self.grid.setSpacing(0); self.grid.setContentsMargins(0, 0, 0, 0)
		for i in range(4):
			if i==3:
				f=QFrame(self); f.setFrameShape(QFrame.VLine); f.setFrameShadow(QFrame.Sunken); f.setFixedWidth(4) # add vertical divider (axis | angle)
				self.grid.addWidget(f)
			w=QLineEdit('')
			self.grid.addWidget(w);
			w.textEdited.connect(self.isHot)
			w.selectionChanged.connect(self.isHot)
			w.editingFinished.connect(self.update)
	def refresh(self):
		val=self.getter(); axis,angle=val.toAxisAngle()
		for i in (0,1,2,4):
			w=self.grid.itemAt(i).widget(); w.setText(str(axis[i] if i<3 else angle));
			if True or not w.hasFocus(): w.home(False)
	def update(self):
		try:
			x=[float((self.grid.itemAt(i).widget().text())) for i in (0,1,2,4)]
		except ValueError: self.refresh()
		q=Quaternion(Vector3(x[0],x[1],x[2]),x[3]); q.normalize() # from axis-angle
		self.trySetter(q)
	def setFocus(self): self.grid.itemAt(0).widget().setFocus()

class AttrEditor_Se3(AttrEditor,QFrame):
	def __init__(self,parent,getter,setter):
		AttrEditor.__init__(self,getter,setter)
		QFrame.__init__(self,parent)
		self.grid=QGridLayout(self); self.grid.setSpacing(0); self.grid.setContentsMargins(0, 0, 0, 0)
		for row,col in itertools.product(range(2),range(5)): # one additional column for vertical line in quaternion
			if (row,col)==(0,3): continue
			if (row,col)==(0,4): self.grid.addWidget(QLabel(u'←<i>pos</i> ↙<i>ori</i>',self),row,col); continue
			if (row,col)==(1,3):
				f=QFrame(self); f.setFrameShape(QFrame.VLine); f.setFrameShadow(QFrame.Sunken); f.setFixedWidth(4); self.grid.addWidget(f,row,col); continue
			w=QLineEdit('')
			self.grid.addWidget(w,row,col);
			w.textEdited.connect(self.isHot)
			w.selectionChanged.connect(self.isHot)
			w.editingFinished.connect(self.update)
	def refresh(self):
		pos,ori=self.getter(); axis,angle=ori.toAxisAngle()
		for i in (0,1,2,4):
			w=self.grid.itemAtPosition(1,i).widget(); w.setText(str(axis[i] if i<3 else angle));
			if True or not w.hasFocus(): w.home(False)
		for i in (0,1,2):
			w=self.grid.itemAtPosition(0,i).widget(); w.setText(str(pos[i]));
			if True or not w.hasFocus(): w.home(False)
	def update(self):
		try:
			q=[float((self.grid.itemAtPosition(1,i).widget().text())) for i in (0,1,2,4)]
			v=[float((self.grid.itemAtPosition(0,i).widget().text())) for i in (0,1,2)]
		except ValueError: self.refresh()
		qq=Quaternion(Vector3(q[0],q[1],q[2]),q[3]); qq.normalize() # from axis-angle
		self.trySetter((v,qq))
	def setFocus(self): self.grid.itemAtPosition(0,0).widget().setFocus()

class AttrEditor_MatrixX(AttrEditor,QFrame):
	def __init__(self,parent,getter,setter,rows,cols,idxConverter):
		'idxConverter converts row,col tuple to either (row,col), (col) etc depending on what access is used for []'
		AttrEditor.__init__(self,getter,setter)
		QFrame.__init__(self,parent)
		self.rows,self.cols=rows,cols
		self.idxConverter=idxConverter
		self.setContentsMargins(0,0,0,0)
		val=self.getter()
		self.grid=QGridLayout(self); self.grid.setSpacing(0); self.grid.setContentsMargins(0, 0, 0, 0)
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			w=QLineEdit('')
			self.grid.addWidget(w,row,col);
			w.textEdited.connect(self.isHot)
			w.selectionChanged.connect(self.isHot)
			w.editingFinished.connect(self.update)
	def refresh(self):
		val=self.getter()
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			w=self.grid.itemAtPosition(row,col).widget()
			w.setText(str(val[self.idxConverter(row,col)]))
			if True or not w.hasFocus: w.home(False) # make the left-most part visible, if the text is wider than the widget
	def update(self):
		try:
			val=self.getter()
			for row,col in itertools.product(range(self.rows),range(self.cols)):
				w=self.grid.itemAtPosition(row,col).widget()
				if w.isModified(): val[self.idxConverter(row,col)]=float(w.text())
			logging.debug('setting'+str(val))
			self.trySetter(val)
		except ValueError: self.refresh()
	def setFocus(self): self.grid.itemAtPosition(0,0).widget().setFocus()

class AttrEditor_MatrixXi(AttrEditor,QFrame):
	def __init__(self,parent,getter,setter,rows,cols,idxConverter):
		'idxConverter converts row,col tuple to either (row,col), (col) etc depending on what access is used for []'
		AttrEditor.__init__(self,getter,setter)
		QFrame.__init__(self,parent)
		self.rows,self.cols=rows,cols
		self.idxConverter=idxConverter
		self.setContentsMargins(0,0,0,0)
		self.grid=QGridLayout(self); self.grid.setSpacing(0); self.grid.setContentsMargins(0, 0, 0, 0)
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			w=QSpinBox()
			w.setRange(int(-1e9),int(1e9)); w.setSingleStep(1);
			self.grid.addWidget(w,row,col);
		self.refresh() # refresh before connecting signals!
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			self.grid.itemAtPosition(row,col).widget().valueChanged.connect(self.update)
	def refresh(self):
		val=self.getter()
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			w=self.grid.itemAtPosition(row,col).widget().setValue(val[self.idxConverter(row,col)])
	def update(self):
		val=self.getter(); modified=False
		for row,col in itertools.product(range(self.rows),range(self.cols)):
			w=self.grid.itemAtPosition(row,col).widget()
			if w.value()!=val[self.idxConverter(row,col)]:
				modified=True; val[self.idxConverter(row,col)]=w.value()
		if not modified: return
		logging.debug('setting'+str(val))
		self.trySetter(val)
	def setFocus(self): self.grid.itemAtPosition(0,0).widget().setFocus()

class AttrEditor_Vector6i(AttrEditor_MatrixXi):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixXi.__init__(self,parent,getter,setter,1,6,lambda r,c:c)
class AttrEditor_Vector3i(AttrEditor_MatrixXi):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixXi.__init__(self,parent,getter,setter,1,3,lambda r,c:c)
class AttrEditor_Vector2i(AttrEditor_MatrixXi):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixXi.__init__(self,parent,getter,setter,1,2,lambda r,c:c)

class AttrEditor_Vector6(AttrEditor_MatrixX):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixX.__init__(self,parent,getter,setter,1,6,lambda r,c:c)
class AttrEditor_Vector3(AttrEditor_MatrixX):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixX.__init__(self,parent,getter,setter,1,3,lambda r,c:c)
class AttrEditor_Vector2(AttrEditor_MatrixX):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixX.__init__(self,parent,getter,setter,1,2,lambda r,c:c)
class AttrEditor_Matrix3(AttrEditor_MatrixX):
	def __init__(self,parent,getter,setter):
		AttrEditor_MatrixX.__init__(self,parent,getter,setter,3,3,lambda r,c:(r,c))

class Se3FakeType: pass

_fundamentalEditorMap={bool:AttrEditor_Bool,str:AttrEditor_Str,int:AttrEditor_Int,float:AttrEditor_Float,Quaternion:AttrEditor_Quaternion,Vector2:AttrEditor_Vector2,Vector3:AttrEditor_Vector3,Vector6:AttrEditor_Vector6,Matrix3:AttrEditor_Matrix3,Vector6i:AttrEditor_Vector6i,Vector3i:AttrEditor_Vector3i,Vector2i:AttrEditor_Vector2i,Se3FakeType:AttrEditor_Se3}
_fundamentalInitValues={bool:True,str:'',int:0,float:0.0,Quaternion:Quaternion((0,1,0),0.0),Vector3:Vector3.Zero,Matrix3:Matrix3.Zero,Vector6:Vector6.Zero,Vector6i:Vector6i.Zero,Vector3i:Vector3i.Zero,Vector2i:Vector2i.Zero,Vector2:Vector2.Zero,Se3FakeType:(Vector3.Zero,Quaternion((0,1,0),0.0))}

class SerQLabel(QLabel):
	def __init__(self,parent,label,tooltip,path):
		QLabel.__init__(self,parent)
		self.path=path
		self.setText(label)
		if tooltip or path: self.setToolTip(('<b>'+path+'</b><br>' if self.path else '')+(tooltip if tooltip else ''))
		self.linkActivated.connect(sudodem.qt.openUrl)
	def mousePressEvent(self,event):
		if event.button()!=Qt.MidButton:
			event.ignore(); return
		# middle button clicked, paste pasteText to clipboard
		cb=QApplication.clipboard()
		cb.setText(self.path,mode=QClipboard.Clipboard)
		cb.setText(self.path,mode=QClipboard.Selection) # X11 global selection buffer
		event.accept()

class SerializableEditor(QFrame):
	"Class displaying and modifying serializable attributes of a sudodem object."
	import collections
	import logging
	# each attribute has one entry associated with itself
	class EntryData:
		def __init__(self,name,T,flags=0):
			self.name,self.T,self.flags=name,T,flags
			self.lineNo,self.widget=None,None
	def __init__(self,ser,parent=None,ignoredAttrs=set(),showType=False,path=None):
		"Construct window, *ser* is the object we want to show."
		QtWidgets.QFrame.__init__(self,parent)
		self.ser=ser
		self.path=(ser.label if (hasattr(ser,'label') and ser.label) else path)
		self.showType=showType
		self.hot=False
		self.entries=[]
		self.ignoredAttrs=ignoredAttrs
		logging.debug('New Serializable of type %s'%ser.__class__.__name__)
		self.setWindowTitle(str(ser))
		self.mkWidgets()
		self.refreshTimer=QTimer(self)
		self.refreshTimer.timeout.connect(self.refreshEvent)
		self.refreshTimer.start(500)
	def getListTypeFromDocstring(self,attr):
		"Guess type of array by scanning docstring for :yattrtype: and parsing its argument; ugly, but works."
		doc=getattr(self.ser.__class__,attr).__doc__
		if doc==None:
			logging.error("Attribute %s has no docstring."%attr)
			return None
		m=re.search(r':yattrtype:`([^`]*)`',doc)
		if not m:
			logging.error("Attribute %s does not contain :yattrtype:`....` (docstring is '%s'"%(attr,doc))
			return None
		cxxT=m.group(1)
		logging.debug('Got type "%s" from :yattrtype:'%cxxT)
		def vecTest(T,cxxT):
			#regexp=r'^\s*(std\s*::)?\s*vector\s*<\s*(std\s*::)?\s*('+T+r')\s*>\s*$'
			regexp=r'^\s*(std\s*::)?\s*vector\s*<\s*(shared_ptr\s*<\s*)?\s*(std\s*::)?\s*('+T+r')(\s*>)?\s*>\s*$'
			m=re.match(regexp,cxxT)
			return m
		vecMap={
			'bool':bool,'int':int,'long':int,'Body::id_t':int,'size_t':int,
			'Real':float,'float':float,'double':float,
			'Vector6r':Vector6,'Vector6i':Vector6i,'Vector3i':Vector3i,'Vector2r':Vector2,'Vector2i':Vector2i,
			'Vector3r':Vector3,'Matrix3r':Matrix3,'Se3r':Se3FakeType,
			'string':str,
			#'BodyCallback':BodyCallback,
			#'IntrCallback':IntrCallback,'BoundFunctor':BoundFunctor,'IGeomFunctor':IGeomFunctor,'IPhysFunctor':IPhysFunctor,'LawFunctor':LawFunctor,'KinematicEngine':KinematicEngine,
            'IntrCallback':IntrCallback,'BoundFunctor':BoundFunctor,'IGeomFunctor':IGeomFunctor,'IPhysFunctor':IPhysFunctor,'LawFunctor':LawFunctor,
			'GlShapeFunctor':GlShapeFunctor,'GlStateFunctor':GlStateFunctor,'GlIGeomFunctor':GlIGeomFunctor,'GlIPhysFunctor':GlIPhysFunctor,'GlBoundFunctor':GlBoundFunctor,'GlExtraDrawer':GlExtraDrawer
		}
		for T,ret in vecMap.items():
			if vecTest(T,cxxT):
				logging.debug("Got type %s from cxx type %s"%(repr(ret),cxxT))
				return (ret,)
		logging.error("Unable to guess python type from cxx type '%s'"%cxxT)
		return None
	def mkAttrEntries(self):
		if self.ser==None: return
		try:
			d=self.ser.dict()
		except TypeError:
			logging.error('TypeError when getting attributes of '+str(self.ser)+',skipping. ')
			import traceback
			traceback.print_exc()
		attrs=list(self.ser.dict().keys()); attrs.sort()
		for attr in attrs:
			val=getattr(self.ser,attr) # get the value using serattr, as it might be different from what the dictionary provides (e.g. Body.blockedDOFs)
			t=None
			doc=getattr(self.ser.__class__,attr).__doc__;
			if '|yhidden|' in doc: continue
			if attr in self.ignoredAttrs: continue
			if isinstance(val,list):
				t=self.getListTypeFromDocstring(attr)
				if not t and len(val)==0: t=(val[0].__class__,) # 1-tuple is list of the contained type
				#if not t: raise RuntimeError('Unable to guess type of '+str(self.ser)+'.'+attr)
			# hack for Se3, which is returned as (Vector3,Quaternion) in python
			elif isinstance(val,tuple) and len(val)==2 and val[0].__class__==Vector3 and val[1].__class__==Quaternion: t=Se3FakeType
			else: t=val.__class__
			match=re.search(':yattrflags:`\s*([0-9]+)\s*`',doc) # non-empty attribute
			flags=int(match.group(1)) if match else 0
			#logging.debug('Attr %s is of type %s'%(attr,((t[0].__name__,) if isinstance(t,tuple) else t.__name__)))
			self.entries.append(self.EntryData(name=attr,T=t))
	def getDocstring(self,attr=None):
		"If attr is *None*, return docstring of the Serializable itself"
		doc=(getattr(self.ser.__class__,attr).__doc__ if attr else self.ser.__class__.__doc__)
		if not doc: return ''
		doc=re.sub(':y(attrtype|default|attrflags):`[^`]*`','',doc)
		statAttr=re.compile('^.. ystaticattr::.*$',re.MULTILINE|re.DOTALL)
		doc=re.sub(statAttr,'',doc) # static classes have their proper docs at the beginning, discard static memeber docs
		# static: attribute of the type is the same object as attribute of the instance
		# in that case, get docstring from the class documentation by parsing it
		if attr and getattr(self.ser.__class__,attr)==getattr(self.ser,attr): doc=self.getStaticAttrDocstring(attr)
		doc=re.sub(':yref:`([^`]*)`','\\1',doc)
		import textwrap
		wrapper=textwrap.TextWrapper(replace_whitespace=False)
		return wrapper.fill(textwrap.dedent(doc))
	def getStaticAttrDocstring(self,attr):
		ret=''; c=self.ser.__class__
		while hasattr(c,attr) and hasattr(c.__base__,attr): c=c.__base__
		start='.. ystaticattr:: %s.%s('%(c.__name__,attr)
		if start in c.__doc__:
			ll=c.__doc__.split('\n')
			for i in range(len(ll)):
				if ll[i].startswith(start): break
			for i in range(i+1,len(ll)):
				if len(ll[i])>0 and ll[i][0] not in ' \t': break
				ret+=ll[i]
			return ret
		else: return '[no documentation found]'
	def mkWidget(self,entry):
		if not entry.T: return None
		# single fundamental object
		Klass=_fundamentalEditorMap.get(entry.T,None)
		getter,setter=lambda: getattr(self.ser,entry.name), lambda x: setattr(self.ser,entry.name,x)
		if Klass:
			widget=Klass(self,getter=getter,setter=setter)
			widget.setFocusPolicy(Qt.StrongFocus)
			if (entry.flags & AttrFlags.readonly): widget.setEnabled(False)
			return widget
		# sequences
		if entry.T.__class__==tuple:
			assert(len(entry.T)==1) # we don't handle tuples of other lenghts
			# sequence of serializables
			T=entry.T[0]
			if (issubclass(T,Serializable) or T==Serializable):
				widget=SeqSerializable(self,getter,setter,T,path=(self.path+'.'+entry.name if self.path else None),shrink=True)
				return widget
			if (T in _fundamentalEditorMap):
				widget=SeqFundamentalEditor(self,getter,setter,T)
				return widget
			return None
		# a serializable
		if issubclass(entry.T,Serializable) or entry.T==Serializable:
			obj=getattr(self.ser,entry.name)
			if hasattr(obj,'label') and obj.label: path=obj.label
			elif self.path: path=self.path+'.'+entry.name
			else: path=None
			widget=SerializableEditor(getattr(self.ser,entry.name),parent=self,showType=self.showType,path=(self.path+'.'+entry.name if self.path else None))
			widget.setFrameShape(QFrame.Box); widget.setFrameShadow(QFrame.Raised); widget.setLineWidth(1)
			return widget
		return None
	def mkWidgets(self):
		self.mkAttrEntries()
		grid=QFormLayout()
		grid.setContentsMargins(2,2,2,2)
		grid.setVerticalSpacing(0)
		grid.setLabelAlignment(Qt.AlignRight)
		if self.showType:
			lab=SerQLabel(self,makeSerializableLabel(self.ser,addr=True,href=True),tooltip=self.getDocstring(),path=self.path)
			lab.setFrameShape(QFrame.Box); lab.setFrameShadow(QFrame.Sunken); lab.setLineWidth(2); lab.setAlignment(Qt.AlignHCenter); lab.linkActivated.connect(sudodem.qt.openUrl)
			grid.setWidget(0,QFormLayout.SpanningRole,lab)
		for entry in self.entries:
			entry.widget=self.mkWidget(entry)
			objPath=(self.path+'.'+entry.name) if self.path else None
			label=SerQLabel(self,serializableHref(self.ser,entry.name),tooltip=self.getDocstring(entry.name),path=objPath)
			grid.addRow(label,entry.widget if entry.widget else QLabel('<i>unhandled type</i>'))
		self.setLayout(grid)
		self.refreshEvent()
	def refreshEvent(self):
		for e in self.entries:
			if e.widget and not e.widget.hot: e.widget.refresh()
	def refresh(self): pass

def makeSerializableLabel(ser,href=False,addr=True,boldHref=True,num=-1,count=-1):
	ret=u''
	if num>=0:
		if count>=0: ret+=u'%d/%d. '%(num,count)
		else: ret+=u'%d. '%num
	if href: ret+=(u' <b>' if boldHref else u' ')+serializableHref(ser)+(u'</b> ' if boldHref else u' ')
	else: ret+=ser.__class__.__name__+' '
	if hasattr(ser,'label') and ser.label: ret+=u' “'+str(ser.label)+u'”'
	# do not show address if there is a label already
	elif addr:
		import re
		ss=str(ser); m=re.match(u'<(.*) instance at (0x.*)>',ss)
		if m: ret+=m.group(2)
		else: logging.warning(u"Serializable converted to str ('%s') does not contain 'instance at 0x…'"%ss)
	return ret

class SeqSerializableComboBox(QFrame):
	def __init__(self,parent,getter,setter,serType,path=None,shrink=False):
		QFrame.__init__(self,parent)
		self.getter,self.setter,self.serType,self.path,self.shrink=getter,setter,serType,path,shrink
		self.layout=QVBoxLayout(self)
		topLineFrame=QFrame(self)
		topLineLayout=QHBoxLayout(topLineFrame);
		for l in self.layout, topLineLayout: l.setSpacing(0); l.setContentsMargins(0,0,0,0)
		topLineFrame.setLayout(topLineLayout)
		buttons=(self.newButton,self.killButton,self.upButton,self.downButton)=[QPushButton(label,self) for label in (u'☘',u'☠',u'↑',u'↓')]
		buttonSlots=(self.newSlot,self.killSlot,self.upSlot,self.downSlot) # same order as buttons
		for b in buttons: b.setStyleSheet('QPushButton { font-size: 15pt; }'); b.setFixedWidth(30); b.setFixedHeight(30)
		self.combo=QComboBox(self)
		self.combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
		for w in buttons[0:2]+[self.combo,]+buttons[2:4]: topLineLayout.addWidget(w)
		self.layout.addWidget(topLineFrame) # nested layout
		self.scroll=QScrollArea(self); self.scroll.setWidgetResizable(True)
		self.layout.addWidget(self.scroll)
		self.seqEdit=None # currently edited serializable
		self.setLayout(self.layout)
		self.hot=None # API compat with SerializableEditor
		self.setFrameShape(QFrame.Box); self.setFrameShadow(QFrame.Raised); self.setLineWidth(1)
		# signals
		for b,slot in zip(buttons,buttonSlots): b.clicked.connect(slot)
		self.combo.currentIndexChanged.connect(self.comboIndexSlot)
		self.refreshEvent()
		# periodic refresh
		self.refreshTimer=QTimer(self)
		self.refreshTimer.timeout.connect(self.refreshEvent)
		self.refreshTimer.start(1000) # 1s should be enough
		#print 'SeqSerializable path is',self.path
	def comboIndexSlot(self,ix): # different seq item selected
		currSeq=self.getter();
		if len(currSeq)==0: ix=-1
		logging.debug('%s comboIndexSlot len=%d, ix=%d'%(self.serType.__name__,len(currSeq),ix))
		self.downButton.setEnabled(ix<len(currSeq)-1)
		self.upButton.setEnabled(ix>0)
		self.combo.setEnabled(ix>=0)
		if ix>=0:
			ser=currSeq[ix]
			self.seqEdit=SerializableEditor(ser,parent=self,showType=seqSerializableShowType,path=(self.path+'['+str(ix)+']') if self.path else None)
			self.scroll.setWidget(self.seqEdit)
			if self.shrink:
				self.sizeHint=lambda: QSize(100,1000)
				self.scroll.sizeHint=lambda: QSize(100,1000)
				self.sizePolicy().setVerticalPolicy(QSizePolicy.Expanding)
				self.scroll.sizePolicy().setVerticalPolicy(QSizePolicy.Expanding)
				self.setMinimumHeight(min(300,self.seqEdit.height()+self.combo.height()+10))
				self.setMaximumHeight(100000)
				self.scroll.setMaximumHeight(100000)
		else:
			self.scroll.setWidget(QFrame())
			if self.shrink:
				self.setMaximumHeight(self.combo.height()+10);
				self.scroll.setMaximumHeight(0)
	def serLabel(self,ser,i=-1):
		return ('' if i<0 else str(i)+'. ')+str(ser)[1:-1].replace('instance at ','')
	def refreshEvent(self,forceIx=-1):
		currSeq=self.getter()
		comboEnabled=self.combo.isEnabled()
		if comboEnabled and len(currSeq)==0: self.comboIndexSlot(-1) # force refresh, otherwise would not happen from the initially empty state
		ix,cnt=self.combo.currentIndex(),self.combo.count()
		# serializable currently being edited (which can be absent) or the one of which index is forced
		ser=(self.seqEdit.ser if self.seqEdit else None) if forceIx<0 else currSeq[forceIx]
		if comboEnabled and len(currSeq)==cnt and (ix<0 or ser==currSeq[ix]): return
		if not comboEnabled and len(currSeq)==0: return
		logging.debug(self.serType.__name__+' rebuilding list from scratch')
		self.combo.clear()
		if len(currSeq)>0:
			prevIx=-1
			for i,s in enumerate(currSeq):
				self.combo.addItem(makeSerializableLabel(s,num=i,count=len(currSeq),addr=False))
				if s==ser: prevIx=i
			if forceIx>=0: newIx=forceIx # force the index (used from newSlot to make the new element active)
			elif prevIx>=0: newIx=prevIx # if found what was active before, use it
			elif ix>=0: newIx=ix         # otherwise use the previous index (e.g. after deletion)
			else: newIx=0                  # fallback to 0
			logging.debug('%s setting index %d'%(self.serType.__name__,newIx))
			self.combo.setCurrentIndex(newIx)
		else:
			logging.debug('%s EMPTY, setting index 0'%(self.serType.__name__))
			self.combo.setCurrentIndex(-1)
		self.killButton.setEnabled(len(currSeq)>0)
	def newSlot(self):
		dialog=NewSerializableDialog(self,self.serType.__name__)
		if not dialog.exec_(): return # cancelled
		ser=dialog.result()
		ix=self.combo.currentIndex()
		currSeq=self.getter(); currSeq.insert(ix,ser); self.setter(currSeq)
		logging.debug('%s new item created at index %d'%(self.serType.__name__,ix))
		self.refreshEvent(forceIx=ix)
	def killSlot(self):
		ix=self.combo.currentIndex()
		currSeq=self.getter(); del currSeq[ix]; self.setter(currSeq)
		self.refreshEvent()
	def upSlot(self):
		i=self.combo.currentIndex()
		assert(i>0)
		currSeq=self.getter();
		prev,curr=currSeq[i-1:i+1]; currSeq[i-1],currSeq[i]=curr,prev; self.setter(currSeq)
		self.refreshEvent(forceIx=i-1)
	def downSlot(self):
		i=self.combo.currentIndex()
		currSeq=self.getter(); assert(i<len(currSeq)-1);
		curr,nxt=currSeq[i:i+2]; currSeq[i],currSeq[i+1]=nxt,curr; self.setter(currSeq)
		self.refreshEvent(forceIx=i+1)
	def refresh(self): pass # API compat with SerializableEditor

SeqSerializable=SeqSerializableComboBox


class NewFundamentalDialog(QDialog):
	def __init__(self,parent,attrName,typeObj,typeStr):
		QDialog.__init__(self,parent)
		self.setWindowTitle('%s (type %s)'%(attrName,typeStr))
		self.layout=QVBoxLayout(self)
		self.scroll=QScrollArea(self)
		self.scroll.setWidgetResizable(True)
		self.buttons=QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel);
		self.buttons.accepted.connect(self.accept)
		self.buttons.rejected.connect(self.reject)
		self.layout.addWidget(self.scroll)
		self.layout.addWidget(self.buttons)
		self.setWindowModality(Qt.WindowModal)
		class FakeObjClass: pass
		self.fakeObj=FakeObjClass()
		self.attrName=attrName
		Klass=_fundamentalEditorMap.get(typeObj,None)
		initValue=_fundamentalInitValues.get(typeObj,typeObj())
		setattr(self.fakeObj,attrName,initValue)
		if Klass:
			self.widget=Klass(None,self.fakeObj,attrName)
			self.scroll.setWidget(self.widget)
			self.scroll.show()
			self.widget.refresh()
		else: raise RuntimeError("Unable to construct new dialog for type %s"%(typeStr))
	def result(self):
		self.widget.update()
		return getattr(self.fakeObj,self.attrName)

class NewSerializableDialog(QDialog):
	def __init__(self,parent,baseClassName,includeBase=True):
		import sudodem.system
		QDialog.__init__(self,parent)
		self.setWindowTitle('Create new object of type %s'%baseClassName)
		self.layout=QVBoxLayout(self)
		self.combo=QComboBox(self)
		childs=list(sudodem.system.childClasses(baseClassName,includeBase=False)); childs.sort()
		if includeBase:
			self.combo.addItem(baseClassName)
			self.combo.insertSeparator(1000)
		self.combo.addItems(childs)
		self.combo.currentIndexChanged.connect(self.comboSlot)
		self.scroll=QScrollArea(self)
		self.scroll.setWidgetResizable(True)
		self.buttons=QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel);
		self.buttons.accepted.connect(self.accept)
		self.buttons.rejected.connect(self.reject)
		self.layout.addWidget(self.combo)
		self.layout.addWidget(self.scroll)
		self.layout.addWidget(self.buttons)
		self.ser=None
		self.combo.setCurrentIndex(0); self.comboSlot(0)
		self.setWindowModality(Qt.WindowModal)
	def comboSlot(self,index):
		item=str(self.combo.itemText(index))
		self.ser=eval(item+'()')
		self.scroll.setWidget(SerializableEditor(self.ser,self.scroll,showType=True))
		self.scroll.show()
	def result(self): return self.ser
	def sizeHint(self): return QSize(180,400)

class SeqFundamentalEditor(QFrame):
	def __init__(self,parent,getter,setter,itemType):
		QFrame.__init__(self,parent)
		self.getter,self.setter,self.itemType=getter,setter,itemType
		self.layout=QVBoxLayout()
		topLineFrame=QFrame(self); topLineLayout=QHBoxLayout(topLineFrame)
		self.form=QFormLayout()
		self.form.setContentsMargins(0,0,0,0)
		self.form.setVerticalSpacing(0)
		self.form.setLabelAlignment(Qt.AlignLeft)
		self.formFrame=QFrame(self); self.formFrame.setLayout(self.form)
		self.layout.addWidget(self.formFrame)
		self.setLayout(self.layout)
		# SerializableEditor API compat
		self.hot=False
		self.rebuild()
		# periodic refresh
		self.refreshTimer=QTimer(self)
		self.refreshTimer.timeout.connect(self.refreshEvent)
		self.refreshTimer.start(1000) # 1s should be enough
	def contextMenuEvent(self, event):
		index=self.localPositionToIndex(event.pos())
		seq=self.getter()
		if len(seq)==0: index=-1
		field=self.form.itemAt(index,QFormLayout.LabelRole).widget() if index>=0 else None
		menu=QMenu(self)
		actNew,actKill,actUp,actDown=[menu.addAction(name) for name in (u'☘ New',u'☠ Remove',u'↑ Up',u'↓ Down')]
		if index<0: [a.setEnabled(False) for a in [actKill,actUp,actDown]]
		if index==len(seq)-1: actDown.setEnabled(False)
		if index==0: actUp.setEnabled(False)
		if field: field.setStyleSheet('QWidget { background: green }')
		act=menu.exec_(self.mapToGlobal(event.pos()))
		if field: field.setStyleSheet('QWidget { background: none }')
		if not act: return
		if act==actNew: self.newSlot(index)
		elif act==actKill: self.killSlot(index)
		elif act==actUp: self.upSlot(index)
		elif act==actDown: self.downSlot(index)
	def localPositionToIndex(self,pos):
		gp=self.mapToGlobal(pos)
		for row in range(self.form.count()/2):
			w,i=self.form.itemAt(row,QFormLayout.FieldRole),self.form.itemAt(row,QFormLayout.LabelRole)
			for wi in w.widget(),i.widget():
				x0,y0,x1,y1=wi.geometry().getCoords(); globG=QRect(self.mapToGlobal(QPoint(x0,y0)),self.mapToGlobal(QPoint(x1,y1)))
				if globG.contains(gp):
					return row
		return -1
	def newSlot(self,i):
		seq=self.getter();
		seq.insert(i,_fundamentalInitValues.get(self.itemType,self.itemType()))
		self.setter(seq)
		self.rebuild()
	def killSlot(self,i):
		seq=self.getter(); assert(i<len(seq)); del seq[i]; self.setter(seq)
		self.refreshEvent()
	def upSlot(self,i):
		seq=self.getter(); assert(i<len(seq));
		prev,curr=seq[i-1:i+1]; seq[i-1],seq[i]=curr,prev; self.setter(seq)
		self.refreshEvent(forceIx=i-1)
	def downSlot(self,i):
		seq=self.getter(); assert(i<len(seq)-1);
		curr,nxt=seq[i:i+2]; seq[i],seq[i+1]=nxt,curr; self.setter(seq)
		self.refreshEvent(forceIx=i+1)
	def rebuild(self):
		currSeq=self.getter()
		# clear everything
		rows=self.form.count()/2
		for row in range(rows):
			logging.trace('counts',self.form.rowCount(),self.form.count())
			for wi in self.form.itemAt(row,QFormLayout.FieldRole),self.form.itemAt(row,QFormLayout.LabelRole):
				self.form.removeItem(wi)
				logging.trace('deleting widget',wi.widget())
				widget=wi.widget(); widget.hide(); del widget # for some reason, deleting does not make the thing disappear visually; hiding does, however
			logging.trace('counts after ',self.form.rowCount(),self.form.count())
		logging.debug('cleared')
		# add everything
		Klass=_fundamentalEditorMap.get(self.itemType,None)
		if not Klass:
			errMsg=QTextEdit(self)
			errMsg.setReadOnly(True); errMsg.setText("Sorry, editing sequences of %s's is not (yet?) implemented."%(self.itemType.__name__))
			self.form.insertRow(0,'<b>Error</b>',errMsg)
			return
		class ItemGetter():
			def __init__(self,getter,index): self.getter,self.index=getter,index
			def __call__(self): return self.getter()[self.index]
		class ItemSetter():
			def __init__(self,getter,setter,index): self.getter,self.setter,self.index=getter,setter,index
			def __call__(self,val): seq=self.getter(); seq[self.index]=val; self.setter(seq)
		for i,item in enumerate(currSeq):
			widget=Klass(self,ItemGetter(self.getter,i),ItemSetter(self.getter,self.setter,i)) #proxy,'value')
			self.form.insertRow(i,'%d. '%i,widget)
			logging.debug('added item %d %s'%(i,str(widget)))
		if len(currSeq)==0: self.form.insertRow(0,'<i>empty</i>',QLabel('<i>(right-click for menu)</i>'))
		logging.debug('rebuilt, will refresh now')
		self.refreshEvent(dontRebuild=True) # avoid infinite recursion it the length would change meanwhile
	def refreshEvent(self,dontRebuild=False,forceIx=-1):
		currSeq=self.getter()
		if len(currSeq)!=self.form.count()/2: #rowCount():
			if dontRebuild: return # length changed behind our back, just pretend nothing happened and update next time instead
			self.rebuild()
			currSeq=self.getter()
		for i in range(len(currSeq)):
			item=self.form.itemAt(i,QFormLayout.FieldRole)
			logging.trace('got item #%d %s'%(i,str(item.widget())))
			widget=item.widget()
			if not widget.hot:
				widget.refresh()
			if forceIx>=0 and forceIx==i: widget.setFocus()
	def refresh(self): pass # SerializableEditor API
