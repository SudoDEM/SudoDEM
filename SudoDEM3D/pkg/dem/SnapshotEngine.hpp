#pragma once

#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/PeriodicEngines.hpp>
#include<sudodem/gui/qt4/OpenGLManager.hpp>


class SnapshotEngine: public PeriodicEngine {
	public:
	virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(SnapshotEngine,PeriodicEngine,"Periodically save snapshots of GLView(s) as .png files. Files are named *fileBase* + *counter* + '.png' (counter is left-padded by 0s, i.e. snap00004.png).",
		((string,format,"PNG",,"Format of snapshots (one of JPEG, PNG, EPS, PS, PPM, BMP) `QGLViewer documentation <http://www.libqglviewer.com/refManual/classQGLViewer.html#abbb1add55632dced395e2f1b78ef491c>`_. File extension will be lowercased *format*. Validity of format is not checked."))
		((string,fileBase,"",,"Basename for snapshots"))
		((int,counter,0,,"Number that will be appended to fileBase when the next snapshot is saved (incremented at every save). |yupdate|"))
		((bool,ignoreErrors,true,,"Only report errors instead of throwing exceptions, in case of timeouts."))
		((vector<string>,snapshots,,,"Files that have been created so far"))
		((int,msecSleep,0,,"number of msec to sleep after snapshot (to prevent 3d hw problems) [ms]"))
		((Real,deadTimeout,3,,"Timeout for 3d operations (opening new view, saving snapshot); after timing out, throw exception (or only report error if *ignoreErrors*) and make myself :yref:`dead<Engine.dead>`. [s]"))
		((string,plot,,,"Name of field in :yref:`sudodem.plot.imgData` to which taken snapshots will be appended automatically."))
	);
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(SnapshotEngine);
