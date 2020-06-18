// Copyright (C) 2004 by Olivier Galizzi, Janek Kozicki                  *
// © 2008 Václav Šmilauer
#pragma once

#ifndef Q_MOC_RUN
#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/OpenGLRenderer.hpp>
#include<sudodem/pkg/common/PeriodicEngines.hpp>
#include<boost/date_time/posix_time/posix_time.hpp>
#endif

#include<sudodem/lib/QGLViewer/qglviewer.h>
#include<sudodem/lib/QGLViewer/manipulatedFrame.h>
#include<sudodem/lib/QGLViewer/constraint.h>

#include <qglobal.h>

using std::setw;
using std::setfill;
using std::setprecision;

/*! Class handling user interaction with the openGL rendering of simulation.
 *
 * Clipping planes:
 * ================
 *
 * Clipping plane is manipulated after hitting F1, F2, .... To end the manipulation, press Escape.
 *
 * Keystrokes during clipping plane manipulation:
 * * space activates/deactives the clipping plane
 * * x,y,z aligns the plane with yz, xz, xy planes
 * * left-double-click aligns the plane with world coordinates system (canonical planes + 45˚ interpositions)
 * * 1,2,... will align the current plane with #1, #2, ... (same orientation)
 * * r reverses the plane (normal*=-1)a
 *
 * Keystrokes that work regardless of whether a clipping plane is being manipulated:
 * * Alt-1,Alt-2,... adds/removes the respective plane to bound group:
 * 	mutual positions+orientations of planes in the group are maintained when one of those planes is manipulated
 *
 * Clip plane number is 3; change SUDODEM_RENDERER_NUM_CLIP_PLANE, complete switches "|| ..." in keyPressEvent
 * and recompile to have more.
 */
class GLViewer : public QGLViewer
{
	Q_OBJECT

	friend class QGLThread;
	protected:
		shared_ptr<OpenGLRenderer> renderer;

	private :

		bool			isMoving;
		bool			wasDynamic;
		float			cut_plane;
		int			cut_plane_delta;
		bool			gridSubdivide;
		long			last;
		int manipulatedClipPlane;
		set<int> boundClipPlanes;
		shared_ptr<qglviewer::LocalConstraint> xyPlaneConstraint;
		string strBoundGroup(){string ret;FOREACH(int i, boundClipPlanes) ret+=" "+boost::lexical_cast<string>(i+1);return ret;}
		boost::posix_time::ptime last_user_event;

     public:
		//virtual void updateGL(void);

		const int viewId;

		void centerMedianQuartile();
		int 	drawGrid;
		bool 	drawScale;
		int timeDispMask;
		enum{TIME_REAL=1,TIME_VIRT=2,TIME_ITER=4};

		GLViewer(int viewId, const shared_ptr<OpenGLRenderer>& renderer, QGLWidget* shareWidget=0);
		virtual ~GLViewer();
		#if 0
			virtual void paintGL();
		#endif
		virtual void draw();
		virtual void drawWithNames();
		void displayMessage(const std::string& s){ QGLViewer::displayMessage(QString(s.c_str()));}
		void centerScene();
		void centerPeriodic();
		void mouseMovesCamera();
		void mouseMovesManipulatedFrame(qglviewer::Constraint* c=NULL);
		void resetManipulation();
		bool isManipulating();
		void startClipPlaneManipulation(int planeNo);
		//! get QGLViewer state as string (XML); QGLViewer normally only supports saving state to file.
		string getState();
		//! set QGLViewer state from string (XML); QGLVIewer normally only supports loading state from file.
		void setState(string);
		//! Load display parameters (QGLViewer and OpenGLRenderer) from Scene::dispParams[n] and use them
		void useDisplayParameters(size_t n);
		//! Save display parameters (QGOViewer and OpenGLRenderer) to Scene::dispParams[n]
		void saveDisplayParameters(size_t n);
		//! Get radius of the part of scene that fits the current view
		float displayedSceneRadius();
		//! Get center of the part of scene that fits the current view
		qglviewer::Vec displayedSceneCenter();

		//! Adds our attributes to the QGLViewer state that can be saved
		QDomElement domElement(const QString& name, QDomDocument& document) const;
		//! Adds our attributes to the QGLViewer state that can be restored
		void initFromDOMElement(const QDomElement& element);

		// if defined, snapshot will be saved to this file right after being drawn and the string will be reset.
		// this way the caller will be notified of the frame being saved successfully.
		string nextFrameSnapshotFilename;
		#ifdef SUDODEM_GL2PS
			// output stream for gl2ps; initialized as needed
			FILE* gl2psStream;
		#endif

		boost::posix_time::ptime getLastUserEvent();


		DECLARE_LOGGER;
	protected :
		virtual void keyPressEvent(QKeyEvent *e);
		virtual void postDraw();
		// overridden in the player that doesn't get time from system clock but from the db
		virtual string getRealTimeString();
		virtual void closeEvent(QCloseEvent *e);
		virtual void postSelection(const QPoint& point);
		virtual void endSelection(const QPoint &point);
		virtual void mouseDoubleClickEvent(QMouseEvent *e);
		virtual void wheelEvent(QWheelEvent* e);
		virtual void mouseMoveEvent(QMouseEvent *e);
		virtual void mousePressEvent(QMouseEvent *e);
};

/*! Get unconditional lock on a GL view.

Use if you need to manipulate GL context in some way.
The ctor doesn't return until the lock has been acquired
and the lock is released when the GLLock object is desctructed;
*/
class GLLock: public boost::try_mutex::scoped_lock{
	GLViewer* glv;
	public:
		GLLock(GLViewer* _glv);
		~GLLock();
};


class SudoDEMCamera : public qglviewer::Camera
{
	Q_OBJECT
	private:
		float cuttingDistance;
        public :
		SudoDEMCamera():cuttingDistance(0){};
		 float zNear() const;
		virtual void setCuttingDistance(float s){cuttingDistance=s;};
};



