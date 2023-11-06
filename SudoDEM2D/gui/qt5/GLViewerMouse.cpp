/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2005 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"GLViewer.hpp"
#include"OpenGLManager.hpp"

#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/core/Body.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/DisplayParameters.hpp>
#include<boost/algorithm/string.hpp>
#include<sstream>
#include<iomanip>
#include<boost/algorithm/string/case_conv.hpp>
#include<sudodem/lib/serialization/ObjectIO.hpp>
#include<sudodem/lib/pyutil/gil.hpp>
#include<QGLViewer/manipulatedCameraFrame.h>

#include<QtGui/qevent.h>

#ifdef SUDODEM_GL2PS
	#include<gl2ps.h>
#endif

void GLViewer::mouseMovesCamera(){
  setWheelBinding(Qt::ShiftModifier , FRAME, ZOOM);
  setWheelBinding(Qt::NoModifier, CAMERA, ZOOM);
  //cout<<"mouse moves camera"<<endl;
#if QGLVIEWER_VERSION>=0x020500
  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, SELECT);

  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, FRAME, ZOOM);
  setMouseBinding(Qt::ShiftModifier, Qt::RightButton, FRAME, ZOOM);
  setMouseBinding(Qt::ShiftModifier, Qt::MidButton, FRAME, TRANSLATE);
  setMouseBinding(Qt::ShiftModifier, Qt::RightButton, FRAME, NO_MOUSE_ACTION);//ROTATE

  setMouseBinding(Qt::NoModifier, Qt::RightButton, CAMERA, ZOOM);
  setMouseBinding(Qt::NoModifier, Qt::MidButton, CAMERA, ZOOM);
  setMouseBinding(Qt::NoModifier, Qt::LeftButton, CAMERA, TRANSLATE);//ROTATE
  setMouseBinding(Qt::NoModifier, Qt::RightButton, CAMERA, TRANSLATE);
#else
  setMouseBinding(Qt::SHIFT + Qt::LeftButton, SELECT);
  setMouseBinding(Qt::SHIFT + Qt::LeftButton + Qt::RightButton, FRAME, ZOOM);
  setMouseBinding(Qt::SHIFT + Qt::MidButton, FRAME, TRANSLATE);
  setMouseBinding(Qt::SHIFT + Qt::RightButton, FRAME, ROTATE);
  cout<<"QGLViewer < 2.5"<<endl;
  setMouseBinding(Qt::LeftButton + Qt::RightButton, CAMERA, ZOOM);
  setMouseBinding(Qt::MidButton, CAMERA, ZOOM);
  setMouseBinding(Qt::LeftButton, CAMERA, ROTATE);
  setMouseBinding(Qt::RightButton, CAMERA, TRANSLATE);
  camera()->frame()->setWheelSensitivity(-1.0f);
#endif
};

void GLViewer::mouseMovesManipulatedFrame(qglviewer::Constraint* c){
	setMouseBinding(Qt::LeftButton + Qt::RightButton, FRAME, ZOOM);
	setMouseBinding(Qt::MidButton, FRAME, ZOOM);
	setMouseBinding(Qt::LeftButton, FRAME, TRANSLATE);
	setMouseBinding(Qt::RightButton, FRAME, TRANSLATE);
	setWheelBinding(Qt::NoModifier , FRAME, ZOOM);
	manipulatedFrame()->setConstraint(c);
}


void GLViewer::mouseMoveEvent(QMouseEvent *e){
	last_user_event = boost::posix_time::second_clock::local_time();
	QGLViewer::mouseMoveEvent(e);
}

void GLViewer::mousePressEvent(QMouseEvent *e){
	last_user_event = boost::posix_time::second_clock::local_time();
	QGLViewer::mousePressEvent(e);
}

/* Handle double-click event; if clipping plane is manipulated, align it with the global coordinate system.
 * Otherwise pass the event to QGLViewer to handle it normally.
 *
 * mostly copied over from ManipulatedFrame::mouseDoubleClickEvent
 */
void GLViewer::mouseDoubleClickEvent(QMouseEvent *event){
	last_user_event = boost::posix_time::second_clock::local_time();

	if(manipulatedClipPlane<0) { /* LOG_DEBUG("Double click not on clipping plane"); */ QGLViewer::mouseDoubleClickEvent(event); return; }
#if QT_VERSION >= 0x040000
	if (event->modifiers() == Qt::NoModifier)
#else
	if (event->state() == Qt::NoButton)
#endif
	switch (event->button()){
		case Qt::LeftButton:  manipulatedFrame()->alignWithFrame(NULL,true); LOG_DEBUG("Aligning cutting plane"); break;
		default: break; // avoid warning
	}
}

void GLViewer::wheelEvent(QWheelEvent* event){
	last_user_event = boost::posix_time::second_clock::local_time();

	if(manipulatedClipPlane<0){ QGLViewer::wheelEvent(event); return; }
	/*assert(manipulatedClipPlane<renderer->numClipPlanes);
	float distStep=1e-3*sceneRadius();
	float dist=event->delta()*manipulatedFrame()->wheelSensitivity()*distStep;
	Vector3r normal=renderer->clipPlaneSe3[manipulatedClipPlane].orientation*Vector3r(0,0,1);
	qglviewer::Vec newPos=manipulatedFrame()->position()+qglviewer::Vec(normal[0],normal[1],normal[2])*dist;
	manipulatedFrame()->setPosition(newPos);
	renderer->clipPlaneSe3[manipulatedClipPlane].position=Vector3r(newPos[0],newPos[1],newPos[2]);
	*/
	updateGL();
	/* in draw, bound cutting planes will be moved as well */
}
