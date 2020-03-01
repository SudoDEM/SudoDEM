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


#include<QtGui/qevent.h>

#ifdef SUDODEM_GL2PS
	#include<gl2ps.h>
#endif

static int last(-1);

void GLViewer::useDisplayParameters(size_t n){
	LOG_DEBUG("Loading display parameters from #"<<n);
	vector<shared_ptr<DisplayParameters> >& dispParams=Omega::instance().getScene()->dispParams;
	if(dispParams.size()<=(size_t)n){ throw std::invalid_argument(("Display parameters #"+boost::lexical_cast<string>(n)+" don't exist (number of entries "+boost::lexical_cast<string>(dispParams.size())+")").c_str());; return;}
	const shared_ptr<DisplayParameters>& dp=dispParams[n];
	string val;
	if(dp->getValue("OpenGLRenderer",val)){ istringstream oglre(val);
		sudodem::ObjectIO::load<decltype(renderer),boost::archive::xml_iarchive>(oglre,"renderer",renderer);
	}
	else { LOG_WARN("OpenGLRenderer configuration not found in display parameters, skipped.");}
	if(dp->getValue("GLViewer",val)){ GLViewer::setState(val); displayMessage("Loaded view configuration #"+boost::lexical_cast<string>(n)); }
	else { LOG_WARN("GLViewer configuration not found in display parameters, skipped."); }
}

void GLViewer::saveDisplayParameters(size_t n){
	LOG_DEBUG("Saving display parameters to #"<<n);
	vector<shared_ptr<DisplayParameters> >& dispParams=Omega::instance().getScene()->dispParams;
	if(dispParams.size()<=n){while(dispParams.size()<=n) dispParams.push_back(shared_ptr<DisplayParameters>(new DisplayParameters));} assert(n<dispParams.size());
	shared_ptr<DisplayParameters>& dp=dispParams[n];
	ostringstream oglre;
	sudodem::ObjectIO::save<decltype(renderer),boost::archive::xml_oarchive>(oglre,"renderer",renderer);
	dp->setValue("OpenGLRenderer",oglre.str());
	dp->setValue("GLViewer",GLViewer::getState());
	displayMessage("Saved view configuration ot #"+boost::lexical_cast<string>(n));
}

void GLViewer::draw()
{
#ifdef SUDODEM_GL2PS
	if(!nextFrameSnapshotFilename.empty() && boost::algorithm::ends_with(nextFrameSnapshotFilename,".pdf")){
		gl2psStream=fopen(nextFrameSnapshotFilename.c_str(),"wb");
		if(!gl2psStream){ int err=errno; throw runtime_error(string("Error opening file ")+nextFrameSnapshotFilename+": "+strerror(err)); }
		LOG_DEBUG("Start saving snapshot to "<<nextFrameSnapshotFilename);
		size_t nBodies=Omega::instance().getScene()->bodies->size();
		int sortAlgo=(nBodies<100 ? GL2PS_BSP_SORT : GL2PS_SIMPLE_SORT);
		gl2psBeginPage(/*const char *title*/"Some title", /*const char *producer*/ "SudoDEM",
			/*GLint viewport[4]*/ NULL,
			/*GLint format*/ GL2PS_PDF, /*GLint sort*/ sortAlgo, /*GLint options*/GL2PS_SIMPLE_LINE_OFFSET|GL2PS_USE_CURRENT_VIEWPORT|GL2PS_TIGHT_BOUNDING_BOX|GL2PS_COMPRESS|GL2PS_OCCLUSION_CULL|GL2PS_NO_BLENDING,
			/*GLint colormode*/ GL_RGBA, /*GLint colorsize*/0,
			/*GL2PSrgba *colortable*/NULL,
			/*GLint nr*/0, /*GLint ng*/0, /*GLint nb*/0,
			/*GLint buffersize*/4096*4096 /* 16MB */, /*FILE *stream*/ gl2psStream,
			/*const char *filename*/NULL);
	}
#endif

	qglviewer::Vec vd=camera()->viewDirection(); renderer->viewDirection=Vector3r(vd[0],vd[1],vd[2]);
	if(Omega::instance().getScene()){
		const shared_ptr<Scene>& scene=Omega::instance().getScene();
		int selection = selectedName();
		if(selection!=-1 && (*(Omega::instance().getScene()->bodies)).exists(selection) && isMoving){
			static Real lastTimeMoved(0);
			qreal v0,v1,v2; manipulatedFrame()->getPosition(v0,v1,v2);//zhswee
			if(last == selection) // delay by one redraw, so the body will not jump into 0,0,0 coords
			{
				Quaternionr& q = (*(Omega::instance().getScene()->bodies))[selection]->state->ori;
				Vector3r&    v = (*(Omega::instance().getScene()->bodies))[selection]->state->pos;
				Vector3r&    vel = (*(Omega::instance().getScene()->bodies))[selection]->state->vel;
				Vector3r&    angVel = (*(Omega::instance().getScene()->bodies))[selection]->state->angVel;
				angVel=Vector3r::Zero();
				Real dt=(scene->time-lastTimeMoved); lastTimeMoved=scene->time;
				if (dt!=0) { vel[0]=-(v[0]-v0)/dt; vel[1]=-(v[1]-v1)/dt; vel[2]=-(v[2]-v2)/dt;}
				else vel[0]=vel[1]=vel[2]=0;
				//FIXME: should update spin like velocity above, when the body is rotated:
				double q0,q1,q2,q3; manipulatedFrame()->getOrientation(q0,q1,q2,q3);	q.x()=q0;q.y()=q1;q.z()=q2;q.w()=q3;
			}
			(*(Omega::instance().getScene()->bodies))[selection]->userForcedDisplacementRedrawHook();
		}
		if(manipulatedClipPlane>=0){
			assert(manipulatedClipPlane<renderer->numClipPlanes);
			qreal v0,v1,v2; manipulatedFrame()->getPosition(v0,v1,v2);//zhswee
			double q0,q1,q2,q3; manipulatedFrame()->getOrientation(q0,q1,q2,q3);
			Se3r newSe3(Vector3r(v0,v1,v2),Quaternionr(q0,q1,q2,q3)); newSe3.orientation.normalize();
			const Se3r& oldSe3=renderer->clipPlaneSe3[manipulatedClipPlane];
			FOREACH(int planeId, boundClipPlanes){
				if(planeId>=renderer->numClipPlanes || !renderer->clipPlaneActive[planeId] || planeId==manipulatedClipPlane) continue;
				Se3r& boundSe3=renderer->clipPlaneSe3[planeId];
				Quaternionr relOrient=oldSe3.orientation.conjugate()*boundSe3.orientation; relOrient.normalize();
				Vector3r relPos=oldSe3.orientation.conjugate()*(boundSe3.position-oldSe3.position);
				boundSe3.position=newSe3.position+newSe3.orientation*relPos;
				boundSe3.orientation=newSe3.orientation*relOrient;
				boundSe3.orientation.normalize();
			}
			renderer->clipPlaneSe3[manipulatedClipPlane]=newSe3;
		}
		scene->renderer=renderer;
		renderer->render(scene, selectedName());
	}
}

void GLViewer::drawWithNames(){
	qglviewer::Vec vd=camera()->viewDirection(); renderer->viewDirection=Vector3r(vd[0],vd[1],vd[2]);
	const shared_ptr<Scene> scene(Omega::instance().getScene());
	scene->renderer=renderer;
	renderer->scene=scene;
	renderer->renderShape();
}


qglviewer::Vec GLViewer::displayedSceneCenter(){
	return camera()->unprojectedCoordinatesOf(qglviewer::Vec(width()/2 /* pixels */ ,height()/2 /* pixels */, /*middle between near plane and far plane*/ .5));
}

float GLViewer::displayedSceneRadius(){
	return (camera()->unprojectedCoordinatesOf(qglviewer::Vec(width()/2,height()/2,.5))-camera()->unprojectedCoordinatesOf(qglviewer::Vec(0,0,.5))).norm();
}

void GLViewer::postDraw(){
	Real wholeDiameter=QGLViewer::camera()->sceneRadius()*2;

	renderer->viewInfo.sceneRadius=QGLViewer::camera()->sceneRadius();
	qglviewer::Vec c=QGLViewer::camera()->sceneCenter();
	renderer->viewInfo.sceneCenter=Vector3r(c[0],c[1],c[2]);

	Real dispDiameter=min(wholeDiameter,max((Real)displayedSceneRadius()*2,wholeDiameter/1e3)); // limit to avoid drawing 1e5 lines with big zoom level
	//qglviewer::Vec center=QGLViewer::camera()->sceneCenter();
	Real gridStep=pow(10,(floor(log10(dispDiameter)-.7)));
	Real scaleStep=pow(10,(floor(log10(displayedSceneRadius()*2)-.7))); // unconstrained
	int nSegments=((int)(wholeDiameter/gridStep))+1;
	Real realSize=nSegments*gridStep;
	//LOG_TRACE("nSegments="<<nSegments<<",gridStep="<<gridStep<<",realSize="<<realSize);
	glPushMatrix();

	nSegments *= 2; // there's an error in QGLViewer::drawGrid(), so we need to mitigate it by '* 2'
	// XYZ grids
	glLineWidth(.5);
	if(drawGrid & 1) {glColor3f(0.6,0.3,0.3); glPushMatrix(); glRotated(90.,0.,1.,0.); QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(drawGrid & 2) {glColor3f(0.3,0.6,0.3); glPushMatrix(); glRotated(90.,1.,0.,0.); QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(drawGrid & 4) {glColor3f(0.3,0.3,0.6); glPushMatrix(); /*glRotated(90.,0.,1.,0.);*/ QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(gridSubdivide){
		if(drawGrid & 1) {glColor3f(0.4,0.1,0.1); glPushMatrix(); glRotated(90.,0.,1.,0.); QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
		if(drawGrid & 2) {glColor3f(0.1,0.4,0.1); glPushMatrix(); glRotated(90.,1.,0.,0.); QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
		if(drawGrid & 4) {glColor3f(0.1,0.1,0.4); glPushMatrix(); /*glRotated(90.,0.,1.,0.);*/ QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
	}

	// scale
	if(drawScale){
		Real segmentSize=scaleStep;
		qglviewer::Vec screenDxDy[3]; // dx,dy for x,y,z scale segments
		int extremalDxDy[2]={0,0};
		for(int axis=0; axis<3; axis++){
			qglviewer::Vec delta(0,0,0); delta[axis]=segmentSize;
			qglviewer::Vec center=displayedSceneCenter();
			screenDxDy[axis]=camera()->projectedCoordinatesOf(center+delta)-camera()->projectedCoordinatesOf(center);
			for(int xy=0;xy<2;xy++)extremalDxDy[xy]=(axis>0 ? min(extremalDxDy[xy],(int)screenDxDy[axis][xy]) : screenDxDy[axis][xy]);
		}
		//LOG_DEBUG("Screen offsets for axes: "<<" x("<<screenDxDy[0][0]<<","<<screenDxDy[0][1]<<") y("<<screenDxDy[1][0]<<","<<screenDxDy[1][1]<<") z("<<screenDxDy[2][0]<<","<<screenDxDy[2][1]<<")");
		int margin=10; // screen pixels
		int scaleCenter[2]; scaleCenter[0]=std::abs(extremalDxDy[0])+margin; scaleCenter[1]=std::abs(extremalDxDy[1])+margin;
		//LOG_DEBUG("Center of scale "<<scaleCenter[0]<<","<<scaleCenter[1]);
		//displayMessage(QString().sprintf("displayed scene radius %g",displayedSceneRadius()));
		startScreenCoordinatesSystem();
			glDisable(GL_LIGHTING);
			glDisable(GL_DEPTH_TEST);
			glLineWidth(3.0);
			for(int axis=0; axis<3; axis++){
				Vector3r color(.4,.4,.4); color[axis]=.9;
				glColor3v(color);
				glBegin(GL_LINES);
				glVertex2f(scaleCenter[0],scaleCenter[1]);
				glVertex2f(scaleCenter[0]+screenDxDy[axis][0],scaleCenter[1]+screenDxDy[axis][1]);
				glEnd();
			}
			glLineWidth(1.);
			glEnable(GL_DEPTH_TEST);
			QGLViewer::drawText(scaleCenter[0],scaleCenter[1],QString().sprintf("%.3g",(double)scaleStep));
		stopScreenCoordinatesSystem();
	}

	// cutting planes (should be moved to OpenGLRenderer perhaps?)
	// only painted if one of those is being manipulated
	if(manipulatedClipPlane>=0){
		for(int planeId=0; planeId<renderer->numClipPlanes; planeId++){
			if(!renderer->clipPlaneActive[planeId] && planeId!=manipulatedClipPlane) continue;
			glPushMatrix();
				const Se3r& se3=renderer->clipPlaneSe3[planeId];
				AngleAxisr aa(se3.orientation);
				glTranslatef(se3.position[0],se3.position[1],se3.position[2]);
				glRotated(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
				Real cff=1;
				if(!renderer->clipPlaneActive[planeId]) cff=.4;
				glColor3f(max((Real)0.,cff*cos(planeId)),max((Real)0.,cff*sin(planeId)),planeId==manipulatedClipPlane); // variable colors
				QGLViewer::drawGrid(realSize,2*nSegments);
				drawArrow(wholeDiameter/6);
			glPopMatrix();
		}
	}

	Scene* rb=Omega::instance().getScene().get();
	#define _W3 setw(3)<<setfill('0')
	#define _W2 setw(2)<<setfill('0')
	if(timeDispMask!=0){
		const int lineHt=13;
		unsigned x=10,y=height()-3-lineHt*2;
		glColor3v(Vector3r(1,1,1));
		if(timeDispMask & GLViewer::TIME_VIRT){
			ostringstream oss;
			const Real& t=Omega::instance().getScene()->time;
			unsigned min=((unsigned)t/60),sec=(((unsigned)t)%60),msec=((unsigned)(1e3*t))%1000,usec=((unsigned long)(1e6*t))%1000,nsec=((unsigned long)(1e9*t))%1000;
			if(min>0) oss<<_W2<<min<<":"<<_W2<<sec<<"."<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (sec>0) oss<<_W2<<sec<<"."<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (msec>0) oss<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (usec>0) oss<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else oss<<_W3<<nsec<<"ns";
			QGLViewer::drawText(x,y,oss.str().c_str());
			y-=lineHt;
		}
		glColor3v(Vector3r(0,.5,.5));
		if(timeDispMask & GLViewer::TIME_REAL){
			QGLViewer::drawText(x,y,getRealTimeString().c_str() /* virtual, since player gets that from db */);
			y-=lineHt;
		}
		if(timeDispMask & GLViewer::TIME_ITER){
			ostringstream oss;
			oss<<"#"<<rb->iter;
			if(rb->stopAtIter>rb->iter) oss<<" ("<<setiosflags(ios::fixed)<<setw(3)<<setprecision(1)<<setfill('0')<<(100.*rb->iter)/rb->stopAtIter<<"%)";
			QGLViewer::drawText(x,y,oss.str().c_str());
			y-=lineHt;
		}
		if(drawGrid){
			glColor3v(Vector3r(1,1,0));
			ostringstream oss;
			oss<<"grid: "<<setprecision(4)<<gridStep;
			if(gridSubdivide) oss<<" (minor "<<setprecision(3)<<gridStep*.1<<")";
			QGLViewer::drawText(x,y,oss.str().c_str());
			y-=lineHt;
		}
	}
	QGLViewer::postDraw();
	if(!nextFrameSnapshotFilename.empty()){
		#ifdef SUDODEM_GL2PS
			if(boost::algorithm::ends_with(nextFrameSnapshotFilename,".pdf")){
				gl2psEndPage();
				LOG_DEBUG("Finished saving snapshot to "<<nextFrameSnapshotFilename);
				fclose(gl2psStream);
			} else
	#endif
		{
			// save the snapshot
			saveSnapshot(QString(nextFrameSnapshotFilename.c_str()),/*overwrite*/ true);
		}
		// notify the caller that it is done already (probably not an atomic op :-|, though)
		nextFrameSnapshotFilename.clear();
	}
}


