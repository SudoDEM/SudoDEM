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

static unsigned initBlocked(State::DOF_NONE);

CREATE_LOGGER(GLViewer);

GLLock::GLLock(GLViewer* _glv): boost::try_mutex::scoped_lock(Omega::instance().renderMutex), glv(_glv){
	glv->makeCurrent();
}
GLLock::~GLLock(){ glv->doneCurrent(); }


#define _W3 setw(3)<<setfill('0')
#define _W2 setw(2)<<setfill('0')
GLViewer::~GLViewer(){ /* get the GL mutex when closing */ GLLock lock(this); /* cerr<<"Destructing view #"<<viewId<<endl;*/ }

void GLViewer::closeEvent(QCloseEvent *e){
	LOG_DEBUG("Will emit closeView for view #"<<viewId);
	OpenGLManager::self->emitCloseView(viewId);
	e->accept();
}

GLViewer::GLViewer(int _viewId, const shared_ptr<OpenGLRenderer>& _renderer, QGLWidget* shareWidget): QGLViewer(/*parent*/(QWidget*)NULL,shareWidget), renderer(_renderer), viewId(_viewId) {
	isMoving=false;
	drawGrid=0;
	drawScale=false;
	timeDispMask=TIME_REAL|TIME_VIRT|TIME_ITER;
	cut_plane = 0;
	cut_plane_delta = -2;
	gridSubdivide = false;
	resize(550,550);
	last=-1;
	if(viewId==0) setWindowTitle("Primary view");
	else setWindowTitle(("Secondary view #"+boost::lexical_cast<string>(viewId)).c_str());

	show();

	mouseMovesCamera();

	 camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
	manipulatedClipPlane=-1;

	if(manipulatedFrame()==0) setManipulatedFrame(new qglviewer::ManipulatedFrame());

	xyPlaneConstraint=shared_ptr<qglviewer::LocalConstraint>(new qglviewer::LocalConstraint());
	// Forbid rotation
  qglviewer::WorldConstraint *constraint = new qglviewer::WorldConstraint();
  constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
	//manipulatedFrame()->setConstraint(NULL);
	manipulatedFrame()->setConstraint(constraint);

	setKeyDescription(Qt::Key_Return,"Run simulation.");
	setKeyDescription(Qt::Key_A,"Toggle visibility of global axes.");
	setKeyDescription(Qt::Key_C,"Set scene center so that all bodies are visible; if a body is selected, center around this body.");
	setKeyDescription(Qt::Key_C & Qt::AltModifier,"Set scene center to median body position (same as space)");
	setKeyDescription(Qt::Key_D,"Toggle time display mask");
	setKeyDescription(Qt::Key_G,"Toggle grid visibility; g turns on and cycles");
	setKeyDescription(Qt::Key_G & Qt::ShiftModifier ,"Hide grid.");
	setKeyDescription(Qt::Key_M, "Move selected object.");
	//setKeyDescription(Qt::Key_X,"Show the xz [shift: xy] (up-right) plane (clip plane: align normal with +x)");
	//setKeyDescription(Qt::Key_Y,"Show the yx [shift: yz] (up-right) plane (clip plane: align normal with +y)");
	//setKeyDescription(Qt::Key_Z,"Show the zy [shift: zx] (up-right) plane (clip plane: align normal with +z)");
	setKeyDescription(Qt::Key_Period,"Toggle grid subdivision by 10");
	setKeyDescription(Qt::Key_S & Qt::AltModifier,"Save QGLViewer state to /tmp/qglviewerState.xml");
	setKeyDescription(Qt::Key_T,"Switch orthographic / perspective camera");
	setKeyDescription(Qt::Key_O,"Set narrower field of view");
	setKeyDescription(Qt::Key_P,"Set wider field of view");
	setKeyDescription(Qt::Key_R,"Revolve around scene center");
	setKeyDescription(Qt::Key_V,"Save PDF of the current view to /tmp/sudodem-snapshot-0001.pdf (whichever number is available first). (Must be compiled with the gl2ps feature.)");
 	//setPathKey(-Qt::Key_F1);
 	//setPathKey(-Qt::Key_F2);
	setKeyDescription(Qt::Key_Escape,"Manipulate scene (default)");
	//setKeyDescription(Qt::Key_F1,"Manipulate clipping plane #1");
	//setKeyDescription(Qt::Key_F2,"Manipulate clipping plane #2");
	//setKeyDescription(Qt::Key_F3,"Manipulate clipping plane #3");
	//setKeyDescription(Qt::Key_1,"Make the manipulated clipping plane parallel with plane #1");
	//setKeyDescription(Qt::Key_2,"Make the manipulated clipping plane parallel with plane #2");
	//setKeyDescription(Qt::Key_2,"Make the manipulated clipping plane parallel with plane #3");
	//setKeyDescription(Qt::Key_1 & Qt::AltModifier,"Add/remove plane #1 to/from the bound group");
	//setKeyDescription(Qt::Key_2 & Qt::AltModifier,"Add/remove plane #2 to/from the bound group");
	//setKeyDescription(Qt::Key_3 & Qt::AltModifier,"Add/remove plane #3 to/from the bound group");
	//setKeyDescription(Qt::Key_0,"Clear the bound group");
	//setKeyDescription(Qt::Key_7,"Load [Alt: save] view configuration #0");
	//setKeyDescription(Qt::Key_8,"Load [Alt: save] view configuration #1");
	//setKeyDescription(Qt::Key_9,"Load [Alt: save] view configuration #2");
	//setKeyDescription(Qt::Key_Space,"Center scene (same as Alt-C); clip plane: activate/deactivate");

	centerScene();
}

bool GLViewer::isManipulating(){
	return isMoving || manipulatedClipPlane>=0;
}

void GLViewer::resetManipulation(){
	mouseMovesCamera();
	setSelectedName(-1);
	isMoving=false;
	manipulatedClipPlane=-1;
}
/*
void GLViewer::startClipPlaneManipulation(int planeNo){
	assert(planeNo<renderer->numClipPlanes);
	resetManipulation();
	mouseMovesManipulatedFrame(xyPlaneConstraint.get());
	manipulatedClipPlane=planeNo;
	const Se3r se3(renderer->clipPlaneSe3[planeNo]);
	manipulatedFrame()->setPositionAndOrientation(qglviewer::Vec(se3.position[0],se3.position[1],se3.position[2]),qglviewer::Quaternion(se3.orientation.x(),se3.orientation.y(),se3.orientation.z(),se3.orientation.w()));
	string grp=strBoundGroup();
	displayMessage("Manipulating clip plane #"+boost::lexical_cast<string>(planeNo+1)+(grp.empty()?grp:" (bound planes:"+grp+")"));
}*/

string GLViewer::getState(){
	QString origStateFileName=stateFileName();
	string tmpFile=Omega::instance().tmpFilename();
	setStateFileName(QString(tmpFile.c_str())); saveStateToFile(); setStateFileName(origStateFileName);
	LOG_DEBUG("State saved to temp file "<<tmpFile);
	// read tmp file contents and return it as string
	// this will replace all whitespace by space (nowlines will disappear, which is what we want)
	ifstream in(tmpFile.c_str()); string ret; while(!in.eof()){string ss; in>>ss; ret+=" "+ss;}; in.close();
	boost::filesystem::remove(boost::filesystem::path(tmpFile));
	return ret;
}

void GLViewer::setState(string state){
	string tmpFile=Omega::instance().tmpFilename();
	std::ofstream out(tmpFile.c_str());
	if(!out.good()){ LOG_ERROR("Error opening temp file `"<<tmpFile<<"', loading aborted."); return; }
	out<<state; out.close();
	LOG_DEBUG("Will load state from temp file "<<tmpFile);
	QString origStateFileName=stateFileName(); setStateFileName(QString(tmpFile.c_str())); restoreStateFromFile(); setStateFileName(origStateFileName);
	boost::filesystem::remove(boost::filesystem::path(tmpFile));
}

void GLViewer::keyPressEvent(QKeyEvent *e)
{
	last_user_event = boost::posix_time::second_clock::local_time();
/*
	if(false){}
	// special keys: Escape and Space
	else if(e->key()==Qt::Key_A){ toggleAxisIsDrawn(); return; }
	else if(e->key()==Qt::Key_Escape){
		if(!isManipulating()){ return; }
		else { resetManipulation(); displayMessage("Manipulating scene."); }
	}
	else if(e->key()==Qt::Key_Space){
		if(manipulatedClipPlane>=0) {displayMessage("Clip plane #"+boost::lexical_cast<string>(manipulatedClipPlane+1)+(renderer->clipPlaneActive[manipulatedClipPlane]?" de":" ")+"activated"); renderer->clipPlaneActive[manipulatedClipPlane]=!renderer->clipPlaneActive[manipulatedClipPlane]; }
		else{ centerMedianQuartile(); }
	}
	// function keys
	else if(e->key()==Qt::Key_F1 || e->key()==Qt::Key_F2 || e->key()==Qt::Key_F3  ){
		int n=0; if(e->key()==Qt::Key_F1) n=1; else if(e->key()==Qt::Key_F2) n=2; else if(e->key()==Qt::Key_F3) n=3; assert(n>0); int planeId=n-1;
		if(planeId>=renderer->numClipPlanes) return;
		//if(planeId!=manipulatedClipPlane) startClipPlaneManipulation(planeId);
	}
	// numbers
	else if(e->key()==Qt::Key_0 && (e->modifiers() & Qt::AltModifier)) { boundClipPlanes.clear(); displayMessage("Cleared bound planes group.");}
	else if(e->key()==Qt::Key_1 || e->key()==Qt::Key_2 || e->key()==Qt::Key_3 ){
		int n=0; if(e->key()==Qt::Key_1) n=1; else if(e->key()==Qt::Key_2) n=2; else if(e->key()==Qt::Key_3) n=3; assert(n>0); int planeId=n-1;
		if(planeId>=renderer->numClipPlanes) return; // no such clipping plane
		if(e->modifiers() & Qt::AltModifier){
			if(boundClipPlanes.count(planeId)==0) {boundClipPlanes.insert(planeId); displayMessage("Added plane #"+boost::lexical_cast<string>(planeId+1)+" to the bound group: "+strBoundGroup());}
			else {boundClipPlanes.erase(planeId); displayMessage("Removed plane #"+boost::lexical_cast<string>(planeId+1)+" from the bound group: "+strBoundGroup());}
		}
		else if(manipulatedClipPlane>=0 && manipulatedClipPlane!=planeId) {
			const Quaternionr& o=renderer->clipPlaneSe3[planeId].orientation;
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(o.x(),o.y(),o.z(),o.w()));
			displayMessage("Copied orientation from plane #1");
		}
	}
	else if(e->key()==Qt::Key_7 || e->key()==Qt::Key_8 || e->key()==Qt::Key_9){
		int nn=-1; if(e->key()==Qt::Key_7)nn=0; else if(e->key()==Qt::Key_8)nn=1; else if(e->key()==Qt::Key_9)nn=2; assert(nn>=0); size_t n=(size_t)nn;
		if(e->modifiers() & Qt::AltModifier) saveDisplayParameters(n);
		else useDisplayParameters(n);
	}
	// letters alphabetically
	else if(e->key()==Qt::Key_C && (e->modifiers() & Qt::AltModifier)){ displayMessage("Median centering"); centerMedianQuartile(); }
	else if(e->key()==Qt::Key_C){
		// center around selected body
		if(selectedName() >= 0 && (*(Omega::instance().getScene()->bodies)).exists(selectedName())) setSceneCenter(manipulatedFrame()->position());
		// make all bodies visible
		else centerScene();
	}
	else if(e->key()==Qt::Key_D &&(e->modifiers() & Qt::AltModifier)){ Body::id_t id; if((id=Omega::instance().getScene()->selectedBody)>=0){ const shared_ptr<Body>& b=Body::byId(id); b->setDynamic(!b->isDynamic()); LOG_INFO("Body #"<<id<<" now "<<(b->isDynamic()?"":"NOT")<<" dynamic"); } }
	else if(e->key()==Qt::Key_D) {timeDispMask+=1; if(timeDispMask>(TIME_REAL|TIME_VIRT|TIME_ITER))timeDispMask=0; }
	else if(e->key()==Qt::Key_G) { if(e->modifiers() & Qt::ShiftModifier){ drawGrid=0; return; } else drawGrid++; if(drawGrid>=8) drawGrid=0; }
	else if (e->key()==Qt::Key_M && selectedName() >= 0){
		if(!(isMoving=!isMoving)){
			displayMessage("Moving done.");
			if (last>=0) {Body::byId(Body::id_t(last))->state->blockedDOFs=initBlocked; last=-1;} mouseMovesCamera();}
		else{ 	displayMessage("Moving selected object");

			long selection = Omega::instance().getScene()->selectedBody;
			initBlocked=Body::byId(Body::id_t(selection))->state->blockedDOFs; last=selection;
			Body::byId(Body::id_t(selection))->state->blockedDOFs=State::DOF_ALL;

			Rotationr& r = Body::byId(selection)->state->ori;
			Quaternionr q = Quaternionr(AngleAxisr(r.angle(),Vector3r::UnitZ()));
			Vector2r&    v = Body::byId(selection)->state->pos;
			manipulatedFrame()->setPositionAndOrientation(qglviewer::Vec(v[0],v[1],0.0),qglviewer::Quaternion(q.x(),q.y(),q.z(),q.w()));
			mouseMovesManipulatedFrame();}
	}
	else if (e->key() == Qt::Key_T) camera()->setType(camera()->type()==qglviewer::Camera::ORTHOGRAPHIC ? qglviewer::Camera::PERSPECTIVE : qglviewer::Camera::ORTHOGRAPHIC);
	else if(e->key()==Qt::Key_O) camera()->setFieldOfView(camera()->fieldOfView()*0.9);
	else if(e->key()==Qt::Key_P) camera()->setFieldOfView(camera()->fieldOfView()*1.1);
	else if(e->key()==Qt::Key_R){ // reverse the clipping plane; revolve around scene center if no clipping plane selected
		if(manipulatedClipPlane>=0 && manipulatedClipPlane<renderer->numClipPlanes){
			// here, we must update both manipulatedFrame orientation and renderer->clipPlaneSe3 orientation in the same way
			Quaternionr& ori=renderer->clipPlaneSe3[manipulatedClipPlane].orientation;
			ori=Quaternionr(AngleAxisr(Mathr::PI,Vector3r(0,1,0)))*ori;
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(qglviewer::Vec(0,1,0),Mathr::PI)*manipulatedFrame()->orientation());
			displayMessage("Plane #"+boost::lexical_cast<string>(manipulatedClipPlane+1)+" reversed.");
		}
		else {
			camera()->setRevolveAroundPoint(sceneCenter());
		}
	}
	else if(e->key()==Qt::Key_S){
		LOG_INFO("Saving QGLViewer state to /tmp/qglviewerState.xml");
		setStateFileName("/tmp/qglviewerState.xml"); saveStateToFile(); setStateFileName(QString::null);
	}
	else if(e->key()==Qt::Key_X || e->key()==Qt::Key_Y || e->key()==Qt::Key_Z){
		int axisIdx=(e->key()==Qt::Key_X?0:(e->key()==Qt::Key_Y?1:2));
		if(manipulatedClipPlane<0){
			qglviewer::Vec up(0,0,0), vDir(0,0,0);
			bool alt=(e->modifiers() && Qt::ShiftModifier);
			up[axisIdx]=1; vDir[(axisIdx+(alt?2:1))%3]=alt?1:-1;
			camera()->setViewDirection(vDir);
			camera()->setUpVector(up);
			centerMedianQuartile();
		}
		else{ // align clipping normal plane with world axis
			// x: (0,1,0),pi/2; y: (0,0,1),pi/2; z: (1,0,0),0
			qglviewer::Vec axis(0,0,0); axis[(axisIdx+1)%3]=1; Real angle=axisIdx==2?0:Mathr::PI/2;
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(axis,angle));
		}
	}
	else if(e->key()==Qt::Key_Period) gridSubdivide = !gridSubdivide;
	else if(e->key()==Qt::Key_Return){
		if (Omega::instance().isRunning()) Omega::instance().pause();
		else Omega::instance().run();
		LOG_INFO("Running...");
	}
#ifdef SUDODEM_GL2PS
	else if(e->key()==Qt::Key_V){
		for(int i=0; ;i++){
			std::ostringstream fss; fss<<"/tmp/sudodem-snapshot-"<<setw(4)<<setfill('0')<<i<<".pdf";
			if(!boost::filesystem::exists(fss.str())){ nextFrameSnapshotFilename=fss.str(); break; }
		}
		LOG_INFO("Will save snapshot to "<<nextFrameSnapshotFilename);
	}
#endif
#if 0
	else if( e->key()==Qt::Key_Plus ){
			cut_plane = std::min(1.0, cut_plane + std::pow(10.0,(double)cut_plane_delta));
			static_cast<SudoDEMCamera*>(camera())->setCuttingDistance(cut_plane);
			displayMessage("Cut plane: "+boost::lexical_cast<std::string>(cut_plane));
	}else if( e->key()==Qt::Key_Minus ){
			cut_plane = std::max(0.0, cut_plane - std::pow(10.0,(double)cut_plane_delta));
			static_cast<SudoDEMCamera*>(camera())->setCuttingDistance(cut_plane);
			displayMessage("Cut plane: "+boost::lexical_cast<std::string>(cut_plane));
	}else if( e->key()==Qt::Key_Slash ){
			cut_plane_delta -= 1;
			displayMessage("Cut plane increment: 1e"+(cut_plane_delta>0?std::string("+"):std::string(""))+boost::lexical_cast<std::string>(cut_plane_delta));
	}else if( e->key()==Qt::Key_Asterisk ){
			cut_plane_delta = std::min(1+cut_plane_delta,-1);
			displayMessage("Cut plane increment: 1e"+(cut_plane_delta>0?std::string("+"):std::string(""))+boost::lexical_cast<std::string>(cut_plane_delta));
	}
#endif

	else if(e->key()!=Qt::Key_Escape && e->key()!=Qt::Key_Space) QGLViewer::keyPressEvent(e);
*/
	updateGL();
}

/* Center the scene such that periodic cell is contained in the view */
void GLViewer::centerPeriodic(){
	Scene* scene=Omega::instance().getScene().get();
	assert(scene->isPeriodic);
	Vector2r center=.5*scene->cell->getSize();
	Vector2r halfSize=.5*scene->cell->getSize();
	float radius=std::max(halfSize[0],halfSize[1]);
	LOG_DEBUG("Periodic scene center="<<center<<", halfSize="<<halfSize<<", radius="<<radius);
	setSceneCenter(qglviewer::Vec(center[0],center[1],0));
	setSceneRadius(radius*1.5);
	showEntireScene();
}

/* Calculate medians for x, y and z coordinates of all bodies;
 *then set scene center to median position and scene radius to 2*inter-quartile distance.
 *
 * This function eliminates the effect of lonely bodies that went nuts and enlarge
 * the scene's Aabb in such a way that fitting the scene to see the Aabb makes the
 * "central" (where most bodies is) part very small or even invisible.
 */
void GLViewer::centerMedianQuartile(){
	Scene* scene=Omega::instance().getScene().get();
	if(scene->isPeriodic){ centerPeriodic(); return; }
	long nBodies=scene->bodies->size();
  //  nBodies +=scene->nodes->size();//including nodes from fem
	if(nBodies<4) {
		LOG_DEBUG("Less than 4 bodies, median makes no sense; calling centerScene() instead.");
		return centerScene();
	}
	std::vector<Real> coords[2];
	for(int i=0;i<2;i++)coords[i].reserve(nBodies);
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if(!b) continue;
		for(int i=0; i<2; i++) coords[i].push_back(b->state->pos[i]);
	}
	/*FOREACH(shared_ptr<Node> b, *scene->nodes){
		if(!b) continue;
		for(int i=0; i<3; i++) coords[i].push_back(b->pos[i]);
	}*/
	Vector2r median,interQuart;
	for(int i=0;i<2;i++){
		sort(coords[i].begin(),coords[i].end());
		median[i]=*(coords[i].begin()+nBodies/2);
		interQuart[i]=*(coords[i].begin()+2*nBodies/4)-*(coords[i].begin()+nBodies/4);//debug for 2d
	}
	LOG_DEBUG("Median position is"<<median<<", inter-quartile distance is "<<interQuart);
	setSceneCenter(qglviewer::Vec(median[0],median[1],0));
	setSceneRadius(2*(interQuart[0]+interQuart[1])/3.);
	showEntireScene();
}

void GLViewer::centerScene(){
	Scene* rb=Omega::instance().getScene().get();
	if (!rb) return;
	if(rb->isPeriodic){ centerPeriodic(); return; }
	LOG_INFO("Select with shift, press 'm' to move.");
	Vector2r min,max;
	if(not(rb->bound)){ rb->updateBound();}

	min=rb->bound->min; max=rb->bound->max;
	bool hasNan=(isnan(min[0])||isnan(min[1])||isnan(max[0])||isnan(max[1]));
	Real minDim=std::min(max[0]-min[0],max[1]-min[1]);
	if(minDim<=0 || hasNan){
		// Aabb is not yet calculated...
		LOG_DEBUG("scene's bound not yet calculated or has zero or nan dimension(s), attempt get that from bodies' positions.");
		Real inf=std::numeric_limits<Real>::infinity();
		min=Vector2r(inf,inf); max=Vector2r(-inf,-inf);
		FOREACH(const shared_ptr<Body>& b, *rb->bodies){
			if(!b) continue;
			max=max.cwiseMax(b->state->pos);
			min=min.cwiseMin(b->state->pos);
		}
      /*  FOREACH(const shared_ptr<Node>& n, *rb->nodes){
			if(!n) continue;
			max=max.cwiseMax(n->pos);
			min=min.cwiseMin(n->pos);
		}*/
		if(isinf(min[0])||isinf(min[1])||isinf(max[0])||isinf(max[1])){ LOG_DEBUG("No min/max computed from bodies either, setting cube (-1,-1)×(1,1)"); min=-Vector2r::Ones(); max=Vector2r::Ones(); }
	} else {LOG_DEBUG("Using scene's Aabb");}

	LOG_DEBUG("Got scene box min="<<min<<" and max="<<max);
	Vector2r center = (max+min)*0.5;
	Vector2r halfSize = (max-min)*0.5;
	float radius=std::max(halfSize[0],halfSize[1]); if(radius<=0) radius=1;
	LOG_DEBUG("Scene center="<<center<<", halfSize="<<halfSize<<", radius="<<radius);
	setSceneCenter(qglviewer::Vec(center[0],center[1],0.0));
	setSceneRadius(radius*1.5);
	showEntireScene();
}

// new object selected.
// set frame coordinates, and isDynamic=false;
void GLViewer::postSelection(const QPoint& point)
{
	LOG_DEBUG("Selection is "<<selectedName());
	int selection = selectedName();
	if(selection<0){
		if (last>=0) {
			Body::byId(Body::id_t(last))->state->blockedDOFs=initBlocked; last=-1; Omega::instance().getScene()->selectedBody = -1;}
		if(isMoving){
			displayMessage("Moving finished"); mouseMovesCamera(); isMoving=false;
			Omega::instance().getScene()->selectedBody = -1;
		}
		return;
	}
	if(selection>=0 && (*(Omega::instance().getScene()->bodies)).exists(selection)){
		resetManipulation();
		if (last>=0) {Body::byId(Body::id_t(last))->state->blockedDOFs=initBlocked; last=-1;}
		if(Body::byId(Body::id_t(selection))->isClumpMember()){ // select clump (invisible) instead of its member
			LOG_DEBUG("Clump member #"<<selection<<" selected, selecting clump instead.");
			selection=Body::byId(Body::id_t(selection))->clumpId;
		}

		setSelectedName(selection);
		LOG_DEBUG("New selection "<<selection);
		displayMessage("Selected body #"+boost::lexical_cast<string>(selection)+(Body::byId(selection)->isClump()?" (clump)":""));
		Omega::instance().getScene()->selectedBody = selection;
			PyGILState_STATE gstate;
			gstate = PyGILState_Ensure();
      boost::python::object main=boost::python::import("__main__");
      boost::python::object global=main.attr("__dict__");
				// the try/catch block must be properly nested inside PyGILState_Ensure and PyGILState_Release
				try{
          boost::python::eval(string("onBodySelect("+boost::lexical_cast<string>(selection)+")").c_str(),global,global);
				} catch (boost::python::error_already_set const &) {
					LOG_DEBUG("unable to call onBodySelect. Not defined?");
				}
			PyGILState_Release(gstate);
			// see https://svn.boost.org/trac/boost/ticket/2781 for exception handling
	}
}

// maybe new object will be selected.
// if so, then set isDynamic of previous selection, to old value
void GLViewer::endSelection(const QPoint &point){
	manipulatedClipPlane=-1;
	QGLViewer::endSelection(point);
}

string GLViewer::getRealTimeString(){
	ostringstream oss;
  boost::posix_time::time_duration t=Omega::instance().getRealTime_duration();
	unsigned d=t.hours()/24,h=t.hours()%24,m=t.minutes(),s=t.seconds();
	oss<<"clock ";
	if(d>0) oss<<d<<"days "<<_W2<<h<<":"<<_W2<<m<<":"<<_W2<<s;
	else if(h>0) oss<<_W2<<h<<":"<<_W2<<m<<":"<<_W2<<s;
	else oss<<_W2<<m<<":"<<_W2<<s;
	return oss.str();
}
#undef _W2
#undef _W3

// cut&paste from QGLViewer::domElement documentation
QDomElement GLViewer::domElement(const QString& name, QDomDocument& document) const{
	QDomElement de=document.createElement("grid");
	string val; if(drawGrid & 1) val+="x"; if(drawGrid & 2)val+="y"; if(drawGrid & 4)val+="z";
	de.setAttribute("normals",val.c_str());
	QDomElement de2=document.createElement("timeDisplay"); de2.setAttribute("mask",timeDispMask);
	QDomElement res=QGLViewer::domElement(name,document);
	res.appendChild(de);
	res.appendChild(de2);
	return res;
}

// cut&paste from QGLViewer::initFromDomElement documentation
void GLViewer::initFromDOMElement(const QDomElement& element){
	QGLViewer::initFromDOMElement(element);
	QDomElement child=element.firstChild().toElement();
	while (!child.isNull()){
		if (child.tagName()=="gridXYZ" && child.hasAttribute("normals")){
			string val=child.attribute("normals").toLower().toStdString();
			drawGrid=0;
			if(val.find("x")!=string::npos) drawGrid+=1; if(val.find("y")!=string::npos)drawGrid+=2; if(val.find("z")!=string::npos)drawGrid+=4;
		}
		if(child.tagName()=="timeDisplay" && child.hasAttribute("mask")) timeDispMask=atoi(child.attribute("mask").toAscii());
		child = child.nextSibling().toElement();
	}
}

boost::posix_time::ptime GLViewer::getLastUserEvent(){return last_user_event;};


float SudoDEMCamera::zNear() const
{
  float z = distanceToSceneCenter() - zClippingCoefficient()*sceneRadius()*(1.f-2*cuttingDistance);

  // Prevents negative or null zNear values.
  const float zMin = zNearCoefficient() * zClippingCoefficient() * sceneRadius();
  if (z < zMin)
/*    switch (type())
      {
      case Camera::PERSPECTIVE  :*/ z = zMin; /*break;
      case Camera::ORTHOGRAPHIC : z = 0.0;  break;
      }*/
  return z;
}
