// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>

#ifdef SUDODEM_OPENGL

#include"OpenGLRenderer.hpp"
#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/lib/opengl/GLUtils.hpp>
#include<sudodem/core/Timing.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/Aabb.hpp>
#include<sudodem/lib/pyutil/gil.hpp>

#ifdef __APPLE__
#  include <OpenGL/glu.h>
#  include <OpenGL/gl.h>
#  include <GLUT/glut.h>
#else
#  include <GL/glu.h>
#  include <GL/gl.h>
#  include <GL/glut.h>
#endif

SUDODEM_PLUGIN((OpenGLRenderer)(GlExtraDrawer));
CREATE_LOGGER(OpenGLRenderer);

void GlExtraDrawer::render(){ throw runtime_error("GlExtraDrawer::render called from class "+getClassName()+". (did you forget to override it in the derived class?)"); }

bool OpenGLRenderer::initDone=false;
const int OpenGLRenderer::numClipPlanes;
OpenGLRenderer::~OpenGLRenderer(){}


void OpenGLRenderer::init(){
	typedef std::pair<string,DynlibDescriptor> strDldPair; // necessary as FOREACH, being macro, cannot have the "," inside the argument (preprocessor does not parse templates)
	FOREACH(const strDldPair& item, Omega::instance().getDynlibsDescriptor()){
		// if (Omega::instance().isInheritingFrom_recursive(item.first,"GlStateFunctor")) stateFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first,"GlBoundFunctor")) boundFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first,"GlShapeFunctor")) shapeFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first,"GlIGeomFunctor")) geomFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first,"GlIPhysFunctor")) physFunctorNames.push_back(item.first);
	}
	initgl(); // creates functor objects in the proper sense

	clipPlaneNormals.resize(numClipPlanes);

	static bool glutInitDone=false;
	if(!glutInitDone){
		glutInit(&Omega::instance().origArgc,Omega::instance().origArgv);
		/* transparent spheres (still not working): glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE | GLUT_ALPHA); glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE); */
		glutInitDone=true;
	}

	initDone=true;
	// glGetError crashes at some machines?! Was never really useful, anyway.
	// reported http://www.mail-archive.com/sudodem-users@lists.launchpad.net/msg01482.html
	#if 0
		int e=glGetError();
		if(e!=GL_NO_ERROR) throw runtime_error((string("OpenGLRenderer::init returned GL error ")+boost::lexical_cast<string>(e)).c_str());
	#endif
}

void OpenGLRenderer::setBodiesRefSe3(){
	LOG_DEBUG("(re)initializing reference positions and orientations.");
	FOREACH(const shared_ptr<Body>& b, *scene->bodies) if(b && b->state) { b->state->refPos=b->state->pos; b->state->refOri=b->state->ori; }
	scene->cell->refHSize=scene->cell->hSize;
}


void OpenGLRenderer::initgl(){
	LOG_DEBUG("(re)initializing GL for gldraw methods.\n");
	#define _SETUP_DISPATCHER(names,FunctorType,dispatcher) dispatcher.clearMatrix(); FOREACH(string& s,names) {shared_ptr<FunctorType> f(boost::static_pointer_cast<FunctorType>(ClassFactory::instance().createShared(s))); f->initgl(); dispatcher.add(f);}
		// _SETUP_DISPATCHER(stateFunctorNames,GlStateFunctor,stateDispatcher);
		_SETUP_DISPATCHER(boundFunctorNames,GlBoundFunctor,boundDispatcher);
		_SETUP_DISPATCHER(shapeFunctorNames,GlShapeFunctor,shapeDispatcher);
		_SETUP_DISPATCHER(geomFunctorNames,GlIGeomFunctor,geomDispatcher);
		_SETUP_DISPATCHER(physFunctorNames,GlIPhysFunctor,physDispatcher);
	#undef _SETUP_DISPATCHER
}

bool OpenGLRenderer::pointClipped(const Vector3r& p){
	if(numClipPlanes<1) return false;
	for(int i=0;i<numClipPlanes;i++) if(clipPlaneActive[i]&&(p-clipPlaneSe3[i].position).dot(clipPlaneNormals[i])<0) return true;
	return false;
}


void OpenGLRenderer::setBodiesDispInfo(){
	if(scene->bodies->size()!=bodyDisp.size()) {
		bodyDisp.resize(scene->bodies->size());
		for (unsigned k=0; k<scene->bodies->size(); k++) bodyDisp[k].hidden=0;}
	bool scaleRotations=(rotScale!=1.0);
	bool scaleDisplacements=(dispScale!=Vector3r::Ones());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->state) continue;
		size_t id=b->getId();
		const Vector3r& pos=b->state->pos; const Vector3r& refPos=b->state->refPos;
		const Quaternionr& ori=b->state->ori; const Quaternionr& refOri=b->state->refOri;
		Vector3r cellPos=(!scene->isPeriodic ? pos : scene->cell->wrapShearedPt(pos)); // inside the cell if periodic, same as pos otherwise
		bodyDisp[id].isDisplayed=!pointClipped(cellPos);
		// if no scaling and no periodic, return quickly
		if(!(scaleDisplacements||scaleRotations||scene->isPeriodic)){ bodyDisp[id].pos=pos; bodyDisp[id].ori=ori; continue; }
		// apply scaling
		bodyDisp[id].pos=cellPos; // point of reference (inside the cell for periodic)
		if(scaleDisplacements) bodyDisp[id].pos+=dispScale.cwiseProduct(Vector3r(pos-refPos)); // add scaled translation to the point of reference
		if(!scaleRotations) bodyDisp[id].ori=ori;
		else{
			Quaternionr relRot=refOri.conjugate()*ori;
			AngleAxisr aa(relRot);
			aa.angle()*=rotScale;
			bodyDisp[id].ori=refOri*Quaternionr(aa);
		}
	}
}

// draw periodic cell, if active
void OpenGLRenderer::drawPeriodicCell(){
	if(!scene->isPeriodic) return;
	glColor3v(cellColor);
	glPushMatrix();
		// Vector3r size=scene->cell->getSize();
		const Matrix3r& hSize=scene->cell->hSize;
		if(dispScale!=Vector3r::Ones()){
			const Matrix3r& refHSize(scene->cell->refHSize);
			Matrix3r scaledHSize;
			for(int i=0; i<3; i++) scaledHSize.col(i)=refHSize.col(i)+dispScale.cwiseProduct(Vector3r(hSize.col(i)-refHSize.col(i)));
			GLUtils::Parallelepiped(scaledHSize.col(0),scaledHSize.col(1),scaledHSize.col(2));
		} else {
			GLUtils::Parallelepiped(hSize.col(0),hSize.col(1),hSize.col(2));
		}
	glPopMatrix();
}

void OpenGLRenderer::resetSpecularEmission(){
	glMateriali(GL_FRONT, GL_SHININESS, 80);
	const GLfloat glutMatSpecular[4]={0.3,0.3,0.3,0.5};
	const GLfloat glutMatEmit[4]={0.2,0.2,0.2,1.0};
	glMaterialfv(GL_FRONT,GL_SPECULAR,glutMatSpecular);
	glMaterialfv(GL_FRONT,GL_EMISSION,glutMatEmit);
}

void OpenGLRenderer::render(const shared_ptr<Scene>& _scene,Body::id_t selection){

	gilLock lockgil;
	if(!initDone) init();
	assert(initDone);
	selId = selection;

	scene=_scene;

	// assign scene inside functors
	boundDispatcher.updateScenePtr();
	geomDispatcher.updateScenePtr();
	physDispatcher.updateScenePtr();
	shapeDispatcher.updateScenePtr();
	// stateDispatcher.updateScenePtr();

	// just to make sure, since it is not initialized by default
	if(!scene->bound) scene->bound=shared_ptr<Aabb>(new Aabb);

	// recompute emissive light colors for highlighted bodies
	Real now=TimingInfo::getNow(/*even if timing is disabled*/true)*1e-9;
	highlightEmission0[0]=highlightEmission0[1]=highlightEmission0[2]=.8*normSquare(now,1);
	highlightEmission1[0]=highlightEmission1[1]=highlightEmission0[2]=.5*normSaw(now,2);

	// clipping
	assert(clipPlaneNormals.size()==(size_t)numClipPlanes);
	for(size_t i=0;i<(size_t)numClipPlanes; i++){
		// someone could have modified those from python and truncate the vectors; fill those here in that case
		if(i==clipPlaneSe3.size()) clipPlaneSe3.push_back(Se3r(Vector3r::Zero(),Quaternionr::Identity()));
		if(i==clipPlaneActive.size()) clipPlaneActive.push_back(false);
		if(i==clipPlaneNormals.size()) clipPlaneNormals.push_back(Vector3r::UnitX());
		// end filling stuff modified from python
		if(clipPlaneActive[i]) clipPlaneNormals[i]=clipPlaneSe3[i].orientation*Vector3r(0,0,1);
		/* glBegin(GL_LINES);glVertex3v(clipPlaneSe3[i].position);glVertex3v(clipPlaneSe3[i].position+clipPlaneNormals[i]);glEnd(); */
	}
	// set displayed Se3 of body (scaling) and isDisplayed (clipping)
	setBodiesDispInfo();

	glClearColor(bgColor[0],bgColor[1],bgColor[2],1.0);

	// set light sources
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1); // important: do lighting calculations on both sides of polygons

	const GLfloat pos[4]	= {(float) lightPos[0], (float) lightPos[1], (float) lightPos[2],1.0};
	const GLfloat ambientColor[4]={0.2,0.2,0.2,1.0};
	const GLfloat specularColor[4]={1,1,1,1.f};
	const GLfloat diffuseLight[4] = { (float) lightColor[0], (float) lightColor[1], (float) lightColor[2], 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION,pos);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	if (light1) glEnable(GL_LIGHT0); else glDisable(GL_LIGHT0);

	const GLfloat pos2[4]	= {(float) light2Pos[0], (float) light2Pos[1], (float) light2Pos[2],1.0};
	const GLfloat ambientColor2[4]={0.0,0.0,0.0,1.0};
	const GLfloat specularColor2[4]={1,1,0.6,1.f};
	const GLfloat diffuseLight2[4] = { (float) light2Color[0], (float) light2Color[1], (float) light2Color[2], 1.0f };
	glLightfv(GL_LIGHT1, GL_POSITION,pos2);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specularColor2);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambientColor2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight2);
	if (light2) glEnable(GL_LIGHT1); else glDisable(GL_LIGHT1);

	glEnable(GL_LIGHTING);

	glEnable(GL_CULL_FACE);
	// http://www.sjbaker.org/steve/omniv/opengl_lighting.html
	glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	//Shared material settings
	resetSpecularEmission();

	drawPeriodicCell();

	if (dof || id) renderDOF_ID();
	if (bound) renderBound();

	if (shape) renderShape();
	if (intrAllWire) renderAllInteractionsWire();
	if (intrGeom) renderIGeom();
	if (intrPhys) renderIPhys();

	FOREACH(const shared_ptr<GlExtraDrawer> d, extraDrawers){
		if(d->dead) continue;
		glPushMatrix();
			d->scene=scene.get();
			d->render();
		glPopMatrix();
	}


}

void OpenGLRenderer::renderAllInteractionsWire(){
	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
		if(!i->functorCache.geomExists) continue;
		glColor3v(i->isReal()? Vector3r(0,1,0) : Vector3r(.5,0,1));
		Vector3r p1=Body::byId(i->getId1(),scene)->state->pos;
		const Vector3r& size=scene->cell->getSize();
		Vector3r shift2(i->cellDist[0]*size[0],i->cellDist[1]*size[1],i->cellDist[2]*size[2]);
		// in sheared cell, apply shear on the mutual position as well
		shift2=scene->cell->shearPt(shift2);
		Vector3r rel=Body::byId(i->getId2(),scene)->state->pos+shift2-p1;
		if(scene->isPeriodic) p1=scene->cell->wrapShearedPt(p1);
		glBegin(GL_LINES); glVertex3v(p1);glVertex3v(Vector3r(p1+rel));glEnd();
	}
}

void OpenGLRenderer::renderDOF_ID(){
	const GLfloat ambientColorSelected[4]={10.0,0.0,0.0,1.0};
	const GLfloat ambientColorUnselected[4]={0.5,0.5,0.5,1.0};
	FOREACH(const shared_ptr<Body> b, *scene->bodies){
		if(!b) continue;
		if(b->shape && ((b->getGroupMask() & mask) || b->getGroupMask()==0)){
			if(!id && b->state->blockedDOFs==0) continue;
			if(selId==b->getId()){glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColorSelected);}
			{ // write text
				glColor3f(1.0-bgColor[0],1.0-bgColor[1],1.0-bgColor[2]);
				unsigned d = b->state->blockedDOFs;
				std::string sDof = std::string()+(((d&State::DOF_X )!=0)?"x":"")+(((d&State::DOF_Y )!=0)?"y":" ")+(((d&State::DOF_Z )!=0)?"z":"")+(((d&State::DOF_RX)!=0)?"X":"")+(((d&State::DOF_RY)!=0)?"Y":"")+(((d&State::DOF_RZ)!=0)?"Z":"");
				std::string sId = boost::lexical_cast<std::string>(b->getId());
				std::string str;
				if(dof && id) sId += " ";
				if(id) str += sId;
				if(dof) str += sDof;
				const Vector3r& h(selId==b->getId() ? highlightEmission0 : Vector3r(1,1,1));
				glColor3v(h);
				GLUtils::GLDrawText(str,bodyDisp[b->id].pos,h);
			}
			if(selId == b->getId()){glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColorUnselected);}
		}
	}
}

void OpenGLRenderer::renderIGeom(){
	geomDispatcher.scene=scene.get(); geomDispatcher.updateScenePtr();
	{
		boost::mutex::scoped_lock lock(scene->interactions->drawloopmutex);
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
			if(!I->geom) continue; // avoid refcount manipulations if the interaction is not real anyway
			shared_ptr<IGeom> ig(I->geom); // keep reference so that ig does not disappear suddenly while being rendered
			if(!ig) continue;
			const shared_ptr<Body>& b1=Body::byId(I->getId1(),scene), b2=Body::byId(I->getId2(),scene);
			if(!(bodyDisp[I->getId1()].isDisplayed||bodyDisp[I->getId2()].isDisplayed)) continue;
			glPushMatrix(); geomDispatcher(ig,I,b1,b2,intrWire); glPopMatrix();
		}
	}
}


void OpenGLRenderer::renderIPhys(){
	physDispatcher.scene=scene.get(); physDispatcher.updateScenePtr();
	{
		boost::mutex::scoped_lock lock(scene->interactions->drawloopmutex);
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
			shared_ptr<IPhys> ip(I->phys);
			if(!ip) continue;
			const shared_ptr<Body>& b1=Body::byId(I->getId1(),scene), b2=Body::byId(I->getId2(),scene);
			Body::id_t id1=I->getId1(), id2=I->getId2();
			if(!(bodyDisp[id1].isDisplayed||bodyDisp[id2].isDisplayed)) continue;
			glPushMatrix(); physDispatcher(ip,I,b1,b2,intrWire); glPopMatrix();
		}
	}
}

void OpenGLRenderer::renderBound(){
	boundDispatcher.scene=scene.get(); boundDispatcher.updateScenePtr();

	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->bound) continue;
		if(!bodyDisp[b->getId()].isDisplayed or bodyDisp[b->getId()].hidden) continue;
		if(b->bound && ((b->getGroupMask()&mask) || b->getGroupMask()==0)){
			glPushMatrix(); boundDispatcher(b->bound,scene.get()); glPopMatrix();
		}
	}
	// since we remove the functor as Scene doesn't inherit from Body anymore, hardcore the rendering routine here
	// for periodic scene, renderPeriodicCell is called separately
	if(!scene->isPeriodic){
		if(!scene->bound) scene->updateBound();
		glColor3v(Vector3r(0,1,0));
		Vector3r size=scene->bound->max-scene->bound->min;
		Vector3r center=.5*(scene->bound->min+scene->bound->max);
		glPushMatrix();
			glTranslatev(center);
			glScalev(size);
			glutWireCube(1);
		glPopMatrix();
	}
}

// this function is called for both rendering as well as
// in the selection mode

// nice reading on OpenGL selection
// http://glprogramming.com/red/chapter13.html

void OpenGLRenderer::renderShape(){
	shapeDispatcher.scene=scene.get(); shapeDispatcher.updateScenePtr();

	// instead of const shared_ptr&, get proper shared_ptr;
	// Less efficient in terms of performance, since memory has to be written (not measured, though),
	// but it is still better than crashes if the body gets deleted meanwile.
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if(!b || !b->shape) continue;
		if(!(bodyDisp[b->getId()].isDisplayed and !bodyDisp[b->getId()].hidden)) continue;
		Vector3r pos=bodyDisp[b->getId()].pos;
		Quaternionr ori=bodyDisp[b->getId()].ori;
		if(!b->shape || !((b->getGroupMask()&mask) || b->getGroupMask()==0)) continue;

		// ignored in non-selection mode, use it always
		glPushName(b->id);
		bool highlight=(b->id==selId || (b->clumpId>=0 && b->clumpId==selId) || b->shape->highlight);

		glPushMatrix();
			AngleAxisr aa(ori);
			glTranslatef(pos[0],pos[1],pos[2]);
			glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
			if(highlight){
				// set hightlight
				// different color for body highlighted by selection and by the shape attribute
				const Vector3r& h((selId==b->id||(b->clumpId>=0 && selId==b->clumpId)) ? highlightEmission0 : highlightEmission1);
				glMaterialv(GL_FRONT_AND_BACK,GL_EMISSION,h);
				glMaterialv(GL_FRONT_AND_BACK,GL_SPECULAR,h);
				shapeDispatcher(b->shape,b->state,wire || b->shape->wire,viewInfo);
				// reset highlight
				resetSpecularEmission();
			} else {
				// no highlight; in case previous functor fiddled with glMaterial
				resetSpecularEmission();
				shapeDispatcher(b->shape,b->state,wire || b->shape->wire,viewInfo);
			}
		glPopMatrix();
		if(highlight){
			if(!b->bound || wire || b->shape->wire) GLUtils::GLDrawInt(b->getId(),pos);
			else {
				// move the label towards the camera by the bounding box so that it is not hidden inside the body
				const Vector3r& mn=b->bound->min; const Vector3r& mx=b->bound->max; const Vector3r& p=pos;
				Vector3r ext(viewDirection[0]>0?p[0]-mn[0]:p[0]-mx[0],viewDirection[1]>0?p[1]-mn[1]:p[1]-mx[1],viewDirection[2]>0?p[2]-mn[2]:p[2]-mx[2]); // signed extents towards the camera
				Vector3r dr=-1.01*(viewDirection.dot(ext)*viewDirection);
				GLUtils::GLDrawInt(b->getId(),pos+dr,Vector3r::Ones());
			}
		}
		// if the body goes over the cell margin, draw it in positions where the bbox overlaps with the cell in wire
		// precondition: pos is inside the cell.
		if(b->bound && scene->isPeriodic && ghosts){
			const Vector3r& cellSize(scene->cell->getSize());
			pos=scene->cell->unshearPt(pos); // remove the shear component
			// traverse all periodic cells around the body, to see if any of them touches
			Vector3r halfSize=b->bound->max-b->bound->min; halfSize*=.5;
			Vector3r pmin,pmax;
			Vector3i i;
			for(i[0]=-1; i[0]<=1; i[0]++) for(i[1]=-1;i[1]<=1; i[1]++) for(i[2]=-1; i[2]<=1; i[2]++){
				if(i[0]==0 && i[1]==0 && i[2]==0) continue; // middle; already rendered above
				Vector3r pos2=pos+Vector3r(cellSize[0]*i[0],cellSize[1]*i[1],cellSize[2]*i[2]); // shift, but without shear!
				pmin=pos2-halfSize; pmax=pos2+halfSize;
				if(pmin[0]<=cellSize[0] && pmax[0]>=0 &&
					pmin[1]<=cellSize[1] && pmax[1]>=0 &&
					pmin[2]<=cellSize[2] && pmax[2]>=0) {
					Vector3r pt=scene->cell->shearPt(pos2);
					if(pointClipped(pt)) continue;
					glLoadName(b->id);
					glPushMatrix();
						glTranslatev(pt);
						glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
						shapeDispatcher(b->shape,b->state,/*wire*/ true, viewInfo);
					glPopMatrix();
				}
			}
		}
		glPopName();
	}
}


#endif /* SUDODEM_OPENGL */
