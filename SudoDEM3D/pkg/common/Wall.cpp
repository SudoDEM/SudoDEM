// © 2009 Václav Šmilauer <eudoxos@arcig.cz>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Aabb.hpp>

SUDODEM_PLUGIN((Wall)(Bo1_Wall_Aabb)
	#ifdef SUDODEM_OPENGL
		(Gl1_Wall)
	#endif
	);

Wall::~Wall(){} // vtable

void Bo1_Wall_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* b){
	Wall* wall=static_cast<Wall*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	if(scene->isPeriodic && scene->cell->hasShear()) throw logic_error(__FILE__ "Walls not supported in sheared cell.");
	const Real& inf=std::numeric_limits<Real>::infinity();
	aabb->min=Vector3r(-inf,-inf,-inf); aabb->min[wall->axis]=se3.position[wall->axis];
	aabb->max=Vector3r( inf, inf, inf); aabb->max[wall->axis]=se3.position[wall->axis];
}


#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	int  Gl1_Wall::div=20;

	void Gl1_Wall::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& pp, bool, const GLViewInfo& glinfo){
		Wall* wall=static_cast<Wall*>(cm.get());
		int ax0=wall->axis,ax1=(wall->axis+1)%3,ax2=(wall->axis+2)%3;
		Vector3r a1,b1,a2,b2; // beginnings (a) and endings (b) of lines in both senses (0,1)
		// compensate for our position, since the functor is called with transformation to the wall se3 already, but we really want to be centered in the middle of the scene
		Real mn1=glinfo.sceneCenter[ax1]-glinfo.sceneRadius-pp->se3.position[ax1];
		Real mn2=glinfo.sceneCenter[ax2]-glinfo.sceneRadius-pp->se3.position[ax2];
		Real step=2*glinfo.sceneRadius/div;
		//cerr<<"center "<<glinfo.sceneCenter<<", radius "<<glinfo.sceneRadius<<", mn["<<ax1<<"]="<<mn1<<", mn["<<ax2<<"]="<<mn2<<endl;

		a1[ax0]=b1[ax0]=a2[ax0]=b2[ax0]=0;
		a1[ax1]=mn1-step; a2[ax2]=mn2-step;
		b1[ax1]=mn1+step*(div+2); b2[ax2]=mn2+step*(div+2);
		glColor3v(cm->color);
		glBegin(GL_LINES);
			for(int i=0; i<=div; i++){
				a1[ax2]=b1[ax2]=mn1+i*step;
				a2[ax1]=b2[ax1]=mn2+i*step;
				glVertex3v(a1); glVertex3v(b1);
				glVertex3v(a2); glVertex3v(b2);
			}
		glEnd();
	}
#endif

