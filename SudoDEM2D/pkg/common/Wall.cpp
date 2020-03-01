// © 2009 Václav Šmilauer <eudoxos@arcig.cz>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Aabb.hpp>

SUDODEM_PLUGIN((Wall)(Bo1_Wall_Aabb) (Fwall) (Bo1_Fwall_Aabb)
	#ifdef SUDODEM_OPENGL
		(Gl1_Wall) (Gl1_Fwall)
	#endif
	);


Wall::~Wall(){} // vtable
Fwall::~Fwall(){}


void Bo1_Wall_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se2r& se2, const Body* b){
	Wall* wall=static_cast<Wall*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	if(scene->isPeriodic && scene->cell->hasShear()) throw logic_error(__FILE__ "Walls not supported in sheared cell.");
	const Real& inf=std::numeric_limits<Real>::infinity();
	aabb->min=Vector2r(-inf,-inf); aabb->min[wall->axis]=se2.position[wall->axis];
	aabb->max=Vector2r( inf, inf); aabb->max[wall->axis]=se2.position[wall->axis];
}


CREATE_LOGGER(Fwall);



void Fwall::postLoad(Fwall&)
{
	// if this fails, it means someone did vertices push_back, but they are resized to 3 at Fwall initialization already
	// in the future, a fixed-size array should be used instead of vector<Vector3r> for vertices
	// this is prevented by sudodem::serialization now IIRC
	//if(vertices.size()!=2){ throw runtime_error(("Fwall must have exactly 2 vertices (not "+boost::lexical_cast<string>(vertices.size())+")").c_str()); }
	if(isnan(vertex1[0]) || isnan(vertex2[0])) return;  // not initialized, nothing to do
	vu = vertex2-vertex1;
	vl = vu.norm();
	normal = Vector2r(-vu[1],vu[0]);
	normal.normalize();
}


void Bo1_Fwall_Aabb::go(	  const shared_ptr<Shape>& cm
				, shared_ptr<Bound>& bv
				, const Se2r& se2
				, const Body* b)
{
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	Fwall* fwall = static_cast<Fwall*>(cm.get());
	const Vector2r& O = se2.position;
	Matrix2r rot =se2.rotation.toRotationMatrix();
	Vector2r vu_new = rot*fwall->vu;
	Vector2r halfsize = 0.5*Vector2r(abs(vu_new[0]),abs(vu_new[1]));
	if(scene->isPeriodic){
			vu_new =scene->cell->unshearPt(O);
			aabb->min = vu_new - halfsize;
			aabb->max = vu_new + halfsize;
	}else{
		aabb->min = O - halfsize;
		aabb->max = O + halfsize;
	}

}



#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	int  Gl1_Wall::div=20;

	void Gl1_Wall::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& pp, bool, const GLViewInfo& glinfo){
		Wall* wall=static_cast<Wall*>(cm.get());
		int ax0=wall->axis,ax1=(wall->axis+1)%2;
		Vector2r a1,b1,a2,b2; // beginnings (a) and endings (b) of lines in both senses (0,1)
		// compensate for our position, since the functor is called with transformation to the wall se3 already, but we really want to be centered in the middle of the scene
		Real mn1=glinfo.sceneCenter[ax1]-glinfo.sceneRadius-pp->se2.position[ax1];
		//Real mn2=glinfo.sceneCenter[ax2]-glinfo.sceneRadius-pp->se3.position[ax2];
		Real step=2*glinfo.sceneRadius/div;
		//cerr<<"center "<<glinfo.sceneCenter<<", radius "<<glinfo.sceneRadius<<", mn["<<ax1<<"]="<<mn1<<", mn["<<ax2<<"]="<<mn2<<endl;

		a1[ax0]=b1[ax0]=0;
		a1[ax1]=mn1-step;
		b1[ax1]=mn1+step*(div+2);
		glColor3v(cm->color);
		//cout<<"color-"<<cm->color<<endl;
		glBegin(GL_LINES);
				glVertex2v(a1); glVertex2v(b1);
		glEnd();
	}

	void Gl1_Fwall::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& ,bool wire,const GLViewInfo&)
	{
		Fwall* fwall = static_cast<Fwall*>(cm.get());
		//const vector<Vector2r>& vertices = fwall->vertices;
		//if(cm->wire || wire){
			// Fwall
			glBegin(GL_LINES);
				glColor3v(cm->color);
				 glVertex2v(fwall->vertex1);
				 glVertex2v(fwall->vertex2);
				glEnd();
		//}
	}
#endif
