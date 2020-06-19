#include <sudodem/pkg/common/Disk.hpp>

SUDODEM_PLUGIN((Bo1_Disk_Aabb));

void Bo1_Disk_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se2r& se2, const Body* b){
	Disk* disk = static_cast<Disk*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	Vector2r halfSize = (aabbEnlargeFactor>0?aabbEnlargeFactor:1.)*Vector2r(disk->radius,disk->radius);
	if(!scene->isPeriodic){
		aabb->min=se2.position-halfSize; aabb->max=se2.position+halfSize;
		return;
	}
	// adjust box size along axes so that disk doesn't stick out of the box even if sheared (i.e. parallelepiped)
  /*
	if(scene->cell->hasShear()) {//debug for 2d
		Vector2r refHalfSize(halfSize);
		const Vector3r& cos=scene->cell->getCos();
		for(int i=0; i<3; i++){
			//cerr<<"cos["<<i<<"]"<<cos[i]<<" ";
			int i1=(i+1)%3,i2=(i+2)%3;
			halfSize[i1]+=.5*refHalfSize[i1]*(1/cos[i]-1);
			halfSize[i2]+=.5*refHalfSize[i2]*(1/cos[i]-1);
		}
	}*/
	if(scene->cell->hasShear()) {//debug for 2d
		Vector2r refHalfSize(halfSize);
		const Real& _cos=scene->cell->getCos();
			//cerr<<"cos["<<i<<"]"<<cos[i]<<" ";
			//halfSize+=.5*refHalfSize*(1/_cos-1);
			halfSize = refHalfSize/_cos;
		}
	//cerr<<" || "<<halfSize<<endl;
	aabb->min = scene->cell->unshearPt(se2.position)-halfSize;
	aabb->max = scene->cell->unshearPt(se2.position)+halfSize;
}


#ifdef SUDODEM_OPENGL
#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
#include <sudodem/core/Scene.hpp>
SUDODEM_PLUGIN((Gl1_Disk));

bool Gl1_Disk::wire;

int  Gl1_Disk::Slices;
int  Gl1_Disk::preSlices=5;

int Gl1_Disk::glStripedDiskList=-1;
int Gl1_Disk::glGlutDiskList=-1;

void Gl1_Disk::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& ,bool wire2, const GLViewInfo&)
{
	glClearDepth(1.0f);
	glEnable(GL_NORMALIZE);

	Real r=(static_cast<Disk*>(cm.get()))->radius;
	glColor3v(cm->color);

	//Check if quality has been modified or if previous lists are invalidated (e.g. by creating a new qt view), then regenerate lists
	bool somethingChanged = (Slices != preSlices || glIsList(glStripedDiskList)!=GL_TRUE);
	if (somethingChanged) {
		initStripedGlList();
		initGlutGlList();
		if(Slices<5){Slices = 5;}
		preSlices=Slices;
	}
	glScalef(r,r,r);
	if(wire||wire2) glCallList(glGlutDiskList);
	else glCallList(glStripedDiskList);

	return;
}


void Gl1_Disk::initStripedGlList() {

	//Generate the list. Only once for each qtView, or more if quality is modified.
	glDeleteLists(glStripedDiskList,1);
	glStripedDiskList = glGenLists(1);
	glNewList(glStripedDiskList,GL_COMPILE);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	// render the disk now
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(1.,0);
  double del_ang = Mathr::TWO_PI/Slices;
	float angle=0.0;
	for (int i=1 ;i<Slices;i++)
	{		angle = del_ang*i;
	    glVertex2f(cos(angle),sin(angle));
	}
	glVertex2f(1.,0);

	glEnd();
	glEndList();

}

void Gl1_Disk::initGlutGlList(){
	//Generate the "no-stripes" display list, each time quality is modified
	glDeleteLists(glGlutDiskList,1);
	glGlutDiskList = glGenLists(1);
	glNewList(glGlutDiskList,GL_COMPILE);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		glBegin(GL_LINE_LOOP);
		glVertex2f(1.,0);
	  double del_ang = Mathr::TWO_PI/Slices;
		float angle=0.0;
		for (int i=1 ;i<Slices;i++)
		{		angle = del_ang*i;
		    glVertex2f(cos(angle),sin(angle));
		}
		//glVertex2f(1.,0);

		glEnd();
	glEndList();
}
#endif
