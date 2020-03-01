#include <sudodem/pkg/common/Box.hpp>

SUDODEM_PLUGIN((Bo1_Box_Aabb));


void Bo1_Box_Aabb::go(	const shared_ptr<Shape>& cm,
				shared_ptr<Bound>& bv,
				const Se3r& se3,
				const Body*	b)
{
	Box* box = static_cast<Box*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());

	if(scene->isPeriodic && scene->cell->hasShear()) throw logic_error(__FILE__ "Boxes not (yet?) supported in sheared cell.");

	Matrix3r r=se3.orientation.toRotationMatrix();
	Vector3r halfSize(Vector3r::Zero());
	for( int i=0; i<3; ++i )
		for( int j=0; j<3; ++j )
			halfSize[i] += std::abs( r(i,j) * box->extents[j] );

	aabb->min = se3.position-halfSize;
	aabb->max = se3.position+halfSize;
}

#ifdef SUDODEM_OPENGL

#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
#include <sudodem/core/Scene.hpp>
SUDODEM_PLUGIN((Gl1_Box));

void Gl1_Box::go(const shared_ptr<Shape>& cg, const shared_ptr<State>&,bool wire,const GLViewInfo&)
{
	glColor3v(cg->color);
	Vector3r &extents = (static_cast<Box*>(cg.get()))->extents;
	glScalef(2*extents[0],2*extents[1],2*extents[2]);
	if (wire) glutWireCube(1);
	else glutSolidCube(1);
}
#endif
