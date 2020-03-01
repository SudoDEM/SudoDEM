#include <sudodem/pkg/common/Aabb.hpp>


#ifdef SUDODEM_OPENGL
#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
#include <sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((Gl1_Aabb));
		void Gl1_Aabb::go(const shared_ptr<Bound>& bv, Scene* scene){
			Aabb* aabb = static_cast<Aabb*>(bv.get());
			glColor3v(bv->color);
			if(!scene->isPeriodic){
				glTranslatev(Vector3r(.5*(aabb->min+aabb->max)));
				glScalev(Vector3r(aabb->max-aabb->min));
			} else {
				glTranslatev(Vector3r(scene->cell->shearPt(scene->cell->wrapPt(.5*(aabb->min+aabb->max)))));
				glMultMatrixd(scene->cell->getGlShearTrsfMatrix());
				glScalev(Vector3r(aabb->max-aabb->min));
			}
			glutWireCube(1);
		}
#endif
