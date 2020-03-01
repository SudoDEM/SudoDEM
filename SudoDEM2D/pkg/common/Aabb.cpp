#include <sudodem/pkg/common/Aabb.hpp>


#ifdef SUDODEM_OPENGL
#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/lib/opengl/GLUtils.hpp>
#include <sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((Gl1_Aabb));
		void Gl1_Aabb::go(const shared_ptr<Bound>& bv, Scene* scene){
			Aabb* aabb = static_cast<Aabb*>(bv.get());
			glColor3v(bv->color);

			if(!scene->isPeriodic){
				glTranslate2v(Vector2r(.5*(aabb->min+aabb->max)));
				glScale2v(Vector2r(aabb->max-aabb->min));
				//glBegin(GL_LINE_STRIP);
		 	 	//glVertex2v(aabb->min); glVertex2(aabb->max[0],aabb->min[1]); glVertex2v(aabb->max);
		    //glVertex2(aabb->min[0],aabb->max[1]); glVertex2v(aabb->min);
		 	  //glEnd();
			} else {
				glTranslate2v(Vector2r(scene->cell->shearPt(scene->cell->wrapPt(.5*(aabb->min+aabb->max)))));
				glMultMatrixd(scene->cell->getGlShearTrsfMatrix());
				glScale2v(Vector2r(aabb->max-aabb->min));//FIXME:not working without reconstrution of Vector2r

				//glutWireCube(1);
				/*Vector2r min,max;
				Vector2r trans = Vector2r(scene->cell->shearPt(scene->cell->wrapPt(.5*(aabb->min+aabb->max))));
				min = aabb->min + trans; max = aabb->max +trans;
				glBegin(GL_LINE_STRIP);
		 	 	glVertex2v(min); glVertex2(max[0],min[1]); glVertex2v(max);
		     glVertex2(min[0],max[1]); glVertex2v(min);
		 	  glEnd();*/
			}
			GLUtils::Square(1);

		}
#endif
