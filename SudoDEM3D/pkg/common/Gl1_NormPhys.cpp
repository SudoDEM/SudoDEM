#ifdef SUDODEM_OPENGL

#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/Gl1_NormPhys.hpp>
#include<sudodem/pkg/common/OpenGLRenderer.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/pkg/dem/Shop.hpp>


	SUDODEM_PLUGIN((Gl1_NormPhys));

	GLUquadric* Gl1_NormPhys::gluQuadric=NULL;
	Real Gl1_NormPhys::maxFn;
	Real Gl1_NormPhys::refRadius;
	Real Gl1_NormPhys::maxRadius;
	int Gl1_NormPhys::signFilter;
	int Gl1_NormPhys::slices;
	int Gl1_NormPhys::stacks;

	Real Gl1_NormPhys::maxWeakFn;
	int Gl1_NormPhys::weakFilter;
	Real Gl1_NormPhys::weakScale;

	void Gl1_NormPhys::go(const shared_ptr<IPhys>& ip, const shared_ptr<Interaction>& i, const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool wireFrame){
		if(!gluQuadric){ gluQuadric=gluNewQuadric(); if(!gluQuadric) throw runtime_error("Gl1_NormPhys::go unable to allocate new GLUquadric object (out of memory?)."); }
		NormPhys* np=static_cast<NormPhys*>(ip.get());
		shared_ptr<IGeom> ig(i->geom); if(!ig) return; // changed meanwhile?
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(ig.get());
		//if(!geom) cerr<<"Gl1_NormPhys: IGeom is not a GenericSpheresContact, but a "<<ig->getClassName()<<endl;
		Real fnNorm=np->normalForce.dot(geom->normal);
		if((signFilter>0 && fnNorm<0) || (signFilter<0 && fnNorm>0)) return;
		int fnSign=fnNorm>0?1:-1;
		fnNorm=std::abs(fnNorm);
		Real radiusScale=1.;
		// weak/strong fabric, only used if maxWeakFn is set
		if(!isnan(maxWeakFn)){
			if(fnNorm*fnSign<maxWeakFn){ // weak fabric
				if(weakFilter>0) return;
				radiusScale=weakScale;
			} else { // strong fabric
				if(weakFilter<0) return;
			}
		}

		maxFn=max(fnNorm,maxFn);
		Real realMaxRadius;
		if(maxRadius<0){
			if(geom->refR1>0) refRadius=min(geom->refR1,refRadius);
			if(geom->refR2>0) refRadius=min(geom->refR2,refRadius);
			realMaxRadius=refRadius;
		}
		else realMaxRadius=maxRadius;
		Real radius=radiusScale*realMaxRadius*(fnNorm/maxFn); // use logarithmic scale here?
		Vector3r color=Shop::scalarOnColorScale(fnNorm*fnSign,-maxFn,maxFn);
		# if 0
			// get endpoints from body positions
			Vector3r p1=b1->state->pos, p2=b2->state->pos;
			Vector3r relPos;
			if(scene->isPeriodic){
				relPos=p2+scene->cell->Hsize*i->cellDist.cast<Real>()-p1;
				p1=scene->cell->wrapShearedPt(p1);
				p2=p1+relPos;
			} else {
				relPos=p2-p1;
			}
			Real dist=relPos.norm();
		#else
			// get endpoints from geom
			// max(r,0) handles r<0 which is the case for "radius" of the facet
			Vector3r cp=scene->isPeriodic? scene->cell->wrapShearedPt(geom->contactPoint) : geom->contactPoint;
			Vector3r p1=cp-max(geom->refR1,(Real)0.)*geom->normal;
			Vector3r p2=cp+max(geom->refR2,(Real)0.)*geom->normal;
			const Vector3r& dispScale=scene->renderer ? scene->renderer->dispScale : Vector3r::Ones();
			if(dispScale!=Vector3r::Ones()){
				// move p1 and p2 by the same amounts as particles themselves would be moved
				p1+=dispScale.cwiseProduct(Vector3r(b1->state->pos-b1->state->refPos));
				p2+=dispScale.cwiseProduct(Vector3r(b2->state->pos-b2->state->refPos));
			}
			Vector3r relPos=p2-p1;
			Real dist=relPos.norm(); //max(geom->refR1,0.)+max(geom->refR2,0.);
	#endif


	glDisable(GL_CULL_FACE);
	glPushMatrix();
		glTranslatef(p1[0],p1[1],p1[2]);
		Quaternionr q(Quaternionr().setFromTwoVectors(Vector3r(0,0,1),relPos/dist /* normalized */));
		// using Transform with OpenGL: http://eigen.tuxfamily.org/dox/TutorialGeometry.html
		//glMultMatrixd(Eigen::Affine3d(q).data());
		glMultMatrix(Eigen::Transform<Real,3,Eigen::Affine>(q).data());
		glColor3v(color);
		gluCylinder(gluQuadric,radius,radius,dist,slices,stacks);
	glPopMatrix();
}

#endif /* SUDODEM_OPENGL */
