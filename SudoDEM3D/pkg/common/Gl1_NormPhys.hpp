#pragma once

#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/lib/opengl/GLUtils.hpp>
#include<sudodem/pkg/common/GLDrawFunctors.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<GL/glu.h>


class Gl1_NormPhys: public GlIPhysFunctor{
		static GLUquadric* gluQuadric; // needed for gluCylinder, initialized by ::go if no initialized yet
	public:
		virtual void go(const shared_ptr<IPhys>&,const shared_ptr<Interaction>&,const shared_ptr<Body>&,const shared_ptr<Body>&,bool wireFrame);
	SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_NormPhys,GlIPhysFunctor,"Renders :yref:`NormPhys` objects as cylinders of which diameter and color depends on :yref:`NormPhys.normalForce` magnitude.",
		// changed doc maxDiameter -> maxRadius ((Real,maxFn,0,,"Value of :yref:`NormPhys.normalForce` corresponding to :yref:`maxDiameter<Gl1_NormPhys.maxDiameter>`. This value will be increased (but *not decreased* ) automatically."))
		((Real,maxFn,0,,"Value of :yref:`NormPhys.normalForce` corresponding to :yref:`maxRadius<Gl1_NormPhys.maxRadius>`. This value will be increased (but *not decreased* ) automatically."))
		((int,signFilter,0,,"If non-zero, only display contacts with negative (-1) or positive (+1) normal forces; if zero, all contacts will be displayed."))
		((Real,refRadius,std::numeric_limits<Real>::infinity(),,"Reference (minimum) particle radius; used only if :yref:`maxRadius<Gl1_NormPhys.maxRadius>` is negative. This value will be decreased (but *not increased* ) automatically. |yupdate|"))
		((Real,maxRadius,-1,,"Cylinder radius corresponding to the maximum normal force. If negative, auto-updated :yref:`refRadius<Gl1_NormPhys.refRadius>` will be used instead."))
		((int,slices,6,,"Number of sphere slices; (see `glutCylinder reference <http://www.opengl.org/sdk/docs/man/xhtml/gluCylinder.xml>`_)")) // FIXME: the link does not exist
		((int,stacks,1,,"Number of sphere stacks; (see `glutCylinder reference <http://www.opengl.org/sdk/docs/man/xhtml/gluCylinder.xml>`_)"))
		// strong/weak fabric attributes
		((Real,maxWeakFn,NaN,,"Value that divides contacts by their normal force into the 'weak fabric' and 'strong fabric'. This value is set as side-effect by :yref:`utils.fabricTensor`."))
		((int,weakFilter,0,,"If non-zero, only display contacts belonging to the 'weak' (-1) or 'strong' (+1) fabric."))
		((Real,weakScale,1.,,"If :yref:`maxWeakFn<Gl1_NormPhys.maxWeakFn>` is set, scale radius of the weak fabric by this amount (usually smaller than 1). If zero, 1 pixel line is displayed. Colors are not affected by this value."))
	);
	RENDERS(NormPhys);
};
REGISTER_SERIALIZABLE(Gl1_NormPhys);

