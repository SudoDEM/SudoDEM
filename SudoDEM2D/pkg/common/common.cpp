// =======================================================
// Some plugins from removed CPP-fiels

//#include <sudodem/pkg/dem/DemXDofGeom.hpp>
#include <sudodem/pkg/common/TorqueEngine.hpp>
#include <sudodem/pkg/common/ForceResetter.hpp>
#include <sudodem/pkg/common/FieldApplier.hpp>
#include <sudodem/pkg/common/Callbacks.hpp>
#include <sudodem/pkg/common/BoundaryController.hpp>
#include <sudodem/pkg/common/NormShearPhys.hpp>
#include <sudodem/pkg/common/Recorder.hpp>
//#include <sudodem/pkg/common/CylScGeom6D.hpp>
//#include <sudodem/pkg/common/Box.hpp>
#include <sudodem/pkg/common/StepDisplacer.hpp>
#include <sudodem/pkg/common/PeriodicEngines.hpp>
#include <sudodem/pkg/common/ElastMat.hpp>
#include <sudodem/pkg/common/PyRunner.hpp>
#include <sudodem/pkg/common/Disk.hpp>
#include <sudodem/pkg/common/Aabb.hpp>

SUDODEM_PLUGIN((IntrCallback)
	#ifdef SUDODEM_BODY_CALLBACK
		(BodyCallback)
	#endif
);

SUDODEM_PLUGIN((ForceResetter)(TorqueEngine)(FieldApplier)(BoundaryController)
            (NormPhys)(NormShearPhys)(Recorder)
            (StepDisplacer)
            (PeriodicEngine)(Disk)(Aabb)(ElastMat)(FrictMat)(PyRunner)
            );

#ifdef SUDODEM_OPENGL
#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
//#include <sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/GLDrawFunctors.hpp>
SUDODEM_PLUGIN((GlBoundFunctor)(GlShapeFunctor)(GlIGeomFunctor)(GlIPhysFunctor)(GlStateFunctor)
            (GlBoundDispatcher)(GlShapeDispatcher)(GlIGeomDispatcher)(GlIPhysDispatcher)
            (GlStateDispatcher));
#endif
