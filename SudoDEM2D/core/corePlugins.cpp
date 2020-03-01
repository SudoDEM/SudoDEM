#include<sudodem/lib/factory/ClassFactory.hpp>
// make core classes known to the class factory
#include<sudodem/core/Body.hpp>
#include<sudodem/core/BodyContainer.hpp>
#include<sudodem/core/Bound.hpp>
#include<sudodem/core/Cell.hpp>
#include<sudodem/core/Dispatcher.hpp>
#include<sudodem/core/EnergyTracker.hpp>
#include<sudodem/core/Engine.hpp>
#include<sudodem/core/FileGenerator.hpp>
#include<sudodem/core/Functor.hpp>
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/InteractionContainer.hpp>
#include<sudodem/core/IGeom.hpp>
#include<sudodem/core/IPhys.hpp>
#include<sudodem/core/Material.hpp>
#include<sudodem/core/PartialEngine.hpp>
#include<sudodem/core/Shape.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/core/TimeStepper.hpp>


// these two are not accessible from python directly (though they should be in the future, perhaps)

#if BOOST_VERSION>=104200
	BOOST_CLASS_EXPORT_IMPLEMENT(BodyContainer);
	BOOST_CLASS_EXPORT_IMPLEMENT(InteractionContainer);
#else
	BOOST_CLASS_EXPORT(BodyContainer);
	BOOST_CLASS_EXPORT(InteractionContainer);
#endif

SUDODEM_PLUGIN((Body)(Bound)(Cell)(Dispatcher)(EnergyTracker)(Engine)(FileGenerator)(Functor)(GlobalEngine)(Interaction)(IGeom)(IPhys)(Material)(PartialEngine)(Shape)(State)(TimeStepper));

EnergyTracker::~EnergyTracker(){} // vtable

//BOOST_CLASS_EXPORT(OpenMPArrayAccumulator<Real>);
//BOOST_CLASS_EXPORT(OpenMPAccumulator<Real>);
