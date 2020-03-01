#pragma once
#include <sudodem/core/GlobalEngine.hpp>
#include <stdexcept>

class FieldApplier: public GlobalEngine{
	virtual void action() {
		throw std::runtime_error("FieldApplier must not be used in simulations directly (FieldApplier::action called).");
	}
	SUDODEM_CLASS_BASE_DOC_ATTRS(FieldApplier,GlobalEngine,"Base for engines applying force files on particles. Not to be used directly.",
		((int,fieldWorkIx,-1,(Attr::hidden|Attr::noSave),"Index for the work done by this field, if tracking energies."))
	);
};
REGISTER_SERIALIZABLE(FieldApplier);

