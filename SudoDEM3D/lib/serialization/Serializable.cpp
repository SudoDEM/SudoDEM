/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Serializable.hpp"


void Serializable::pyRegisterClass(boost::python::object _scope) {
	checkPyClassRegistersItself("Serializable");
	boost::python::scope thisScope(_scope);
	boost::python::class_<Serializable, shared_ptr<Serializable>, boost::noncopyable >("Serializable")
		.def("__str__",&Serializable::pyStr).def("__repr__",&Serializable::pyStr)
		.def("dict",&Serializable::pyDict,"Return dictionary of attributes.")
		.def("updateAttrs",&Serializable::pyUpdateAttrs,"Update object attributes from given dictionary")
		#if 1
			/* boost::python pickling support, as per http://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/pickle.html */
			.def("__getstate__",&Serializable::pyDict).def("__setstate__",&Serializable::pyUpdateAttrs)
			.add_property("__safe_for_unpickling__",&Serializable::getClassName,"just define the attr, return some bogus data")
			.add_property("__getstate_manages_dict__",&Serializable::getClassName,"just define the attr, return some bogus data")
		#endif
		// constructor with dictionary of attributes
		.def("__init__",boost::python::raw_constructor(Serializable_ctor_kwAttrs<Serializable>))
		// comparison operators
		.def(boost::python::self == boost::python::self)
		.def(boost::python::self != boost::python::self)
		;
}

void Serializable::checkPyClassRegistersItself(const std::string& thisClassName) const {
	if(getClassName()!=thisClassName) throw std::logic_error(("Class "+getClassName()+" does not register with SUDODEM_CLASS_BASE_DOC_ATTR*, would not be accessible from python.").c_str());
}


void Serializable::pyUpdateAttrs(const boost::python::dict& d){
	boost::python::list l=d.items(); size_t ll=boost::python::len(l); if(ll==0) return;
	for(size_t i=0; i<ll; i++){
		boost::python::tuple t=boost::python::extract<boost::python::tuple>(l[i]);
		string key=boost::python::extract<string>(t[0]);
		pySetAttr(key,t[1]);
	}
	callPostLoad();
}
