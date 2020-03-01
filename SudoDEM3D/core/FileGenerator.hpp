/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/base/Logging.hpp>

#include "Scene.hpp"
#include "ThreadWorker.hpp"

class FileGenerator: public Serializable
{
	protected:
		shared_ptr<Scene>	 scene;
	public:
		bool generateAndSave(const string& outFile, string& message);
	protected :
	//! Returns whether the generation was successful; message for user is in FileGenerator::message
	virtual bool generate(std::string& msg);

	void pyGenerate(const string& out);
	void pyLoad();

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(FileGenerator,Serializable,"Base class for scene generators, preprocessors.",
		/*attrs*/,
		/*ctor*/
		,
		.def("generate",&FileGenerator::pyGenerate,(boost::python::arg("out")),"Generate scene, save to given file")
		.def("load",&FileGenerator::pyLoad,"Generate scene, save to temporary file and load immediately");
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FileGenerator);


