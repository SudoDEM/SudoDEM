// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/pkg/common/PeriodicEngines.hpp>
class Recorder: public PeriodicEngine{
	void openAndCheck() {
		assert(!out.is_open());

		std::string fileTemp = file;
		if (addIterNum) fileTemp+="-" + boost::lexical_cast<string>(scene->iter);

		if(fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
		out.open(fileTemp.c_str(), truncate ? fstream::trunc : fstream::app);
		if(!out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `"+fileTemp+"'.");
	}
	protected:
		//! stream object that derived engines should write to
		std::ofstream out;
	public:
		virtual ~Recorder() {};
		virtual bool isActivated(){
			if(PeriodicEngine::isActivated()){
				if(!out.is_open()) openAndCheck();
				return true;
			}
			return false;
		}
	SUDODEM_CLASS_BASE_DOC_ATTRS(Recorder,PeriodicEngine,"Engine periodically storing some data to (one) external file. In addition PeriodicEngine, it handles opening the file as needed. See :yref:`PeriodicEngine` for controlling periodicity.",
		((std::string,file,,,"Name of file to save to; must not be empty."))
		((bool,truncate,false,,"Whether to delete current file contents, if any, when opening (false by default)"))
		((bool,addIterNum,false,,"Adds an iteration number to the file name, when the file was created. Useful for creating new files at each call (false by default)"))
	);
};
REGISTER_SERIALIZABLE(Recorder);
