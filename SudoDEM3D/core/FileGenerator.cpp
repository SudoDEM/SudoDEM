/*************************************************************************
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<sudodem/core/Omega.hpp>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<sudodem/lib/pyutil/gil.hpp>
#include<sudodem/lib/serialization/ObjectIO.hpp>

#include"FileGenerator.hpp"

CREATE_LOGGER(FileGenerator);


bool FileGenerator::generate(std::string& msg){ throw invalid_argument("Calling abstract FileGenerator::generate() does not make sense."); }


bool FileGenerator::generateAndSave(const string& outputFileName, string& message)
{
	bool status;
	message="";
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	try {
		status=generate(message); // will modify message
	}
	catch(std::exception& e){
		LOG_FATAL("Unhandled exception: "<<typeid(e).name()<<" : "<<e.what());
		//abort(); // use abort, since we may want to inspect core
		message = message + "Unhandled exception: " + typeid(e).name() + " : " + e.what();
		return false;
	}
	// generation wasn't successful
	if(status==false) return false;

	else {
		boost::posix_time::ptime now2 = boost::posix_time::second_clock::local_time();
		boost::posix_time::time_duration generationTime = now2 - now; // generation time, without save time
		try
		{
			sudodem::ObjectIO::save(outputFileName,"scene",scene);
		}
		catch(const std::runtime_error& e)
		{
			message+=std::string("File "+outputFileName+" cannot be saved: "+e.what());
			return false;
		}
		boost::posix_time::ptime now3 = boost::posix_time::second_clock::local_time();
		boost::posix_time::time_duration saveTime = now3 - now2; // save time
		message=std::string("File "+outputFileName+" generated successfully."
				+ "\ngeneration time: " + boost::posix_time::to_simple_string(generationTime)
				+ "\nsave time: "       + boost::posix_time::to_simple_string(saveTime)
				+"\n\n")+message;
		return true;
	}
}

void FileGenerator::pyGenerate(const string& out){
	string message;
	bool ret=generateAndSave(out,message);
	LOG_INFO((ret?"SUCCESS:\n":"FAILURE:\n")<<message);
	if(ret==false) throw runtime_error(getClassName()+" reported error: "+message);
}
void FileGenerator::pyLoad(){
	string xml(Omega::instance().tmpFilename()+".xml.bz2");
	// LOG_DEBUG("Using temp file "<<xml);
	pyGenerate(xml);
	//this is ugly hack, yes...
	pyRunString("sudodem.wrapper.Omega().load('"+xml+"')");
}
