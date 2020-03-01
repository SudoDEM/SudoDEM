/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ClassFactory.hpp"

#include<boost/algorithm/string/regex.hpp>
#include<sudodem/lib/base/Logging.hpp>

CREATE_LOGGER(ClassFactory);
SINGLETON_SELF(ClassFactory);

class Factorable;

bool ClassFactory::registerFactorable( std::string name 			   , CreateFactorableFnPtr create,
					 CreateSharedFactorableFnPtr createShared, CreatePureCustomFnPtr createPureCustom)
{

	bool tmp = map.insert( FactorableCreatorsMap::value_type( name , FactorableCreators(create,createShared, createPureCustom) )).second;

	#if 0
		if (tmp)
			std::cout << "registering factorable: " << name << " OK\n";
		else
			std::cout << "registering factorable: " << name << " FAILED\n";
	#endif

	return tmp;
}

shared_ptr<Factorable> ClassFactory::createShared( std::string name )
{
#if 0
	cerr<<"Creating shared lib: "<<name<<"\n";
	cerr<<"Available libs: ";
		for(FactorableCreatorsMap::iterator i=map.begin(); i!=map.end(); i++){
			cerr<<i->first<<" ";
		}
	cerr<<"\n";
#endif

	FactorableCreatorsMap::const_iterator i = map.find( name );
	if( i == map.end() )
	{
		dlm.load(name);
		if (dlm.isLoaded(name))
		{
			if( map.find( name ) == map.end() )
			{
				// Well, exception are also a way to return value, right?
				// This throws at startup for every .so that doesn't contain class named the same as the library.
				// I.e. almost everything in sudodem-libs and some others installed locally...
				// Don't log that, it would confuse users.
				//LOG_FATAL("Can't create class "<<name<<" - was never registered.");
				throw std::runtime_error(("Class "+name+" not registered in the ClassFactory.").c_str());
			}
			return createShared(name);
		}
		throw std::runtime_error(("Class "+name+" could not be factored in the ClassFactory.").c_str());
	}
	return ( i -> second.createShared ) ();
}

Factorable* ClassFactory::createPure( std::string name )
{
	FactorableCreatorsMap::const_iterator i = map.find( name );
	if( i == map.end() )
	{
		//cerr << "------------ going to load something" << endl;
		dlm.load(name);
		if (dlm.isLoaded(name))
		{
			if( map.find( name ) == map.end() )
			{
				throw std::runtime_error(("Class "+name+" not registered in the ClassFactory.").c_str());
			}
			return createPure(name);
		}
		throw std::runtime_error(("Class "+name+" could not be factored in the ClassFactory.").c_str());
	}
	return ( i -> second.create ) ();
}

void * ClassFactory::createPureCustom( std::string name )
{
	FactorableCreatorsMap::const_iterator i = map.find( name );
	if( i == map.end() ) throw std::runtime_error(("Class "+name+" could not be factored in the ClassFactory.").c_str());
	return ( i -> second.createPureCustom ) ();
}

bool ClassFactory::load(const string& name)
{
        return dlm.load(name);
}

string ClassFactory::lastError()
{
        return dlm.lastError();
}


void ClassFactory::registerPluginClasses(const char* fileAndClasses[]){
	assert(fileAndClasses[0]!=NULL); // must be file name
	// only filename given, no classes names explicitly
	if(fileAndClasses[1]==NULL){
		/* strip leading path (if any; using / as path separator) and strip one suffix (if any) to get the contained class name */
		string heldClass=boost::algorithm::replace_regex_copy(string(fileAndClasses[0]),boost::regex("^(.*/)?(.*?)(\\.[^.]*)?$"),string("\\2"));
		#ifdef SUDODEM_DEBUG
			if(getenv("SUDODEM_DEBUG")) cerr<<__FILE__<<":"<<__LINE__<<": Plugin "<<fileAndClasses[0]<<", class "<<heldClass<<" (deduced)"<<endl;
		#endif
		pluginClasses.push_back(heldClass); // last item with everything up to last / take off and .suffix strip
	}
	else {
		for(int i=1; fileAndClasses[i]!=NULL; i++){
			#ifdef SUDODEM_DEBUG
				if(getenv("SUDODEM_DEBUG")) cerr<<__FILE__<<":"<<__LINE__<<": Plugin "<<fileAndClasses[0]<<", class "<<fileAndClasses[i]<<endl;
			#endif
			pluginClasses.push_back(fileAndClasses[i]);
		}
	}
}

