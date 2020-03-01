#include<sudodem/core/Omega.hpp>
#include<sudodem/lib/base/Logging.hpp>

#include<signal.h>

#ifdef SUDODEM_DEBUG
	void crashHandler(int sig){
	switch(sig){
		case SIGABRT:
		case SIGSEGV:
			signal(SIGSEGV,SIG_DFL); signal(SIGABRT,SIG_DFL); // prevent loops - default handlers
			cerr<<"SIGSEGV/SIGABRT handler called; gdb batch file is `"<<Omega::instance().gdbCrashBatch<<"'"<<endl;
			std::system((string("gdb -x ")+Omega::instance().gdbCrashBatch).c_str());
			raise(sig); // reemit signal after exiting gdb
			break;
		}
	}
#endif

/* Initialize sudodem, loading given plugins */
void sudodemInitialize(boost::python::list& pp, const std::string& confDir){

	PyEval_InitThreads();

	Omega& O(Omega::instance());
	O.init();
	O.origArgv=NULL; O.origArgc=0; // not needed, anyway
	O.confDir=confDir;
	O.initTemps();
	#ifdef SUDODEM_DEBUG
		ofstream gdbBatch;
		O.gdbCrashBatch=O.tmpFilename();
		gdbBatch.open(O.gdbCrashBatch.c_str()); gdbBatch<<"attach "<<boost::lexical_cast<string>(getpid())<<"\nset pagination off\nthread info\nthread apply all backtrace\ndetach\nquit\n"; gdbBatch.close();
		signal(SIGABRT,crashHandler);
		signal(SIGSEGV,crashHandler);
	#endif
	vector<string> ppp; for(int i=0; i<boost::python::len(pp); i++) ppp.push_back(boost::python::extract<string>(pp[i]));
	Omega::instance().loadPlugins(ppp);
}
void sudodemFinalize(){ Omega::instance().cleanupTemps(); }

BOOST_PYTHON_MODULE(boot){
  boost::python::scope().attr("initialize")=&sudodemInitialize;
  boost::python::scope().attr("finalize")=&sudodemFinalize; //,"Finalize sudodem (only to be used internally).")
}
