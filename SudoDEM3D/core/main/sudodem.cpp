#include <Python.h>
#include <iostream>

#include <stdexcept>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <boost/process/search_path.hpp>

//#include <omp.h>

#include"sudodemcfg.h"

using namespace boost;
using namespace std;

void sudodem_print_help();
void sudodem_print_version();

int main( int argc, char **argv )
{
  cout<<"Welcome to SudoDEM!"<<endl;
  //setPathEnv(argc, argv);
  //set environment Variables
  filesystem::path exePath;
  if (argc > 0)
  {
    std::string exeName = argv[0];
    if (exeName.size() > 0)
    {
      std::string origPathValue = "";
      const char *getenvVal = getenv("PATH");
      if (getenvVal != NULL)
      {
        origPathValue = getenvVal;
        //cout<<"original PATH="<<origPathValue<<endl;
      }
      std::string origLDPathValue = "";
      const char *getenvLDVal = getenv("LD_LIBRARY_PATH");
      if (getenvLDVal != NULL)
      {
        origLDPathValue = getenvLDVal;
        //cout<<"original PATH="<<origPathValue<<endl;
      }
      //std::string fullpath = dll::program_location().parent_path().string();
      std::vector<boost::filesystem::path> searchpath;
      searchpath.push_back(filesystem::initial_path<filesystem::path>());
      //std::string currentpath = filesystem::initial_path<filesystem::path>().string();//current path at the terminal
      exePath = boost::process::search_path(filesystem::path(exeName), searchpath);//first search the current path
      //cout<<"exepath"<<exePath<<endl;
      if(exePath.empty()){
        //cout<<"searching the PATH."<<endl;
        exePath = boost::process::search_path(filesystem::path(exeName));
      }
      //filesystem::path exePath = filesystem::system_complete(filesystem::path(exeName));
      std::string prefix= (exePath.branch_path().string())+"/../'";
      std::string newPathValue = origPathValue + ":" + (exePath.branch_path().string());
      std::string libPathValue = (exePath.branch_path().string()) + "/../lib/sudodem/";
      libPathValue = libPathValue + ":" + (exePath.branch_path().string()) + "/../lib/sudodem/py";
      libPathValue = libPathValue + ":" + (exePath.branch_path().string()) + "/../lib/sudodem/py/sudodem";
      libPathValue = libPathValue + ":" + (exePath.branch_path().string()) + "/../lib/sudodem/py/sudodem/qt";
      libPathValue = libPathValue + ":" + (exePath.branch_path().string()) + "/../lib/3rdlibs/py";
      std::string libLDPathValue = (exePath.branch_path().string()) + "/../lib/3rdlibs/";
      //cout<<"the path = "<<exePath.string()<<endl;
      setenv("PATH", newPathValue.c_str(), 1);
      //cout<<"Prefix: "<<prefix<<endl;
      setenv("PYTHONPATH",libPathValue.c_str(),1);
      setenv("LD_LIBRARY_PATH",libLDPathValue.c_str(),1);
      setenv("SUDODEM_PREFIX",prefix.c_str(),1);
      //cout<<"updated PATH="<<newPathValue<<endl;
    }
  }
  //options
  setenv("OMP_NUM_THREADS", "1", 1);
  std :: vector< const char * >modulesArgs;
  bool isGUI=true;
  if(argc != 1){
    for(int i=1;i < argc; i++){
      if ( strcmp(argv[i], "-v") == 0 ){//version
        sudodem_print_version();
        if(argc == 2){
          //just for enquiring version
          exit(1);
        }
      }else if (strcmp(argv[i], "-h") == 0 ){//help
        sudodem_print_help();
        if(argc == 2){
          //just for enquiring version
          exit(1);
        }
      }
      else if( strncmp(argv[i], "-j",2) == 0 ){//this first two characters are '-j'
        int len = strlen(argv[i]);
        if(len > 2){//i.e., the user may type '-j4' instead of a standard way '-j 4'
          char cores[2];
          cores[0]='1';cores[1]=' ';
          switch (len) {
            case 3:
              strncpy(cores,argv[i]+2,1);//the core number should not exceed 99 here by default.
              break;
            case 4:
              strncpy(cores,argv[i]+2,2);//the core number should not exceed 99 here by default.
              break;
            default:
              cout<<"Warning: wrong number of threads was assigned to this sim. Please check the command input. Only one thread will be used for continuing running."<<endl;
          }
          //cout<<"cores="<<cores<<"len="<<len<<endl;
          setenv("OMP_NUM_THREADS", cores, 1);
        }else if ( i + 1 < argc ) {
          i++;
          //putenv("OMP_NUM_THREADS=2");
          setenv("OMP_NUM_THREADS", argv[i], 1);
        }else{
          cout<<"The number of thread has been set incorrectly!"<<endl;
          setenv("OMP_NUM_THREADS", "1", 1);
        }
      }
      else if(strcmp(argv[i], "-n") == 0){
        //no gui
        isGUI = false;
      }
      else if(strcmp(argv[i], "-cores") == 0){
        //todo
      }
      else{//
        modulesArgs.push_back(argv[i]);
      }
    }
  }else{
    sudodem_print_help();
  }
  //program features
  std::string feat = PRG_FEATURE;

  std::istringstream iss(feat);
  std::vector<std::string> prg_feats((std::istream_iterator<std::string>(iss)),
                                 std::istream_iterator<std::string>());

  wchar_t* program = Py_DecodeLocale(argv[0], NULL);
  Py_SetProgramName(program);  /* optional but recommended */


  //cout<<"omp_get_max_threads="<<omp_get_max_threads()<<endl;
  Py_Initialize();
  
  PyRun_SimpleString("import sys");

  std::string cmd = "sysArgv =[";
  for(int i=0;i<argc;i++){
    //cout<<argv[i]<<endl;
    cmd = cmd +"'"+argv[i]+"',";
  }
  cmd +="]";
  //cout<<"cmd="<<cmd<<endl;
  PyRun_SimpleString(cmd.c_str());
  
PyRun_SimpleString(
  //"opts = None\n"
  "exitAfter=False\n"
);
cmd = "args =[";
for(int i=0;i<modulesArgs.size();i++){
  //cout<<argv[i]<<endl;
  cmd = cmd +"'"+modulesArgs[i]+"',";
}
cmd +="]";
//cout<<"cmd="<<cmd<<endl;
PyRun_SimpleString(cmd.c_str());

//putenv("OMP_NUM_THREADS=2");
PyRun_SimpleString(
  //initialization and c++ plugins import
  "import sudodem\n"
  //other parts we will need soon
  //"import sudodem.config\n"
  "import sudodem.wrapper\n"
  "import sudodem.system\n"
  "import sudodem.runtime\n"
  "sys.argv=sudodem.runtime.argv=args\n"
  "sudodem.runtime.exitAfter=exitAfter\n"
  "from sudodem import utils\n"
  "from sudodem.utils import *\n"
  "from math import *\n"
  "gui=None\n"
  "sudodem.runtime.hasDisplay=False \n"//# this is the default initialized in the module, anyway
);
if(!isGUI){
  PyRun_SimpleString(
    "gui=None\n"
  );
}else{
  for(int i=0;i<prg_feats.size();i++){
    //cout<<"features"<<prg_feats[i]<<endl;
     if(strcmp(prg_feats[i].c_str(), "GUI") == 0){
       PyRun_SimpleString(
         "gui='qt4'\n"
       );
       break;
     }
  }
  PyRun_SimpleString(
  "if gui:\n"
      "\timport Xlib.display\n"
      // PyQt4's QApplication does exit(1) if it is unable to connect to the display
      // we however want to handle this gracefully, therefore
      // we test the connection with bare xlib first, which merely raises DisplayError
      "\ttry:\n"
          // contrary to display.Display, _BaseDisplay does not check for extensions and that avoids spurious message "Xlib.protocol.request.QueryExtension" (bug?)
          "\t\tXlib.display._BaseDisplay();\n"
          "\t\tsudodem.runtime.hasDisplay=True\n"
          //"\t\tprint(sudodem.runtime.hasDisplay)\n"
      "\texcept:\n"
          // usually Xlib.error.DisplayError, but there can be Xlib.error.XauthError etc as well
          // let's just pretend any exception means the display would not work
          "\t\tgui=None\n"
          "\t\tprint('Warning: no X rendering available')\n"
  );
}
PyRun_SimpleString(
"def userSession(qt4=False,qapp=None):\n"
    //# prepare nice namespace for users
    "\timport sudodem.runtime,sys\n"
    //import sys
    "\tglobal gui\n"
    "\tif len(sys.argv)>0:\n"
        "\t\targ0=sys.argv[0]\n"
        "\t\tif qt4: sudodem.qt.Controller();\n"
        //"\t\tif sum(bool(arg0.endswith(ext)) for ext in ('.xml','.xml.bz2','.xml.gz','.sudodem','.sudodem.gz','.sudodem.bz2','.bin','.bin.gz','.bin.bz2'))>0:\n"
        //    "\t\t\tif len(sys.argv)>1: raise RuntimeError('Extra arguments to saved simulation to run: '+' '.join(sys.argv[1:]))\n"
        //    "\t\t\tsys.stderr.write('Running simulation '+arg0)\n"
        "\t\tif arg0.endswith('.py'):\n"
            "\t\t\tdef runScript(script):\n"
                "\t\t\t\tsys.stderr.write('Running script '+arg0)\n"
                "\t\t\t\ttry: execfile(script,globals())\n"
                "\t\t\t\texcept SystemExit: raise\n"
                "\t\t\t\texcept: # all other exceptions\n"
                    "\t\t\t\t\timport traceback\n"
                    "\t\t\t\t\ttraceback.print_exc()\n"
                    "\t\t\t\t\tif sudodem.runtime.exitAfter: sys.exit(1)\n"
                "\t\t\t\tif sudodem.runtime.exitAfter: sys.exit(0)\n"
            "\t\t\trunScript(arg0)\n"
    "\tif sudodem.runtime.exitAfter: sys.exit(0)\n"
    //# common ipython configuration
    //"\tbanner='[[ Ctrl+L clears screen, Ctrl+U kills line. '+', '.join((['F12 controller','F11 3d view (use h-key for showing help)','F10 both','F9 generator'] if (qt4) else [])+['F8 plot'])+'. ]]'\n"
    //# ipython options, see e.g. http://www.cv.nrao.edu/~rreid/casa/tips/ipy_user_conf.py
    //#execfile=[prefix+'/lib/sudodem'+suffix+'/py/sudodem/ipython.py'],
    "\tipconfig=dict(prompt_in1='SudoDEM [\\#]: ',prompt_in2='     .\\D.: ',prompt_out=' ->  [\\#]: ',\
    separate_in='',separate_out='',separate_out2='',\
    readline_parse_and_bind=['tab: complete',] +\
    (['\"\\e[24~\": \"\\C-Usudodem.qt.Controller();\\C-M\"',\
    '\"\\e[23~\": \"\\C-Usudodem.qt.View();\\C-M\"',\
    '\"\\e[21~\": \"\\C-Usudodem.qt.Controller(), sudodem.qt.View();\\C-M\"',\
    '\"\\e[20~\": \"\\C-Usudodem.qt.Generator();\\C-M\"'] if (qt4) else []) +\
    ['\"\\e[19~\": \"\\C-Uimport sudodem.plot; sudodem.plot.plot();\\C-M\"', \
     '\"\\e[A\": history-search-backward', '\"\\e[B\": history-search-forward'])\n"
     //#F12,F11,F10,F9,F8
  //"\tipconfig=dict(prompt_in1='SudoDEM [\\#]: ',prompt_in2='     .\\D.: ',prompt_out=' ->  [\\#]: ')\n"
  //"\tprint ipconfig\n"
  //# show python console IPYTHON3.0 embeded
	"\tfrom IPython.terminal.embed import InteractiveShellEmbed\n"
	//"\tfrom IPython.config.loader import Config\n"
    "\tfrom traitlets.config.loader import Config\n"
	"\tcfg = Config()\n"
	"\tprompt_config = cfg.PromptManager\n"
	"\tprompt_config.in_template = ipconfig['prompt_in1']\n"
	"\tprompt_config.in2_template = ipconfig['prompt_in2']\n"
	"\tprompt_config.out_template = ipconfig['prompt_out']\n"
	"\timport readline\n"
	"\tfor k in ipconfig['readline_parse_and_bind']: readline.parse_and_bind(k)\n"
	"\tInteractiveShellEmbed.config=cfg\n"
	"\tInteractiveShellEmbed.banner1=''\n"
	"\tipshell=InteractiveShellEmbed()\n"
	"\tipshell()\n"
);
cout<<"Tips: Ctrl+L clears screen, Ctrl+U kills line. F12 controller, use h-key for 3d view help,F9 generator, F8 plot"<<endl;
PyRun_SimpleString(
"if gui==None:\n"
        "\tsys.exit(userSession())\n"
    "elif gui=='qt4':\n"
        // we already tested that DISPLAY is available and can be opened
        // otherwise Qt4 might crash at this point
        "\timport PyQt5\n"
        "\tfrom PyQt5 import QtGui,QtCore,QtWidgets\n"
        "\timport sudodem.qt\n"
        "\tqapp=QtWidgets.QApplication(sys.argv)\n"
        "\tsys.exit(userSession(qt4=True,qapp=qapp))\n"
      );

    Py_Finalize();
    return 0;
}


void sudodem_print_help(){
    cout<<"\nUsage:\n"<<endl;
    cout<<"sudodemexe [options] user-script.py"<<endl;
    cout<<"\nOptions:\n"<<endl;
    cout<<"  -v  prints SudoDEM version"<<endl;
    cout<<"  -j  (int) Number of OpenMP threads to run (default 1). Equivalent to setting OMP_NUM_THREADS environment variable."<<endl;
    cout<<"  -n  Run without graphical interface (equivalent to unsetting the DISPLAY environment variable)"<<endl;
}

void sudodem_print_version(){
   //printf("\n%s (%s, %s)\nof %s on %s\n", PRG_VERSION, HOST_TYPE, MODULE_LIST, __DATE__, HOST_NAME);
   printf("\n%s\n",PRG_VERSION);
}
