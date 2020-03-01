#include<sudodem/lib/pyutil/gil.hpp>
void pyRunString(const std::string& cmd){
	gilLock lock; PyRun_SimpleString(cmd.c_str());
};

