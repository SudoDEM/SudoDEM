#ifdef SUDODEM_OPENGL

#include"SnapshotEngine.hpp"

SUDODEM_PLUGIN((SnapshotEngine));
CREATE_LOGGER(SnapshotEngine);

void SnapshotEngine::action(){
	if(!OpenGLManager::self) throw logic_error("No OpenGLManager instance?!");
	if(OpenGLManager::self->views.size()==0){
		int viewNo=OpenGLManager::self->waitForNewView(deadTimeout);
		if(viewNo<0){
			if(!ignoreErrors) throw runtime_error("SnapshotEngine: Timeout waiting for new 3d view.");
			else {
				LOG_WARN("Making myself Engine::dead, as I can not live without a 3d view (timeout)."); dead=true; return;
			}
		}
	}
	const shared_ptr<GLViewer>& glv=OpenGLManager::self->views[0];
	ostringstream fss; fss<<fileBase<<setw(5)<<setfill('0')<<counter++<<"."<<boost::algorithm::to_lower_copy(format);
	LOG_DEBUG("GL view → "<<fss.str())
	glv->setSnapshotFormat(QString(format.c_str()));
	glv->nextFrameSnapshotFilename=fss.str();
	// wait for the renderer to save the frame (will happen at next postDraw)
	timespec t1,t2; t1.tv_sec=0; t1.tv_nsec=10000000; /* 10 ms */
	long waiting=0;
	while(!glv->nextFrameSnapshotFilename.empty()){
		nanosleep(&t1,&t2); waiting++;
		if(((waiting) % 1000)==0) LOG_WARN("Already waiting "<<waiting/100<<"s for snapshot to be saved. Something went wrong?");
		if(waiting/100.>deadTimeout){
			if(ignoreErrors){ LOG_WARN("Timeout waiting for snapshot to be saved, making byself Engine::dead"); dead=true; return; }
			else throw runtime_error("SnapshotEngine: Timeout waiting for snapshot to be saved.");
		}
	}
	snapshots.push_back(fss.str());
	usleep((long)(msecSleep*1000));
	//if(!plot.empty()){ pyRunString("import sudodem.plot; sudodem.plot.addImgData("+plot+"='"+fss.str()+"')"); }
}


#endif /* SUDODEM_OPENGL */
