#include"OpenGLManager.hpp"

CREATE_LOGGER(OpenGLManager);

OpenGLManager* OpenGLManager::self=NULL;

OpenGLManager::OpenGLManager(QObject* parent): QObject(parent){
	if(self) throw runtime_error("OpenGLManager instance already exists, uses OpenGLManager::self to retrieve it.");
	self=this;
	renderer=shared_ptr<OpenGLRenderer>(new OpenGLRenderer);
	renderer->init();
	connect(this,SIGNAL(createView()),this,SLOT(createViewSlot()));
	connect(this,SIGNAL(resizeView(int,int,int)),this,SLOT(resizeViewSlot(int,int,int)));
	connect(this,SIGNAL(closeView(int)),this,SLOT(closeViewSlot(int)));
	connect(this,SIGNAL(startTimerSignal()),this,SLOT(startTimerSlot()),Qt::QueuedConnection);
}

void OpenGLManager::timerEvent(QTimerEvent* event){
	//cerr<<".";
	boost::mutex::scoped_lock lock(viewsMutex);
	// when sharing the 0th view widget, it should be enough to update the primary view only
	//if(views.size()>0) views[0]->updateGL();
	#if 1
		FOREACH(const shared_ptr<GLViewer>& view, views){ if(view) view->updateGL(); }
	#endif
}

void OpenGLManager::createViewSlot(){
	boost::mutex::scoped_lock lock(viewsMutex);
	if(views.size()==0){
		views.push_back(shared_ptr<GLViewer>(new GLViewer(0,renderer,/*shareWidget*/(QGLWidget*)0)));
	} else {
		throw runtime_error("Secondary views not supported");
		//views.push_back(shared_ptr<GLViewer>(new GLViewer(views.size(),renderer,views[0].get())));
	}	
}

void OpenGLManager::resizeViewSlot(int id, int wd, int ht){
	views[id]->resize(wd,ht);
}

void OpenGLManager::closeViewSlot(int id){
	boost::mutex::scoped_lock lock(viewsMutex);
	for(size_t i=views.size()-1; (!views[i]); i--){ views.resize(i); } // delete empty views from the end
	if(id<0){ // close the last one existing
		assert(*views.rbegin()); // this should have been sanitized by the loop above
		views.resize(views.size()-1); // releases the pointer as well
	}
	if(id==0){
		LOG_DEBUG("Closing primary view.");
		if(views.size()==1) views.clear();
		else{ LOG_INFO("Cannot close primary view, secondary views still exist."); }
	}
}
void OpenGLManager::centerAllViews(){
	boost::mutex::scoped_lock lock(viewsMutex);
	FOREACH(const shared_ptr<GLViewer>& g, views){ if(!g) continue; g->centerScene(); }
}
void OpenGLManager::startTimerSlot(){
	startTimer(50);
}

int OpenGLManager::waitForNewView(float timeout,bool center){
	size_t origViewCount=views.size();
	emitCreateView();
	float t=0;
	while(views.size()!=origViewCount+1){
		usleep(50000); t+=.05;
		// wait at most 5 secs
		if(t>=timeout) {
			LOG_ERROR("Timeout waiting for the new view to open, giving up."); return -1;
		}
	}
	if(center)(*views.rbegin())->centerScene();
	return (*views.rbegin())->viewId; 
}
