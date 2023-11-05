#pragma once
//#include

#ifndef Q_MOC_RUN 
#include"GLViewer.hpp"
#endif 

#include<QObject>
#include<boost/thread/mutex.hpp>

/*
Singleton class managing OpenGL views,
a renderer instance and timer to refresh the display.
*/
class OpenGLManager: public QObject{
	Q_OBJECT
	DECLARE_LOGGER;
	public:
		static OpenGLManager* self;
		OpenGLManager(QObject *parent=0);
		// manipulation must lock viewsMutex!
		std::vector<shared_ptr<GLViewer> > views;
		shared_ptr<OpenGLRenderer> renderer;
		// signals are protected, emitting them is therefore wrapped with such funcs
		void emitResizeView(int id, int wd, int ht){ emit resizeView(id,wd,ht); }
		void emitCreateView(){ emit createView(); }
		void emitStartTimer(){ emit startTimerSignal(); }
		void emitCloseView(int id){ emit closeView(id); }
		// create a new view and wait for it to become available; return the view number
		// if timout (in seconds) elapses without the view to come up, reports error and returns -1
		int waitForNewView(float timeout=5., bool center=true);
	signals:
		void createView();
		void resizeView(int id, int wd, int ht);
		void closeView(int id);
		// this is used to start timer from the main thread via postEvent (ugly)
		void startTimerSignal();
	public slots:
		virtual void createViewSlot();
		virtual void resizeViewSlot(int id, int wd, int ht);
		virtual void closeViewSlot(int id=-1);
		virtual void timerEvent(QTimerEvent* event);
		virtual void startTimerSlot();
		void centerAllViews();
	private:
		boost::mutex viewsMutex;
};
