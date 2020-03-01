// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include <boost/thread/mutex.hpp>

#define FRIEND_SINGLETON(Class) friend class Singleton<Class>;					
// use to instantiate the self static member.
#define SINGLETON_SELF(Class) template<> Class* Singleton<Class>::self=NULL;
namespace { boost::mutex singleton_constructor_mutex; }
template <class T>
class Singleton{
	protected:
		static T* self; // must not be method-local static variable, since it gets created in different translation units multiple times.
	public:
		static T& instance(){
			if(!self) {
				boost::mutex::scoped_lock lock(singleton_constructor_mutex);
				if(!self) self=new T;
			}
			return *self;
		}
};


