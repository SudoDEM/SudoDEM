/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <boost/scoped_ptr.hpp>
#include<stdexcept>
#include<string>

/*! \brief Abstract interface for all Indexable class.
	An indexable class is a class that will be managed by a MultiMethodManager.
	The index the function getClassIndex() returns, corresponds to the index in the matrix where the class will be handled.
*/

#define _THROW_NOT_OVERRIDDEN  throw std::logic_error(std::string("Derived class did not override ")+__PRETTY_FUNCTION__+", use REGISTER_INDEX_COUNTER and REGISTER_CLASS_INDEX.")

class Indexable
{
	protected :
		void createIndex() {
			int& index = getClassIndex();
			if(index == -1)				// assign new index
			{
				index = getMaxCurrentlyUsedClassIndex()+1;
				// so that other dispatchers will not fall in conflict with this index
				incrementMaxCurrentlyUsedClassIndex();
			}
		}

	public :
		Indexable () {};
		virtual ~Indexable () {};

		/// Returns the id of the current class. This id is set by a multimethod manager
		virtual int& getClassIndex()                             { _THROW_NOT_OVERRIDDEN;};
		virtual const int& getClassIndex() const                 { _THROW_NOT_OVERRIDDEN;};
		virtual int& getBaseClassIndex(int )                     { _THROW_NOT_OVERRIDDEN;};
		virtual const int& getBaseClassIndex(int ) const         { _THROW_NOT_OVERRIDDEN;};
		virtual const int& getMaxCurrentlyUsedClassIndex() const { _THROW_NOT_OVERRIDDEN;};
		virtual void incrementMaxCurrentlyUsedClassIndex()       { _THROW_NOT_OVERRIDDEN;};

};

#undef _THROW_NOT_OVERRIDDEN

// this macro is used by classes that are a dimension in multimethod matrix

#define REGISTER_CLASS_INDEX(SomeClass,BaseClass)                                      \
	public: static int& getClassIndexStatic() { static int index = -1; return index; } \
	public: virtual int& getClassIndex()       { return getClassIndexStatic(); }        \
	public: virtual const int& getClassIndex() const { return getClassIndexStatic(); }  \
	public: virtual int& getBaseClassIndex(int depth) {              \
		static boost::scoped_ptr<BaseClass> baseClass(new BaseClass); \
		if(depth == 1) return baseClass->getClassIndex();             \
		else           return baseClass->getBaseClassIndex(--depth);  \
	}                                                                \
	public: virtual const int& getBaseClassIndex(int depth) const {  \
		static boost::scoped_ptr<BaseClass> baseClass(new BaseClass); \
		if(depth == 1) return baseClass->getClassIndex();             \
		else           return baseClass->getBaseClassIndex(--depth);  \
	}

// this macro is used by base class for classes that are a dimension in multimethod matrix
// to keep track of maximum number of classes of their kin. Multimethod matrix can't
// count this number (ie. as a size of the matrix), as there are many multimethod matrices

#define REGISTER_INDEX_COUNTER(SomeClass) \
	private: static int& getClassIndexStatic()       { static int index = -1; return index; }\
	public: virtual int& getClassIndex()             { return getClassIndexStatic(); }       \
	public: virtual const int& getClassIndex() const { return getClassIndexStatic(); }       \
	public: virtual int& getBaseClassIndex(int)             { throw std::logic_error("One of the following errors was detected:\n(1) Class " #SomeClass " called createIndex() in its ctor (but it shouldn't, being a top-level indexable; only use REGISTER_INDEX_COUNTER, but not createIndex()).\n(2) Some DerivedClass deriving from " #SomeClass " forgot to use REGISTER_CLASS_INDEX(DerivedClass," #SomeClass ").\nPlease fix that and come back again." ); } \
	public: virtual const int& getBaseClassIndex(int) const { throw std::logic_error("One of the following errors was detected:\n(1) Class " #SomeClass " called createIndex() in its ctor (but it shouldn't, being a top-level indexable; only use REGISTER_INDEX_COUNTER, but not createIndex()).\n(2) Some DerivedClass deriving from " #SomeClass " forgot to use REGISTER_CLASS_INDEX(DerivedClass," #SomeClass ").\nPlease fix that and come back again." ); } \
	private: static int& getMaxCurrentlyUsedIndexStatic() { static int maxCurrentlyUsedIndex = -1; return maxCurrentlyUsedIndex; } \
	public: virtual const int& getMaxCurrentlyUsedClassIndex() const {  \
		SomeClass * Indexable##SomeClass = 0;                            \
		Indexable##SomeClass = dynamic_cast<SomeClass*>(const_cast<SomeClass*>(this)); \
		if (Indexable##SomeClass) {                                      \
			assert(Indexable##SomeClass);                                  \
		}                                                                \
		return getMaxCurrentlyUsedIndexStatic();                         \
	}                                                                   \
	public: virtual void incrementMaxCurrentlyUsedClassIndex() {        \
		SomeClass * Indexable##SomeClass = 0;                            \
		Indexable##SomeClass = dynamic_cast<SomeClass*>(this);           \
		if (Indexable##SomeClass) {                                      \
			assert(Indexable##SomeClass);                                  \
		}                                                                \
		int& max = getMaxCurrentlyUsedIndexStatic();                     \
		max++;                                                           \
	}                                                                  \

// macro that should be passed in the 4th argument of SUDODEM_CLASS_BASE_ATTR_PY in the top-level indexable
#define SUDODEM_PY_TOPINDEXABLE(className) .add_property("dispIndex",&Indexable_getClassIndex<className>,"Return class index of this instance.").def("dispHierarchy",&Indexable_getClassIndices<className>,(boost::python::arg("names")=true),"Return list of dispatch classes (from down upwards), starting with the class instance itself, top-level indexable at last. If names is true (default), return class names rather than numerical indices.")


