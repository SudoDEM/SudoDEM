#pragma once
#include<boost/python/raw_function.hpp>
// many thanks to http://markmail.org/message/s4ksg6nfspw2wxwd
namespace boost { namespace python { namespace detail {
	template <class F> struct raw_constructor_dispatcher{
		raw_constructor_dispatcher(F f): f(make_constructor(f)) {}
		PyObject* operator()(PyObject* args, PyObject* keywords)
		{
			 borrowed_reference_t* ra = borrowed_reference(args); object a(ra);
			 return incref(object(f(object(a[0]),object(a.slice(1,len(a))),keywords ? dict(borrowed_reference(keywords)) : dict())).ptr() );
		}
		private: object f;
	};
	}
	template <class F> object raw_constructor(F f, std::size_t min_args = 0){
		return detail::make_raw_function(objects::py_function(detail::raw_constructor_dispatcher<F>(f),mpl::vector2<void, object>(),min_args+1,(std::numeric_limits<unsigned>::max)()));
	}
}} // namespace boost::python


