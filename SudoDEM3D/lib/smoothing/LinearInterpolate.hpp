// 2008 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<vector>
#include<iostream>
#include<cassert>

/* Linear interpolation function suitable only for sequential interpolation.
 *
 * Template parameter T must support: addition, subtraction, scalar multiplication.
 *
 * @param t "time" at which the interpolated variable should be evaluated.
 * @param tt "time" points at which values are given; must be increasing.
 * @param values values at "time" points specified by tt
 * @param pos holds lookup state
 *
 * @return value at "time" t; out of range: t<t0 → value(t0), t>t_last → value(t_last)
 */
template<typename T, typename timeT>
T linearInterpolate(const Real t, const std::vector<timeT>& tt, const std::vector<T>& values, size_t& pos){
	assert(tt.size()==values.size());
	if(t<=tt[0]){ pos=0; return values[0];}
	if(t>=*tt.rbegin()){ pos=tt.size()-2; return *values.rbegin();}
	pos=std::min(pos,tt.size()-2);
	while((tt[pos]>t) || (tt[pos+1]<t)){
		assert(tt[pos]<tt[pos+1]);
		if(tt[pos]>t) pos--;
		else pos++;
	}
	const Real& t0=tt[pos], t1=tt[pos+1]; const T& v0=values[pos], v1=values[pos+1];
	return v0+(v1-v0)*((t-t0)/(t1-t0));
}

#if 0
	// test program
	int main(void){
		Real t,v;
		std::vector<Real> tt,vv;
		while(std::cin){
			std::cin>>t>>v;
			tt.push_back(t); vv.push_back(v);	
		}
		size_t pos;
		for(Real t=0; t<10; t+=0.1){
			std::cout<<t<<" "<<linearInterpolate<Real,Real>(t,tt,vv,pos)<<std::endl;
		}
	}
#endif 
