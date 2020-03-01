// 2009 © Václav Šmilauer <eudoxos@arcig.cz
#pragma once
#include"numpy_boost.hpp"

// helper macro do assign Vector3r and Matrix3r values to subarrays
#define VECTOR3R_TO_NUMPY(vec,arr) arr[0]=vec[0]; arr[1]=vec[1]; arr[2]=vec[2]
#define MATRIX3R_TO_NUMPY(mat,arr) arr[0]=mat(0,0);arr[1]=mat(0,1);arr[2]=mat(0,2);arr[3]=mat(1,0);arr[4]=mat(1,1);arr[5]=mat(1,2);arr[6]=mat(2,0);arr[7]=mat(2,1);arr[8]=mat(2,2);
