/**
 *       @file  utlCoreMKL.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlCoreMKL_h
#define __utlCoreMKL_h

#include <cmath>

#ifdef UTL_USE_MKL
#include "utlMKL.h"
#endif

/** @addtogroup utlMath
@{ */
namespace utl
{

// interfaces to a few functions from the intel Vector Mathematical Library
/// interface to v*Add
template <typename T> void vAdd( int n,  T* vecIn,  T* vecIn2, T* vecOut);
/// interface to v*Sub
template <typename T> void vSub( int n,  T* vecIn,  T* vecIn2, T* vecOut);
/// interface to v*Mul
template <typename T> void vMul( int n,  T* vecIn,  T* vecIn2, T* vecOut);
/// interface to v*Div
template <typename T> void vDiv( int n,  T* vecIn,  T* vecIn2, T* vecOut);
/// interface to v*Sqr
template <typename T> void vSqr( int n,  T* vecIn, T* vecOut);
/// interface to v*Abs
template <typename T> void vAbs( int n,  T* vecIn, T* vecOut);
/// interface to v*Exp
template <typename T> void vExp( int n,  T* vecIn, T* vecOut);
/// interface to v*Inv
template <typename T> void vInv( int n,  T* vecIn, T* vecOut);
/// interface to v*Sqrt
template <typename T> void vSqrt( int n,  T* vecIn, T* vecOut);
/// interface to v*InvSqrt
template <typename T> void vInvSqrt( int n,  T* vecIn, T* vecOut);
/// interface to v*Cos
template <typename T> void vCos( int n,  T* vecIn, T* vecOut);
/// interface to v*Sin
template <typename T> void vSin( int n,  T* vecIn, T* vecOut);


#ifdef UTL_USE_MKL 

/// Implemenation of the interface for vdSqr
template <> inline void 
vSqr<double>( int n,  double* vecIn, double* vecOut) 
{
   vdSqr(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsSqr
template <> inline void 
vSqr<float>( int n,  float* vecIn, float* vecOut) 
{
   vsSqr(n,vecIn,vecOut);
}
template <> inline void 
vSqrt<double>( int n,  double* vecIn, double* vecOut) 
{
   vdSqrt(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsSqr
template <> inline void 
vSqrt<float>( int n,  float* vecIn, float* vecOut) 
{
   vsSqrt(n,vecIn,vecOut);
}
template <> inline void 
vInvSqrt<double>( int n,  double* vecIn, double* vecOut) 
{
   vdInvSqrt(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsSqr
template <> inline void 
vInvSqrt<float>( int n,  float* vecIn, float* vecOut) 
{
   vsInvSqrt(n,vecIn,vecOut);
}

/// Implemenation of the interface for vdSub
template <> inline void 
vSub<double>( int n,  double* vecIn, double* vecIn2, double* vecOut) 
{
   vdSub(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vsSub
template <> inline void 
vSub<float>( int n,  float* vecIn, float* vecIn2, float* vecOut) 
{
   vsSub(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vdDiv
template <> inline void 
vDiv<double>( int n,  double* vecIn, double* vecIn2, double* vecOut) 
{
   vdDiv(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vsDiv
template <> inline void 
vDiv<float>( int n,  float* vecIn, float* vecIn2, float* vecOut) 
{
   vsDiv(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vdExp
template <> inline void 
vExp<double>( int n,  double* vecIn, double* vecOut) 
{
   vdExp(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsExp
template <> inline void 
vExp<float>( int n,  float* vecIn, float* vecOut) 
{
   vsExp(n,vecIn,vecOut);
}
/// Implemenation of the interface for vdCos
template <> inline void 
vCos<double>( int n,  double* vecIn, double* vecOut) 
{
   vdCos(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsCos
template <> inline void 
vCos<float>( int n,  float* vecIn, float* vecOut) 
{
   vsCos(n,vecIn,vecOut);
}
/// Implemenation of the interface for vdSin
template <> inline void 
vSin<double>( int n,  double* vecIn, double* vecOut) 
{
   vdSin(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsSin
template <> inline void 
vSin<float>( int n,  float* vecIn, float* vecOut) 
{
   vsSin(n,vecIn,vecOut);
}
/// Implemenation of the interface for vdInv
template <> inline void 
vInv<double>( int n,  double* vecIn, double* vecOut) 
{
   vdInv(n,vecIn,vecOut);
}
/// Implemenation of the interface for vsInv
template <> inline void 
vInv<float>( int n,  float* vecIn, float* vecOut) 
{
   vsInv(n,vecIn,vecOut);
}
/// Implemenation of the interface for vdAdd
template <> inline void 
vAdd<double>( int n,  double* vecIn, double* vecIn2, double* vecOut) 
{
   vdAdd(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vsAdd
template <> inline void 
vAdd<float>( int n,  float* vecIn, float* vecIn2, float* vecOut) 
{
   vsAdd(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vdMul
template <> inline void 
vMul<double>( int n,  double* vecIn, double* vecIn2, double* vecOut) 
{
   vdMul(n,vecIn,vecIn2,vecOut);
}
/// Implemenation of the interface for vsMul
template <> inline void 
vMul<float>( int n,  float* vecIn, float* vecIn2, float* vecOut) 
{
   vsMul(n,vecIn,vecIn2,vecOut);
}

/// interface to vdAbs
template <> inline void 
vAbs( int n,  double* vecIn, double* vecOut) 
{
   vdAbs(n,vecIn,vecOut);
}
/// interface to vdAbs
template <> inline void 
vAbs( int n,  float* vecIn, float* vecOut) 
{
   vsAbs(n,vecIn,vecOut);
}

#else

template <typename T> inline void 
vSqr( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=vecIn[i]*vecIn[i];
}
template <typename T> inline void 
vSqrt( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=std::sqrt(vecIn[i]);
}
template <typename T> inline void 
vInvSqrt( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=T(1.0)/std::sqrt(vecIn[i]);
}

/// Slow implementation of vdSub and vsSub
template <typename T> inline void 
vSub( int n,  T* vecIn1, T* vecIn2, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=vecIn1[i]-vecIn2[i];
}
/// Slow implementation of vdInv and vsInv
template <typename T> inline void 
vInv( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=1.0/vecIn[i];
}
/// Slow implementation of vdExp and vsExp
template <typename T> inline void 
vExp( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=std::exp(vecIn[i]);
}
/// Slow implementation of vdAdd and vsAdd
template <typename T> inline void 
vAdd( int n,  T* vecIn1, T* vecIn2, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=vecIn1[i]+vecIn2[i];
}
/// Slow implementation of vdMul and vsMul
template <typename T> inline void 
vMul( int n,  T* vecIn1, T* vecIn2, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=vecIn1[i]*vecIn2[i];
}
/// Slow implementation of vdDiv and vsDiv
template <typename T> inline void 
vDiv( int n,  T* vecIn1, T* vecIn2, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=vecIn1[i]/vecIn2[i];
}
/// Slow implementation of vAbs
template <typename T> inline void 
vAbs( int n,  T* vecIn, T* vecOut) 
{
   for (int i = 0; i<n; ++i) 
     vecOut[i]=std::fabs(vecIn[i]);
}

#endif

}

/** @}  */

#endif 
