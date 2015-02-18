/*!
 * \file mexutils.h
 *
 * \brief Contains miscellaneous functions for mex files.  
 *       utl functions for mex code. 
 *       Some codes are from spams by Julien Mairal julien.mairal@inria.fr
 *
 * */

#ifndef MEXUTILS_H
#define MEXUTILS_H

#include <mex.h>
#include <mat.h>
#include "utlCoreMacro.h"

#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#ifndef MATLAB_MEX_FILE
#define MATLAB_MEX_FILE
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <complex>

// #include <utils.h>
// #include <misc.h>


//#ifndef EM64T
//#define mwSize int
//#endif

namespace utl
{


/// Check the type of an array
template <typename T>
bool mexCheckType(const mxArray* array);

/// Check the type of an array (double)
template <> inline bool mexCheckType<double>(const mxArray* array) {
   return mxGetClassID(array) == mxDOUBLE_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (float)
template <> inline bool mexCheckType<float>(const mxArray* array) {
   return mxGetClassID(array) == mxSINGLE_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (int)
template <> inline bool mexCheckType<int>(const mxArray* array) {
   return mxGetClassID(array) == mxINT32_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (int)
template <> inline bool mexCheckType<bool>(const mxArray* array) {
   return mxGetClassID(array) == mxLOGICAL_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (complex double)
template <> inline bool mexCheckType<std::complex<double> >(const mxArray* array) {
   return mxGetClassID(array) == mxDOUBLE_CLASS && mxIsComplex(array);
};

/// Check the type of an array (complex float)
template <> inline bool mexCheckType<std::complex<float> >(const mxArray* array) {
   return mxGetClassID(array) == mxSINGLE_CLASS && mxIsComplex(array);
};


/// Check the size of a 2D-array
inline bool 
CheckSize(const mxArray* array, const int m, const int n) 
{
   const mwSize* dims=mxGetDimensions(array);
   int _m=static_cast<int>(dims[0]);
   int _n=static_cast<int>(dims[1]);
   return _n==n && _m==m;
};

/// Create a sparse copy of an array. Useful to deal with non-standard 
/// sparse matlab matrices
template <typename T>
void 
CreateCopySparse(T*& alpha_v2, int*& alpha_r2, int*& alpha_pB2, int*& alpha_pE2,
      double* alpha_v, mwSize* alpha_r, mwSize* alpha_pB, mwSize* alpha_pE, int M) 
{
   if (typeid(alpha_v) == typeid(alpha_v2)) {
      alpha_v2=reinterpret_cast<T*>(alpha_v);
   } else {
      alpha_v2 = new T[alpha_pB[M]];
      for (mwSize i = 0; i<alpha_pB[M]; ++i) alpha_v2[i] = static_cast<T>(alpha_v[i]);
   }
   if (typeid(alpha_r2) == typeid(alpha_r)) {
      alpha_r2=reinterpret_cast<int*>(alpha_r);
      alpha_pB2=reinterpret_cast<int*>(alpha_pB);
      alpha_pE2=reinterpret_cast<int*>(alpha_pE);
   } else {
      alpha_r2= new int[alpha_pB[M]];
      for (mwSize i = 0; i<alpha_pB[M]; ++i) alpha_r2[i]=static_cast<int>(alpha_r[i]);
      alpha_pB2= new int[M+1];
      for (int i = 0; i<=M; ++i) alpha_pB2[i]=static_cast<int>(alpha_pB[i]);
      alpha_pE2=alpha_pB2+1;
   }
};

/// Delete a sparse matrix which has been created using createCopySparse
template <typename T>
inline void 
DeleteCopySparse(T*& alpha_v2, int*& alpha_r2, int*& alpha_pB2, int*& alpha_pE2,
      double* alpha_v, mwSize* alpha_r) 
{
   if (typeid(alpha_v) != typeid(alpha_v2)) {
      delete[](alpha_v2);
   }
   if (typeid(alpha_r2) != typeid(alpha_r)) {
      delete[](alpha_r2);
      delete[](alpha_pB2);
   }
   alpha_v2=NULL;
   alpha_r2=NULL;
   alpha_pB2=NULL;
   alpha_pE2=NULL;
};

/// Create a m x n matrix
template <typename T>
inline mxArray* CreateMatrix(int m, int n);

/// Create a m x n double matrix
template <> 
inline mxArray* 
CreateMatrix<double>(int m, int n) 
{
   return mxCreateNumericMatrix(static_cast<mwSize>(m), static_cast<mwSize>(n),mxDOUBLE_CLASS,mxREAL);
}

/// Create a m x n float matrix
template <> 
inline mxArray* 
CreateMatrix<float>(int m, int n) 
{
   return mxCreateNumericMatrix(static_cast<mwSize>(m),static_cast<mwSize>(n),mxSINGLE_CLASS,mxREAL);
}

/// Create a h x w x V image
template <typename T>
inline mxArray* CreateImage(int h, int w, int V);

/// Create a h x w x V double image
template <>
inline mxArray* 
CreateImage<double>(int h, int w, int V) 
{
   if (V ==1) {
      return CreateMatrix<double>(h,w);
   } else {
      mwSize dims[3];
      dims[0]=h;
      dims[1]=w;
      dims[2]=V;
      return mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
   }
}

/// Create a h x w x V float image
template <>
inline mxArray* 
CreateImage<float>(int h, int w, int V) 
{
   if (V ==1) {
      return CreateMatrix<float>(h,w);
   } else {
      mwSize dims[3];
      dims[0]=h;
      dims[1]=w;
      dims[2]=V;
      return mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
   }
}

/// Create a h x w x V x dim image
template <typename T>
inline mxArray* 
Create4DImage(int h, int w, int V, int dim);

template <>
inline mxArray* 
Create4DImage<double>(int h, int w, int V, int dim) 
{
   if (dim ==1) {
      return CreateImage<double>(h,w,V);
   } else {
      mwSize dims[4];
      dims[0]=h;
      dims[1]=w;
      dims[2]=V;
      dims[3]=dim;
      return mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
   }
}

template <>
inline mxArray* 
Create4DImage<float>(int h, int w, int V, int dim) 
{
   if (dim ==1) {
      return CreateImage<float>(h,w,V);
   } else {
      mwSize dims[4];
      dims[0]=h;
      dims[1]=w;
      dims[2]=V;
      dims[3]=dim;
      return mxCreateNumericArray(4,dims,mxSINGLE_CLASS,mxREAL);
   }
}

/// Create a scalar
template <typename T> 
inline mxArray* 
CreateScalar() 
{
   return CreateMatrix<T>(1,1);
}

/// convert sparse matrix to Matlab sparse matrix
template <typename T> 
inline void 
ConvertSpMatrix(mxArray*& matlab_mat, int K,
      int M, int n, int nzmax, const T* v, const int* r, const int* pB) 
{
   matlab_mat=mxCreateSparse(K,M,nzmax,mxREAL);
   double* Pr=mxGetPr(matlab_mat);
   for (int i = 0; i<nzmax; ++i) Pr[i]=static_cast<double>(v[i]);
   mwSize* Ir=mxGetIr(matlab_mat);
   for (int i = 0; i<nzmax; ++i) Ir[i]=static_cast<mwSize>(r[i]);
   mwSize* Jc=mxGetJc(matlab_mat);
   if (n == 0) return;
   for (int i = 0; i<=n; ++i) Jc[i]=static_cast<mwSize>(pB[i]);
}

/// get a scalar from a struct
template <typename T> 
inline T 
GetScalarStruct(const mxArray* pr_struct,const char* name) 
{
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   if (!pr_field) {
      mexPrintf("Missing field: ");
      mexErrMsgTxt(name);
   }
   return static_cast<T>(mxGetScalar(pr_field));
}

/// get a scalar from a struct
inline void 
GetStringStruct(const mxArray* pr_struct, const char* name, char* field, const mwSize length) 
{
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   if (!pr_field) {
      mexPrintf("Missing field: ");
      mexErrMsgTxt(name);
   }
   mxGetString(pr_field,field,length);
}

inline void 
GetString(const mxArray* pr, std::string& str) 
{
  mwSize buflen = mxGetNumberOfElements(pr) + 1;
  char* buf = reinterpret_cast<char*>(mxCalloc(buflen, sizeof(char)));
  
  /* Copy the string data from string_array_ptr and place it into buf. */ 
  if (mxGetString(pr, buf, buflen) != 0)
    mexErrMsgTxt( "Could not convert string data.");
  str = std::string(buf);
}

inline std::string 
GetString(const mxArray* pr) 
{
  mwSize buflen = mxGetNumberOfElements(pr) + 1;
  char* buf = reinterpret_cast<char*>(mxCalloc(buflen, sizeof(char)));
  
  /* Copy the string data from string_array_ptr and place it into buf. */ 
  if (mxGetString(pr, buf, buflen) != 0)
    mexErrMsgTxt( "Could not convert string data.");
  return std::string(buf);
}

/// get a scalar from a struct
inline bool 
CheckField(const mxArray* pr_struct,const char* name) 
{
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   if (!pr_field) {
      mexPrintf("Missing field: ");
      mexPrintf(name);
      return false;
   }
   return true;
};

/// get a scalar from a struct  and provide a default value
template <typename T> 
inline T 
GetScalarStructDef(const mxArray* pr_struct, const char* name, const T def) 
{
  if (!pr_struct)
    return def;
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   return pr_field ? (T)(mxGetScalar(pr_field)) : def;
}

template <> 
inline std::string 
GetScalarStructDef(const mxArray* pr_struct, const char* name, const std::string def) 
{
  if (!pr_struct)
    return def;
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   return pr_field ? GetString(pr_field) : def;
}

inline mxArray* 
GetArrayStruct(const mxArray* pr_struct,const char* name) 
{
  if (!pr_struct)
    return NULL;
   mxArray *pr_field = mxGetField(pr_struct,0,name);
   return pr_field;
}

inline void 
super_flush(std::ostream& stream) 
{
    std::flush(stream);   
    mexEvalString("pause(0.0000000001);"); // to dump string.
}

template <class TMatrixType>
void
SaveMatrixToMatFile ( const TMatrixType& matrix, const int NumberRows, const int NumberColumns, const std::string fileName, const std::string varibleName )
{
  MATFile *pmat;
  pmat = matOpen(fileName.c_str(), "w");
  utlGlobalException(!pmat, "Error creating file " << fileName);
  int status; 

  mxArray *pa1= mxCreateDoubleMatrix(NumberRows,NumberColumns,mxREAL);
  utlGlobalException(!pa1, "Out of memory. Unable to create mxArray.");

  double * data = (double*)mxGetData(pa1);
  for ( int i = 0; i < NumberRows; i += 1 ) 
    for ( int j = 0; j < NumberColumns; j += 1 ) 
      data[j*NumberRows+i] = matrix(i,j);

  status = matPutVariableAsGlobal(pmat, varibleName.c_str(), pa1);
  utlGlobalException(status != 0, "Error using matPutVariableAsGlobal");

  mxDestroyArray(pa1);
  status = matClose(pmat);
  utlGlobalException(status!=0, "Error closing file " << fileName);
}

template <class TMatrixType>
void
ReadMatFileToMatrix ( const std::string fileName, const std::string varibleName, TMatrixType& matrix )
{
  MATFile *pmat;
  pmat = matOpen(fileName.c_str(), "r");
  utlGlobalException(!pmat, "Error creating file " << fileName);
  int status; 

  mxArray *pa1= matGetVariable(pmat, varibleName.c_str());
  utlGlobalException(!pa1, "Out of memory. Unable to create mxArray.");
  const mwSize* dims = mxGetDimensions(pa1);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  double * data = (double*)mxGetData(pa1);
  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      matrix(i,j) = data[j*row+i];

  mxDestroyArray(pa1);
  status = matClose(pmat);
  utlGlobalException(status!=0, "Error closing file " << fileName);
}

}

#endif
