/**
 *       @file  utlMEX.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-07-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlMEX_h
#define __utlMEX_h

#include "mexutils.h"
#include "mexSTD.h"
#include "mexVNL.h"
#include "mexITK.h"

#include "utlNDArray.h"

namespace utl
{

template <class T>
inline void 
GetUtlMatrixFromMXArray ( const mxArray* pr, NDArray<T,2>* mat )
{
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  double * data = mxGetPr(pr);
  if (mat->Rows()!=row || mat->Columns()!=column )
    mat->ReSize(row, column);

  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      (*mat)(i,j) = data[j*row+i];
}

template <class T>
inline void 
GetMXArrayFromUtlMatrix ( const NDArray<T,2>* mat, mxArray*& pr )
{
  int row = mat->Rows();
  int column = mat->Columns();
  pr = CreateMatrix<double>(row, column);
  
  double * data = mxGetPr(pr);
  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      data[j*row+i] = (*mat)(i,j);
}

template <class T>
inline void 
GetUtlVectorFromMXArray ( const mxArray* pr, NDArray<T,1>* mat )
{
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  utlException(row>1 && column>1, "mxArray is not a vector");
  int size = utl::max(row, column);
  double * data = mxGetPr(pr);
  if (mat->Size()!=size )
    mat->ReSize(size);

  for ( int i = 0; i < size; i += 1 ) 
    (*mat)[i] = data[i];
}

template <class T>
inline void 
GetMXArrayFromUtlVector ( const NDArray<T,1>* mat, mxArray*& pr )
{
  int row = mat->Size();
  pr = CreateMatrix<double>(row, 1);
  
  double * data = mxGetPr(pr);
  for ( int i = 0; i < row; i += 1 ) 
      data[i] = (*mat)[i];
}


}


#endif 
