/**
 *       @file  mexSTD.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "09-27-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __mexSTD_h
#define __mexSTD_h

#include <mex.h>
#include "mexutils.h"
#include "utlCoreMacro.h"

namespace utl
{

template <class T>
inline void 
GetSTDVectorFromMXArray ( const mxArray* pr, std::vector<T>* vec )
{
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  double * data = mxGetPr(pr);

  if (row==1)
    {
    vec->resize(column);
    for ( int i = 0; i < column; i += 1 ) 
      (*vec)[i] = data[i];
    }
  else if (column==1)
    {
    vec->resize(row);
    for ( int i = 0; i < row; i += 1 ) 
      (*vec)[i] = data[i];
    }
  else
    utlException(true, "the matrix should be a vector");
}

/** convert a std::vector to a colume vector  */
template <class T>
inline void 
GetMXArrayFromSTDVector ( const std::vector<T>* vec, mxArray*& pr  )
{
  utlException (!vec || vec && vec->size()==0, "the vector is null");

  int row = vec->size();
  int column = 1;
  
  pr = CreateMatrix<double>(row, 1);
  double * data = mxGetPr(pr);
  for ( int i = 0; i < row; i += 1 ) 
    data[i] = (*vec)[i];
}

}

#endif 

