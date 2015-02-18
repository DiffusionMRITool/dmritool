/**
 *       @file  mexVNL.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-04-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __mexVNL_h
#define __mexVNL_h

#include <mex.h>
#include <vnl/vnl_matrix.h>
#include "mexutils.h"

#include "utlCore.h"

namespace utl 
{

template <class T>
inline void 
GetVNLMatrixFromMXArray ( const mxArray* pr, vnl_matrix<T>* mat )
{
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  double * data = mxGetPr(pr);
  if (mat->rows()!=row || mat->columns()!=column )
    mat->set_size(row, column);

  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      (*mat)[i][j] = data[j*row+i];
}

template <class T>
inline void 
GetMXArrayFromVNLMatrix ( const vnl_matrix<T>* mat, mxArray*& pr )
{
  int row = mat->rows();
  int column = mat->columns();
  pr = CreateMatrix<double>(row, column);
  
  double * data = mxGetPr(pr);
  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      data[j*row+i] = (*mat)[i][j];
}

template <class T>
inline void 
GetVNLVectorFromMXArray ( const mxArray* pr, vnl_vector<T>* mat )
{
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  utlException(row>1 && column>1, "mxArray is not a vector");
  int size = utl::max(row, column);
  double * data = mxGetPr(pr);
  if (mat->size()!=size )
    mat->set_size(size);

  for ( int i = 0; i < size; i += 1 ) 
    (*mat)[i] = data[i];
}

template <class T>
inline void 
GetMXArrayFromVNLVector ( const vnl_vector<T>* mat, mxArray*& pr )
{
  int row = mat->size();
  pr = CreateMatrix<double>(row, 1);
  
  double * data = mxGetPr(pr);
  for ( int i = 0; i < row; i += 1 ) 
      data[i] = (*mat)[i];
}

template <class T>
inline void
SaveVNLMatrixToMatFile ( const vnl_matrix<T>* mat, const std::string fileName, const std::string varibleName )
{
  SaveMatrixToMatFile<vnl_matrix<T> >(*mat, mat->rows(), mat->columns(), fileName, varibleName);
}

template <class T>
inline void
ReadMatFileToVNLMatrix ( const std::string fileName, const std::string varibleName, vnl_matrix<T>* mat )
{
  MATFile *pmat;
  pmat = matOpen(fileName.c_str(), "r");
  utlGlobalException(!pmat, "Error creating file " << fileName);
  int status; 

  mxArray *pr= matGetVariable(pmat, varibleName.c_str());
  utlGlobalException(!pr, "Out of memory. Unable to create mxArray.");
  const mwSize* dims = mxGetDimensions(pr);
  int row = static_cast<int>(dims[0]);
  int column = static_cast<int>(dims[1]);
  
  if (mat->rows()!=row || mat->columns()!=column )
    mat->set_size(row, column);

  double * data = (double*)mxGetData(pr);
  for ( int i = 0; i < row; i += 1 ) 
    for ( int j = 0; j < column; j += 1 ) 
      (*mat)(i,j) = data[j*row+i];
}

}

#endif 

