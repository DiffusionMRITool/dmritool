/**
 *       @file  GetSPFBasisMatrix.cxx
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


#include <cstdio>
#include <iostream>
#include <cstdlib>

#include "mex.h" 
#include "utl.h"
#include "mexutils.h"
#include "mexSTD.h"
#include "utlMEX.h"

#include "itkSphericalPolarFourierImageFilter.h"



template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(!utl::mexCheckType<T>(prhs[2]), "type of argument 2 is not consistent");

  int shRank = mxGetScalar(prhs[0]);
  int raRank = mxGetScalar(prhs[1]);

  const mwSize* dimsOrientation = mxGetDimensions(prhs[2]);
  int row = static_cast<int>(dimsOrientation[0]);
  int column = static_cast<int>(dimsOrientation[1]);

  typedef utl::NDArray<double,2> MatrixType;

  utlException (column!=3, "orientation matrix should have 3 columns (x,y,z)");
  
  utl_shared_ptr<MatrixType > qOrientationMatrix( new MatrixType(row, column));
  utl::GetUtlMatrixFromMXArray(prhs[2], qOrientationMatrix.get());
  *qOrientationMatrix = utl::CartesianToSpherical(*qOrientationMatrix); // spherical format is used 

  const mwSize* dimsBVectors = mxGetDimensions(prhs[3]);
  int rowB = static_cast<int>(dimsBVectors[0]);
  int columnB = static_cast<int>(dimsBVectors[1]);

  // utlPrintVar2(true, rowB, columnB);
  utlException(columnB!=1 && rowB!=1, "orientation matrix should have 1 column or 1 row");
  utlException((columnB==1 && rowB!=row) || (rowB==1 && columnB!=column), "orientation matrix and B vector should have the same number of rows");

  utl_shared_ptr<std::vector<double> > bVector (new std::vector<double>());
  utl::GetSTDVectorFromMXArray( prhs[3], bVector.get());
  

  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3>  VectorImageType;
  typedef itk::SphericalPolarFourierImageFilter<VectorImageType, VectorImageType> SPFIFilterBaseType;
  typename SPFIFilterBaseType::Pointer spfiFilter = SPFIFilterBaseType::New();
  spfiFilter->GetSamplingSchemeQSpace()->SetBVector(bVector);
  spfiFilter->GetSamplingSchemeQSpace()->SetOrientationsSpherical(qOrientationMatrix);
  spfiFilter->SetSHRank(shRank);
  spfiFilter->SetRadialRank(raRank);
  bool originalBasis = true;
  if (nrhs>4)
    {
    utlGlobalException(mxGetClassID(prhs[4])!=mxSTRUCT_CLASS, "the fifth input has to be a struct");
    T tau = utl::GetScalarStructDef<T>(prhs[4],"tau", ONE_OVER_4_PI_2);  //  4*pi*pi*tau=1
    spfiFilter->GetSamplingSchemeQSpace()->SetTau(tau);

    T scale = utl::GetScalarStructDef<T>(prhs[4],"scale", -1.0);
    if (scale>0)
      spfiFilter->SetBasisScale(scale);

    originalBasis = utl::GetScalarStructDef<bool>(prhs[4],"original", true);  // 0: original SPF matrix, 1: independent SPF basis matrix
    }


  utl_shared_ptr<MatrixType > basisMatrix (new MatrixType());
  // spfiFilter->SetDebug(true);
  spfiFilter->ComputeBasisMatrix();
  basisMatrix = spfiFilter->GetBasisMatrix();
  // if (basisMatrix)
  //   {
  //   utl::PrintUtlMatrix(*basisMatrix, "basisMatrix");
  //   }

  if (originalBasis)
    {
    utl::GetMXArrayFromUtlMatrix(basisMatrix.get(), plhs[0]);
    }
  else
    {

    int n_b_sh = (shRank+1)*(shRank+2)/2;
    int n_b_ra = raRank + 1;

    MatrixType basisMatrixIndependent;

    spfiFilter->ComputeRadialVectorForE0InBasis();
    utl_shared_ptr<std::vector<double> > R_n_0_vec = spfiFilter->GetGn0();
    basisMatrixIndependent.ReSize(basisMatrix->Rows(), n_b_sh*(n_b_ra-1));
    // utlException(R_n_0_vec[0]==0, "it should be not zero!, R_n_0_vec[0]="<< R_n_0_vec[0]);

    for ( int i = 0; i < raRank; i += 1 ) 
      {
      for ( int j = 0; j < n_b_sh; j += 1 ) 
        {
        for ( int ss = 0; ss < basisMatrixIndependent.Rows(); ss += 1 ) 
          {
          basisMatrixIndependent(ss,i*n_b_sh+j) = (*basisMatrix)(ss,(i+1)*n_b_sh+j) - (*R_n_0_vec)[i+1]/(*R_n_0_vec)[0] * (*basisMatrix)(ss,j);
          }
        }
      }

    utl::GetMXArrayFromUtlMatrix(&basisMatrixIndependent, plhs[0]);
    }


  return;
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=4 && nrhs!=5, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=1, "Bad number of outputs arguments");

  // if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS) 
  //   callFunction<float>(plhs,prhs,nlhs,nrhs);
  // else
    callFunction<double>(plhs,prhs,nlhs,nrhs);
} 

