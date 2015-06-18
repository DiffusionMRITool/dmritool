/**
 *       @file  SamplingSchemeDistance.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-28-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "SamplingSchemeDistanceCLP.h"
#include "utl.h"

#include <vnl/vnl_matrix.h>

/**
 * \brief  calculate the distances between to sampling schemes
 */
int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;

  typedef utl::NDArray<double,2> MatrixType;
  typedef utl_shared_ptr<MatrixType> MatrixPointer;
  MatrixPointer grad1, grad2;
  grad1 = utl::ReadGrad<double>(_InputFile1, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  grad2 = utl::ReadGrad<double>(_InputFile2, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);

  utlGlobalException(grad1->Rows()!=grad2->Rows(), "two file should have the same number of rows");

  // utl::PrintUtlMatrix(*grad1, "grad1");
  // utl::PrintUtlMatrix(*grad2, "grad2");
  int N = grad1->Rows();
  double maxDist=-1, distTmp=-1;
  int iMax=-1, jMax=-1;
  MatrixType dist(N,N,0.0);
  utl::NDArray<double,1> v1,v2;
  for ( int i = 0; i < N; i += 1 ) 
    {
    v1 = grad1->GetRow(i);
    // std::cout << "v1 = " << v1 << std::endl << std::flush;
    distTmp = std::numeric_limits<double>::max();
    int j = 0, jj=-1;
    for ( ; j < N; j += 1 ) 
      {
      v2 = grad2->GetRow(j);
      double tmp = utl::min( utl::ToVector<double>(v1-v2)->GetTwoNorm(), utl::ToVector<double>(v1+v2)->GetTwoNorm() );
      dist(i,j) = tmp;
      dist(j,i) = tmp;
      if (tmp < distTmp)
        {
        distTmp = tmp;
        jj = j;
        }
      // utlPrintVar4(true, i, j, tmp, distTmp);
      }
    if (distTmp > maxDist)
      {
      maxDist = distTmp;
      iMax = i;
      jMax = jj;
      }
    utlPrintVar3(_DebugArg.isSet(), i, distTmp, maxDist);
    utlPrintVar2(_DebugArg.isSet(), iMax, jMax);
    }

  std::cout << "maximal distance = " << maxDist << ", between " << iMax+1 << "-th sample ("<< grad1->GetRow(iMax) <<") in file " << _InputFile1 << 
    " and " << jMax+1 << "-th sample (" << grad2->GetRow(jMax) << ") in file " << _InputFile2 << std::endl << std::flush;
  
  return 0;
}
