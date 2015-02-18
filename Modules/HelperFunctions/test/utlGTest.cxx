/**
 *       @file  utlGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-11-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "utl.h"

typedef utl::NDArray<double,2>         MatrixType;
typedef utl_shared_ptr<MatrixType> MatrixPointer;

const MatrixPointer dirctions = utl::ReadGrad<double>(3);



ITK_THREAD_RETURN_TYPE __ThreadedMethod(void* arg)
{
  const itk::MultiThreader::ThreadInfoStruct* threadInfo = static_cast<itk::MultiThreader::ThreadInfoStruct*>(arg);
  if (threadInfo)
    {
    const unsigned int threadId = threadInfo->ThreadID;
    std::vector<MatrixPointer>* dataVec = static_cast<std::vector<MatrixPointer>*>(threadInfo->UserData);
    if (dataVec)
      {
      (*dataVec)[threadId] = utl::ComputeSHMatrix<double>(4, *dirctions, SPHERICAL_TO_SPHERICAL);
      MatrixPointer mat = (*dataVec)[threadId];
      } 
    else 
      {
      std::cerr << "ERROR: UserData was not of type ThreadDataVec*" << std::endl;
      return ITK_THREAD_RETURN_VALUE;
      }
    } 
  else 
    {
    std::cerr << "ERROR: arg was not of type itk::MultiThreader::ThreadInfoStruct*" << std::endl;
    return ITK_THREAD_RETURN_VALUE;
    }
  return ITK_THREAD_RETURN_VALUE;
}


TEST(utl, ComputeSHMatrix_MultiThreads)
{
  itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
  int numberOfThreads = 8;
  std::vector<MatrixPointer> matVec(numberOfThreads);

  threader->SetGlobalMaximumNumberOfThreads(numberOfThreads+10);
  threader->SetNumberOfThreads(numberOfThreads);
  threader->SetSingleMethod(__ThreadedMethod, &matVec);
  threader->SingleMethodExecute();
  for ( int i = 0; i < numberOfThreads-1; i += 1 ) 
    {
    MatrixType diff = *matVec.back()- *matVec[i];
    EXPECT_NEAR(0.0, diff.AbsoluteMaxValue(), 1e-15 );
    }
}


