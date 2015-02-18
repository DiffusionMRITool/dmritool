/**
 *       @file  spamsMatrixTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-02-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlVNL.h"
#include "linalg.h"
#include "utlITKSpams.h"



/**
 * \brief  test Vector Matrix in spams
 */
int 
main (int argc, char const* argv[])
{
  
  // test matrix
  double mmatrix [] =
    {
   2.2700 ,  -1.5400  ,  1.1500  , -1.9400,
   0.2800 ,  -1.6700  ,  0.9400  , -0.7800,
  -0.4800 ,  -3.0900  ,  0.9900  , -0.2100,
   1.0700 ,   1.2200  ,  0.7900  ,  0.6300,
  -2.3500 ,   2.9300  , -1.4500  ,  2.3000,
   0.6200 ,  -7.3900  ,  1.0300  , -2.5700
    };


  vnl_matrix<double> M (mmatrix, 6, 4);
  M.print(std::cout<<"M:\n");

  spams::Matrix<double> DM(mmatrix, 4,6);
  DM.print("DM");

  DM.setData(mmatrix, 4,6);
  DM.print("DM2");

  spams::VnlMatrixToMatrix(M, DM);
  DM.print("DM3");

  spams::MatrixToVnlMatrix(DM, M);
  M.print(std::cout<<"M:\n");

  spams::Vector<double> v;
  DM.copyCol(2, v);
  v.print("v");

  vnl_vector<double> v2;
  spams::VectorToVnlVector(v, v2);
  utl::PrintVnlVector(v2, "v2");

  return 0;
}

