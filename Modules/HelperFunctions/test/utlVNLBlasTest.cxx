/**
 *       @file  utlVNLBlasTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlCore.h"
#include "utlVNL.h"
#include "utlVNLBlas.h"


/**
 * \brief  test dgemm in blas
 */
int 
main (int argc, char const* argv[])
{
    {
    vnl_matrix<double> mat1, mat2;
    vnl_vector<double> vec1, vec2;

    const double array_1 [] =
      {
      1,  1, 
      -1, 2,
      2,  1
      };
    const double array_2 [] =
      {
      1,  1, 3, 4 ,
      -1, 2, 5, 6
      };
    vec1.set_size(2);
    vec1[0]=1; vec1[1]=4;

    mat1 = vnl_matrix<double>(array_1, 3, 2);
    mat2 = vnl_matrix<double>(array_2, 2, 4);

      {
      std::cout << "\nM * M" << std::endl << std::flush;
      utl::PrintVnlMatrix(mat1, "mat1");
      utl::PrintVnlMatrix(mat2, "mat2");
      vnl_matrix<double> mat_mul = mat1*mat2;
      utl::PrintVnlMatrix(mat_mul, "mat1*mat2");
      utl::ProductVnlMM(mat1, mat2, mat_mul);
      utl::PrintVnlMatrix(mat_mul, "ProductMM");
      vnl_matrix<double> mat2t = mat2.transpose();
      mat_mul = mat1*mat2t.transpose();
      utl::PrintVnlMatrix(mat_mul, "mat1 * mat2^T");
      utl::ProductVnlMMt(mat1, mat2t, mat_mul);
      utl::PrintVnlMatrix(mat_mul, "ProductMMt");
      }

      {
      std::cout << "\nM * v" << std::endl << std::flush;
      utl::PrintVnlMatrix(mat1, "mat1");
      utl::PrintVnlVector(vec1, "vec1");
      vnl_vector<double> vec_mul = mat1*vec1;
      utl::PrintVnlVector(vec_mul, "mat1*vec1");
      utl::ProductVnlMv(mat1, vec1, vec_mul);
      utl::PrintVnlVector(vec_mul, "ProductMv");
      }

      {
      std::cout << "\nv * M" << std::endl << std::flush;
      utl::PrintVnlMatrix(mat2, "mat2");
      utl::PrintVnlVector(vec1, "vec1");
      vnl_vector<double> vec_mul = vec1*mat2;
      utl::PrintVnlVector(vec_mul, "vec1*mat2");
      utl::ProductVnlvM(vec1, mat2, vec_mul);
      utl::PrintVnlVector(vec_mul, "ProductbA");
      vnl_matrix<double> mat2t = mat2.transpose();
      utl::ProductVnlvMt(vec1, mat2t, vec_mul);
      utl::PrintVnlVector(vec_mul, "ProductvMt");
      }

    }

    {
      std::cout << "\nM * M time cost" << std::endl << std::flush;
    vnl_matrix<double> mat1(514,180), mat2(180,254);
    for ( int i = 0; i < mat1.rows(); i += 1 ) 
      for ( int j = 0; j < mat1.columns(); j += 1 ) 
        mat1(i,j) = utl::Random<double>(-2.0,2.0);
    for ( int i = 0; i < mat2.rows(); i += 1 ) 
      for ( int j = 0; j < mat2.columns(); j += 1 ) 
        mat2(i,j) = utl::Random<double>(-2.0,2.0);

    int N = 100;
      {
      utl::Tic(std::cout<<"mat1*mat2 start \n");
      for ( int i = 0; i < N; i += 1 ) 
        {
        vnl_matrix<double> tmp;
        tmp = mat1*mat2;
        }
      utl::Toc();
      }
      {
      utl::Tic(std::cout<<"ProductVnlMM \n");
      for ( int i = 0; i < N; i += 1 ) 
        {
        vnl_matrix<double> tmp;
        utl::ProductVnlMM(mat1, mat2, tmp);
        }
      utl::Toc();
      }
    }

  return 0;
}

