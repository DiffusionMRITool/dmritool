/**
 *       @file  utl.h
 *      @brief  helper functions specifically used in dmritool
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "10-31-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utl_h
#define __utl_h


#include "DMRITOOLConfigure.h"

#ifndef UTL_MAX_THREADS
#define UTL_MAX_THREADS 8
#endif

#ifdef UTL_USE_OPENMP
#include <omp.h>
#include "utlOpenMP.h"
#endif

#ifdef UTL_USE_MKL
#include "utlMKL.h"
#endif

#include "utlSTDHeaders.h"

#include "utlCore.h"
#include "utlMath.h"
#include "utlVNL.h"
#include "utlVNLIO.h"
#include "utlITK.h"
#include "utlVNLBlas.h"
#include "utlVNLLapack.h"
#include "utlNDArray.h"

#include "utlDMRI.h"

#include "itkSphericalHarmonicsGenerator.h"

#include "itkMultiVolumeImageToVectorImageFilter.h"

namespace itk
{


/** if it is a VectorImage, directly read it. If it is a 4D image, read it as Image, then convert it to a VectorImage  */
template <class PixelType>
inline void 
ReadVectorImage ( const std::string filename, SmartPointer<VectorImage<PixelType,3> > &image, const std::string printInfo="Reading Image:" )
{
  typedef VectorImage<PixelType,3> VectorImageType;
  if (itk::IsVectorImage(filename))
    itk::ReadImage<VectorImageType>(filename, image, printInfo);
  else
    {
    typedef itk::Image<PixelType, 4> MultiVolumeImageType;
    typename MultiVolumeImageType::Pointer MultiVolumeImage;
    itk::ReadImage<MultiVolumeImageType>(filename, MultiVolumeImage, printInfo);

    typedef itk::MultiVolumeImageToVectorImageFilter<PixelType, PixelType> MultiVolumeToVectorFilter;
    typename MultiVolumeToVectorFilter::Pointer filter = MultiVolumeToVectorFilter::New();
    filter->SetInput(MultiVolumeImage);
    filter->Update();
    image = filter->GetOutput();
    }
}

}



namespace utl 
{

// [>* The table of integration of triple SH basis (real, thesis), genegrated by print_sh_integration  <]
// const static std::string SH3Itegralhdr = std::string(SH3Integral_HDR);

// [>* gradients file for tess=1 with 6 directions  <]
// const static std::string DirectionsT1 = std::string(Path_Directions_T1);
// [>* gradients file for tess=2 with 21 directions  <]
// const static std::string DirectionsT2 = std::string(Path_Directions_T2);
// [>* gradients file for tess=3 with 81 directions  <]
// const static std::string DirectionsT3 = std::string(Path_Directions_T3);
// [>* gradients file for tess=4 with 321 directions  <]
// const static std::string DirectionsT4 = std::string(Path_Directions_T4);
// [>* gradients file for tess=5 with 1281 directions  <]
// const static std::string DirectionsT5 = std::string(Path_Directions_T5);
// [>* gradients file for tess=6 with 5121 directions  <]
// const static std::string DirectionsT6 = std::string(Path_Directions_T6);
// [>* gradients file for tess=7 with 20481 directions  <]
// const static std::string DirectionsT7 = std::string(Path_Directions_T7);


template <class T>
itk::VariableLengthVector<T>
UtlVectorToVariableLengthVector ( const NDArray<T,1>& vec )
{
  itk::VariableLengthVector<T> v(vec.Size());
  for ( int i = 0; i < vec.Size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
NDArray<T,1> 
VariableLengthVectorToUtlVector ( const itk::VariableLengthVector<T>& vec )
{
  NDArray<T,1> v(vec.GetSize());
  for ( int i = 0; i < vec.GetSize(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
inline utl_shared_ptr<NDArray<T,2> >
ReadGrad(const int tess=3, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode= CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
{
  switch ( tess )
    {
  case 1 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT1), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 2 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT2), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 3 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT3), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 4 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT4), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 5 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT5), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 6 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT6), NoSymmetricDuple, mode, flipx, flipy, flipz);
  case 7 :
    return ReadGrad<T>(CreateExpandedPath(DirectionsT7), NoSymmetricDuple, mode, flipx, flipy, flipz);
  default :
    utlAssert(false, "tess should be 1, 2, 3, 4, 5, 6, 7");
    return utl_shared_ptr<NDArray<T,2> > ();
    }
}

template <class T>
inline utl_shared_ptr<NDArray<T,2> >
ReadGradElectricRepulsion(const int num, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode= CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
{
  char buf[255];
  if (num<10)
    sprintf(buf, "00%d", num);
  else if (num<100)
    sprintf(buf, "0%d", num);
  else
    sprintf(buf, "%d", num);
  std::string index (buf);
  std::string file= CreateExpandedPath(GradientsElec) + std::string("/Elec") + index + std::string(".txt");
  return ReadGrad<T>(file, NoSymmetricDuple, mode, flipx, flipy, flipz, need_normalize);
}

// [>* generate SH basis with given rank and the gradients <]
// template <class T>
// inline utl_shared_ptr<vnl_matrix<T> >
// ComputeSHMatrix ( const unsigned int rank, const vnl_matrix<T>& grad, const int mode)
// {
//   utlException(grad.rows()==0 || grad.columns()!=3, "wrong size of gradients!");
//   int numberOfBasisFunctions = (rank + 1)*(rank + 2)/2;
//   int numberOfDirections = grad.rows();
//   vnl_matrix<T> gradSpherical;
//   if (mode==CARTESIAN_TO_SPHERICAL)
//     gradSpherical = CartesianToSpherical(grad);
//   else if (mode==SPHERICAL_TO_SPHERICAL)
//     ;
//   else
//     utlGlobalException(true, "wrong mode");

//   utl_shared_ptr<vnl_matrix<T> > BMatrix (new vnl_matrix<T>(numberOfDirections, numberOfBasisFunctions));

//   // itk::SphericalHarmonicsGenerator<double>::Pointer sh = itk::SphericalHarmonicsGenerator<double>::New();
//   int jj=0;
//   for ( unsigned int k = 0; k < numberOfDirections; k++ )
//     {
//     jj = 0;
//     for ( int l = 0; l <= rank; l += 2 )
//       for ( int m = -l; m <= l; m++ )
//         {
//         if (mode==SPHERICAL_TO_SPHERICAL)
//           (*BMatrix)(k,jj) = itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,grad(k,1),grad(k,2));
//         else
//           (*BMatrix)(k,jj) = itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,gradSpherical(k,1),gradSpherical(k,2));
//         jj++;
//         }
//     }
//   return BMatrix;
// }

/** generate SH basis with given rank and the gradients */
template <class T>
inline utl_shared_ptr<NDArray<T,2> >
ComputeSHMatrix ( const unsigned int rank, const NDArray<T,2>& grad, const int mode)
{
  utlException(grad.Rows()==0 || grad.Columns()!=3, "wrong size of gradients!");
  int numberOfBasisFunctions = (rank + 1)*(rank + 2)/2;
  int numberOfDirections = grad.Rows();
  NDArray<T,2> gradSpherical;
  if (mode==CARTESIAN_TO_SPHERICAL)
    gradSpherical = CartesianToSpherical(grad);
  else if (mode==SPHERICAL_TO_SPHERICAL)
    ;
  else
    utlGlobalException(true, "wrong mode");

  utl_shared_ptr<NDArray<T,2> > BMatrix (new NDArray<T,2>(numberOfDirections, numberOfBasisFunctions));

  // itk::SphericalHarmonicsGenerator<double>::Pointer sh = itk::SphericalHarmonicsGenerator<double>::New();
  int jj=0;
  for ( unsigned int k = 0; k < numberOfDirections; k++ )
    {
    jj = 0;
    for ( int l = 0; l <= rank; l += 2 )
      for ( int m = -l; m <= l; m++ )
        {
        if (mode==SPHERICAL_TO_SPHERICAL)
          (*BMatrix)(k,jj) = itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,grad(k,1),grad(k,2));
        else
          (*BMatrix)(k,jj) = itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,gradSpherical(k,1),gradSpherical(k,2));
        jj++;
        }
    }
  return BMatrix;
}

// template < class T >
// void
// MatchBVectorAndGradientMatrix (const T& br, std::vector<T>& vec, const vnl_matrix<T>& grad )
// {
//   utlException(grad.rows()==0, "grad.size()==0");
//   vec = std::vector<T>(grad.rows(), br);
// }

// template < class T >
// void
// MatchBVectorAndGradientMatrix ( std::vector<T>& vec, vnl_matrix<T>& grad )
// {
//   utlException(vec.size()==0 || grad.rows()==0, "wrong size! vec.size()="<<vec.size()<<", grad.rows()="<<grad.rows());
  
//   if (grad.rows()==vec.size())
//     return;

//   std::vector<T> vec_result;
//   if (vec.size()==1)
//     {
//     MatchBVectorAndGradientMatrix(vec[0], vec_result, grad);
//     vec = vec_result;
//     return;
//     }

//   vec_result.resize(vec.size()* grad.rows());
//   vnl_matrix<T> grad_result(vec.size()* grad.rows(), grad.columns());
//   for (int i = 0; i<vec.size(); i++)
//     {
//     for ( int j = 0; j < grad.rows(); j += 1 ) 
//       {
//       vec_result[i*grad.rows()+j] = vec[i];
//       for ( int kk = 0; kk < 3; kk += 1 ) 
//         grad_result(i*grad.rows()+j,kk) = grad(j,kk);
//       }
//     }
//   vec = vec_result;
//   grad = grad_result;
//   return;
// }

template <class PointsContainer, class VnlValueType>
void 
PointsContainerToUtlMatrix ( const PointsContainer& points, utl::NDArray<VnlValueType,2>& matrix )
{
  matrix.Clear();
  typedef typename PointsContainer::ConstIterator PointsIterator;
  typedef typename PointsContainer::Element PointType;

  const unsigned int pointDimension = PointType::PointDimension;
  unsigned int numberOfPoints = points.Size();
  if (numberOfPoints==0)
    return;

  matrix.ReSize(numberOfPoints, pointDimension);
  PointsIterator iterator = points.Begin();
  PointsIterator end = points.End();
  
  unsigned int count=0;
  while( iterator != end )
    {
    PointType orientation = iterator.Value();
    for (unsigned int k=0; k<pointDimension; k++)
      matrix(count,k) = orientation[k];
    iterator++;
    count++;
    }
}

template <class VnlValueType, class PointsContainer>
void 
UtlMatrixToPointsContainer ( const NDArray<VnlValueType,2>& matrix, PointsContainer& points )
{
  points.Initialize();
  typedef typename PointsContainer::Element PointType;

  const unsigned int pointDimension = PointType::PointDimension;
  utlGlobalException(pointDimension!=matrix.Columns(), "wrong size of matrix or point dimension. pointDimension="<< pointDimension << ", matrix.Columns()="<<matrix.Columns());

  for ( int i = 0; i < matrix.Rows(); i += 1 ) 
    {
    PointType point;
    for ( int j = 0; j < matrix.Columns(); j += 1 ) 
      {
      point[j] = matrix(i,j);
      }
    points.InsertElement(i, point);
    }
}



static inline int 
InitializeMKL(const int numThreads) 
{
  int NUM_THREADS;
#ifdef UTL_USE_MKL
  NUM_THREADS = (numThreads <= 0 ) ? std::min(UTL_MAX_THREADS,mkl_get_max_threads()) : numThreads;
  mkl_domain_set_num_threads ( NUM_THREADS, MKL_DOMAIN_ALL );
  mkl_set_num_threads(NUM_THREADS);
  mkl_set_dynamic(0);
#else
  NUM_THREADS = 1;
#endif
  return NUM_THREADS;
}

static inline int 
InitializeOpenMP(const int numThreads) 
{
  int NUM_THREADS;
#ifdef UTL_USE_OPENMP
  NUM_THREADS = (numThreads <= 0) ? std::min(UTL_MAX_THREADS,omp_get_num_procs()) : numThreads;
  omp_set_nested(0);
  omp_set_dynamic(0);
  omp_set_num_threads(NUM_THREADS);
#else
  NUM_THREADS = 1;
#endif
  return NUM_THREADS;
}

static inline void 
InitializeThreadedLibraries(const int numThreads)
{
  utl::InitializeMKL((numThreads>1 || numThreads<=0)?1:-1);
  utl::InitializeOpenMP((numThreads>1 || numThreads<=0)?1:-1);
}

}

#endif 
