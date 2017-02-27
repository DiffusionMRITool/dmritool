/**
 *       @file  utlDMRIStoredTables.h
 *      @brief  
 *     Created  "11-12-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlDMRIStoredTables_h
#define __utlDMRIStoredTables_h

#include "itkUnaryFunctorLookUpTable.h"

#include "DMRITOOLConfigure.h"
#include "utlDMRI.h"
#include "utlCoreMacro.h"
#include "utlITK.h"
#include "utlCore.h"

// #include "itkSphericalHarmonicsGenerator.h"
#include "utlNDArray.h"



namespace itk
{
template<typename T> class SphericalHarmonicsGenerator;


typedef UnaryFunctorLookUpTable<utl::Functor::Exp<double> >       LUTExpType;
typedef LUTExpType::Pointer                                       LUTExpPointer;
/** use lookup table to approximate exp(-x)  */
static  LUTExpPointer lutExp = LUTExpType::New();

inline void 
InitializeLUTExp()
{
  if (!lutExp->IsTableBuilt())
    {
    lutExp = LUTExpType::New();
    lutExp->SetVariableMax(0);
    lutExp->SetVariableMin(-30);
    lutExp->SetNumberOfBins(30*5e3);
    lutExp->BuildTable();
    }
}

inline double 
lutExpValue(const double x)
{
  if (x<=0)
    return lutExp->GetFunctionValue(x);
  else 
    return std::exp(x);
}

}

namespace utl
{

static utl_shared_ptr<NDArray<double,3> > SH3IntegralTable(new NDArray<double,3>);


/**
 *   \class   GradientTable
 *   \brief   gradient table related functions and stored tables.
 *   \author  Jian Cheng
 *   \date    10-03-2016
 */
template < class T=double >
class GradientTable 
  {
public:

  typedef GradientTable<T>                           Self;
  typedef utl_shared_ptr<GradientTable<T> >          Pointer;
  
  typedef utl::NDArray<T,1>                          VectorType;
  typedef utl::NDArray<T,2>                          MatrixType;
  typedef utl_shared_ptr<MatrixType>                 MatrixPointer;

  /** GradientTable constructor  */
  GradientTable ()
    {
    }

  /** GradientTable deconstructor  */
  virtual ~GradientTable ()
    {}

  static Pointer GetInstance()
    {
    if (!m_Instance)
      {
      m_Instance = utl_shared_ptr<Self>(new Self());
      if (m_Instance->m_GradTable.size()!=7)
        {
        m_Instance->m_GradTable.resize(7);
        for ( int i = 0; i < 7; ++i ) 
          m_Instance->m_GradTable[i] = MatrixPointer(new MatrixType());
        }
      }
    return m_Instance;
    }

  /** Initially read stored gradient tables  */
  static void 
  Initialize(const int tessorder)
    {
    utlSAGlobalException(tessorder<1 || tessorder>7)(tessorder).msg("wrong tess order. tessorder should be in [1,7]");
    Pointer instance = GetInstance();
    MatrixPointer mat = instance->m_GradTable[tessorder-1];
    switch ( tessorder )
      {
    case 1 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT1), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 2 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT2), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 3 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT3), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 4 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT4), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 5 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT5), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 6 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT6), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    case 7 : {  if (mat->Rows()==0) mat = utl::ReadGrad<double>(CreateExpandedPath(DirectionsT7), DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); break;    }
    default : utlSAGlobalException(true)(tessorder).msg("wrong tess order. tessorder should be in [1,7]");   break;
      }
    instance->m_GradTable[tessorder-1] = mat;
    }

  /** Get a copy of the stored gradient table, considering the output may be modified later.  */
  static MatrixPointer
  GetGrad(const int tessorder, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode= CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
    {
    utlSAGlobalException(tessorder<1 || tessorder>7)(tessorder).msg("wrong tess order. tessorder should be in [1,7]");
    utlSAException(mode==SPHERICAL_TO_CARTESIAN || mode==SPHERICAL_TO_SPHERICAL)(mode).msg("wrong mode. The stored gradient tables are in cartesian mode");

    Pointer instance = GetInstance();
    MatrixPointer grad = instance->m_GradTable[tessorder-1];
    utlException(grad->Size()==0, "call Initialize first to read grad");
    MatrixPointer mat(new MatrixType());
    if (NoSymmetricDuple==DIRECTION_DUPLICATE)
      mat->ReSize(2*grad->Rows(),3);
    else
      mat->ReSize(grad->Rows(),3);

    for ( int i=0, j=0; i < grad->Rows(); ++i, ++j ) 
      {
      (*mat)(j,0) = flipx==DIRECTION_FLIP ? -(*grad)(i,0) : (*grad)(i,0);
      (*mat)(j,1) = flipy==DIRECTION_FLIP ? -(*grad)(i,1) : (*grad)(i,1);
      (*mat)(j,2) = flipz==DIRECTION_FLIP ? -(*grad)(i,2) : (*grad)(i,2);

      if (NoSymmetricDuple==DIRECTION_DUPLICATE)
        {
        j++;
        (*mat)(j,0) = -(*mat)(j-1,0);
        (*mat)(j,1) = -(*mat)(j-1,1);
        (*mat)(j,2) = -(*mat)(j-1,2);
        }
      }

    if (mode==CARTESIAN_TO_SPHERICAL)
      *mat = utl::CartesianToSpherical(*mat);

    return mat;
    }


protected:

  std::vector<MatrixPointer> m_GradTable;

  static Pointer m_Instance;

private:
  /** GradientTable constructor  */
  GradientTable ( const GradientTable &other );
  /** GradientTable operator = , assignment operator */
  GradientTable & operator = ( const GradientTable &other );

  }; 

template<>
typename GradientTable<double>::Pointer GradientTable<double>::m_Instance = NULL;  


template <class T>
inline utl_shared_ptr<NDArray<T,2> >
ReadGrad(const int tess, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode= CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
{
  utl::GradientTable<double>::Initialize(tess);
  utl_shared_ptr< utl::Matrix<double> > grad = utl::GradientTable<double>::GetGrad(tess, NoSymmetricDuple, mode, flipx, flipy, flipz);
  return grad;
}

template <>
inline utl_shared_ptr<NDArray<float,2> >
ReadGrad<float>(const int tess, const int NoSymmetricDuple, const int mode, 
  const int flipx, const int flipy, const int flipz, const bool need_normalize) 
{
  utl_shared_ptr< utl::Matrix<double> > grad = utl::ReadGrad<double>(tess, NoSymmetricDuple, mode, flipx, flipy, flipz);
  utl_shared_ptr<utl::Matrix<float> > matOut(new utl::Matrix<float>());
  *matOut = *grad;
  return matOut;
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

/** Initialize SH3IntegralTable with 3 dimensional sizes from 3 SH rank arguments. 
 * 
 * \param useExactSize If it is false, SH3IntegralTable may have larger size than the input three ranks based on pre-computed table. If it is true, the size is determined by the input SH ranks. 
 * */
inline void
InitializeSHTripleIntegrationTable(const int rank0=-1, const int rank1=-1, const int rank2=-1, const bool useExactSize=false)
{
  int dim[3];
  dim[0]=rank0>0? utl::RankToDimSH(rank0) : 0;
  dim[1]=rank1>0? utl::RankToDimSH(rank1) : 0;
  dim[2]=rank2>0? utl::RankToDimSH(rank2) : 0;

  const unsigned* shapeOld = SH3IntegralTable->GetShape();

  // no useExactSize, SH3IntegralTable is larger than what is requested
  if (!useExactSize && SH3IntegralTable->Size()>0 && dim[0]<=shapeOld[0] && dim[1]<=shapeOld[1] && dim[2]<=shapeOld[2])
    return;

  // SH3IntegralTable is the same as what is requested
  if (SH3IntegralTable->Size()>0 && dim[0]==shapeOld[0] && dim[1]==shapeOld[1] && dim[2]==shapeOld[2])
    return; 

  typedef itk::Image<double,3> ImageType;
  static ImageType::Pointer image = ImageType::New();
  if (itk::IsImageEmpty(image))
    itk::ReadImage<ImageType>(utl::SH3Itegralhdr, image);
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType pixelIndex;

  unsigned shape[3];
  for ( int i = 0; i < 3; ++i ) 
    shape[i] = useExactSize? dim[i] : utl::max<int>(dim[i], size[i]);
  
  if (SH3IntegralTable->Size()==0 
    || useExactSize 
    || (!useExactSize && SH3IntegralTable->Size()>0 && (dim[0]>SH3IntegralTable->GetShape()[0] || dim[1]>SH3IntegralTable->GetShape()[1] || dim[2]>SH3IntegralTable->GetShape()[2])) )
    {
    SH3IntegralTable->ReSize(shape);

    unsigned index[3];
    for ( int i = 0; i < shape[0]; ++i ) 
      for ( int j = 0; j < shape[1]; ++j ) 
        for ( int k = 0; k < shape[2]; ++k ) 
          {
          index[0]=i, index[1]=j, index[2]=k;
          if (i<size[0] && j<size[1] && k<size[2] )
            {
            pixelIndex[0]=i, pixelIndex[1]=j, pixelIndex[2]=k;
            (*SH3IntegralTable)(index) = image->GetPixel(pixelIndex);
            }
          else
            {
            std::vector<int> v1, v2, v3;
            v1 = utl::GetIndexSHlm(i);
            v2 = utl::GetIndexSHlm(j);
            v3 = utl::GetIndexSHlm(k);
            (*SH3IntegralTable)(index) = itk::SphericalHarmonicsGenerator<double>::RealTripleIntegration(v1[0],v1[1],v2[0],v2[1],v3[0],v3[1],false);
            }
          }
    }
}

}


#endif 
