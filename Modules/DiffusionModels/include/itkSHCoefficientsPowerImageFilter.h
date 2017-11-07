/**
 *       @file  itkSHCoefficientsPowerImageFilter.h
 *      @brief  
 *     Created  "02-04-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkSHCoefficientsPowerImageFilter_h
#define __itkSHCoefficientsPowerImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include <cmath>

#include "utlNDArray.h"
#include "utl.h"
#include "itkSphericalHarmonicsGenerator.h"


namespace itk
{
/** \class SHCoefficientsPowerImageFilter
 *
 * \brief In each vxoel, the input VectorImage is a SH coefficient vector which represents a spherical function. 
 * The output VectorImage is another SH coefficient vector whose spherical function is a power of the input function. 
 * In formula outputFunction = (inputFunction)^k, where k is a given order.
 * When k==2, then we have closed form of the output SH coefficient vector.
 * When k is not 2, we resample the input spherical function, power it, then refit it using SH basis. 
 *
 * \ingroup DiffusionModels
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <class TInput, class TOutput>
class SHCoefficientsPower
{
public:

  typedef utl::NDArray<double, 1>         VectorType;
  typedef utl::NDArray<double, 2>         MatrixType;
  typedef utl_shared_ptr<MatrixType>      MatrixPointer;

  SHCoefficientsPower() :
    m_SHBasisMatrix(new MatrixType())
  {
  m_Power = 1.0;
  m_SHRank = -1;
  }

  ~SHCoefficientsPower()
  {
  }

  bool operator!=(const SHCoefficientsPower & other) const
  {
    return !( *this == other );
  }

  bool operator==(const SHCoefficientsPower & other) const
  {
    return false;
  }

  inline TOutput operator()(const TInput & shCoef)
  {
  if (utl::IsSame<double>(m_Power, 1))
    return shCoef;
  else if (utl::IsSame<double>(m_Power, 2))
    {
    int shInputDim = shCoef.GetSize();
    int shOutputDim = utl::RankToDimSH( 2* utl::DimToRankSH(shInputDim) );

    TOutput out(shOutputDim);
    out.Fill(0.0);
    if (shCoef.GetSquaredNorm()>0)
      {
      unsigned index[3];
      for ( int k = 0; k < shOutputDim; k += 1 ) 
        {
        index[2]=k;
        for ( int i = 0; i < shInputDim; i += 1 ) 
          {
          index[1]=i;
          for ( int j = 0; j < shInputDim; j += 1 ) 
            {
            index[0]=j;
            out[k] += shCoef[i]*shCoef[j]* (*utl::SH3IntegralTable)(index);
            }
          }
        }
      }
    return out;
    }
  else
    {
    int shOutputDim = utl::RankToDimSH( m_SHRank );
    TOutput out(shOutputDim);
    out.Fill(0.0);
    if (shCoef.GetSquaredNorm()>0)
      {
      VectorType shInput(shCoef.GetSize());
      shInput.SetData(const_cast<double* const>(shCoef.GetDataPointer()), shCoef.GetSize());
      VectorType sf = (*m_SHBasisMatrix) * shInput;
      utl::PowerVector(sf.Begin(), sf.End(), m_Power);
      VectorType shOutput = *m_SHBasisMatrixInverse * sf;
      utl::VectorToVector(shOutput, out, shOutput.Size());
      }
    return out;
    }
  }
  
  utlSetGetMacro(Power, double);
  utlSetGetMacro(SHRank, int);
  utlSetGetMacro(SHBasisMatrix, MatrixPointer);
  utlSetGetMacro(SHBasisMatrixInverse, MatrixPointer);

private:

  MatrixPointer m_SHBasisMatrix;
  MatrixPointer m_SHBasisMatrixInverse;

  /** power of spherical function  */
  double m_Power;

  /** it is used when m_Power is not 2 or 1 */
  int m_SHRank;
};
}

template <class TInputImage, class TOutputImage=itk::Image<double,3> >
class ITK_EXPORT SHCoefficientsPowerImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::SHCoefficientsPower<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef SHCoefficientsPowerImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::SHCoefficientsPower<
        typename TInputImage::PixelType,  typename TOutputImage::PixelType> 
      >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SHCoefficientsPowerImageFilter, UnaryFunctorImageFilter );
  
  typedef Functor::SHCoefficientsPower< typename TInputImage::PixelType,
    typename TOutputImage::PixelType>   FunctorType;
  
  typedef typename FunctorType::MatrixType      MatrixType;
  typedef typename FunctorType::MatrixPointer   MatrixPointer;
  
  itkFunctorSetGetMacro(Power, double);
  itkFunctorSetGetMacro(SHRank, int);

protected:
  SHCoefficientsPowerImageFilter()
  {
  }

  virtual ~SHCoefficientsPowerImageFilter()
  {
  }

  void BeforeThreadedGenerateData () ITK_OVERRIDE
    {  
    typename TInputImage::Pointer inputSH = const_cast<TInputImage *>(this->GetInput());
    int shRank = utl::DimToRankSH(inputSH->GetNumberOfComponentsPerPixel());
    if (utl::IsSame<double>(this->GetPower(), 2))
      {
      utl::InitializeSHTripleIntegrationTable(shRank, shRank, 2*shRank);
      }
    else if (!utl::IsSame<double>(this->GetPower(), 1))
      {
      utlGlobalException(this->GetSHRank()<=0, "need to set SHRank to fit the SH coefficients when power is not 1 or 2");
      MatrixPointer grad(new MatrixType()), shMatrixInput(new MatrixType()), shMatrixOutput(new MatrixType()), shMatrixOutputInverse(new MatrixType());
      grad = utl::ReadGrad<double>(4, DIRECTION_NODUPLICATE, CARTESIAN_TO_SPHERICAL);

      shMatrixInput = utl::ComputeSHMatrix(shRank, *grad, SPHERICAL_TO_SPHERICAL);
      this->GetFunctor().SetSHBasisMatrix(shMatrixInput);

      shMatrixOutput = utl::ComputeSHMatrix(this->GetSHRank(), *grad, SPHERICAL_TO_SPHERICAL);
      *shMatrixOutputInverse = utl::PInverseMatrix(*shMatrixOutput);
      this->GetFunctor().SetSHBasisMatrixInverse(shMatrixOutputInverse);
      }
    }

  void GenerateOutputInformation() ITK_OVERRIDE
    {
    Superclass::GenerateOutputInformation();
    typename TOutputImage::Pointer outputPtr = this->GetOutput();
    typename TInputImage::Pointer inputSH = const_cast<TInputImage *>(this->GetInput());

    int numberOfComponentsPerPixel=0; 
    int dimInput = inputSH->GetNumberOfComponentsPerPixel(); 
    if (this->GetPower()==1)
      numberOfComponentsPerPixel = dimInput;
    else if (this->GetPower()==2)
      numberOfComponentsPerPixel = utl::RankToDimSH ( 2*utl::DimToRankSH(dimInput));
    else
      {
      utlGlobalException(this->GetSHRank()<=0, "need to set shRank when power is not 1 or 2");
      numberOfComponentsPerPixel = utl::RankToDimSH (this->GetSHRank() );
      }
    outputPtr->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);
    }
  

private:
  SHCoefficientsPowerImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};
} // end namespace itk


#endif 
