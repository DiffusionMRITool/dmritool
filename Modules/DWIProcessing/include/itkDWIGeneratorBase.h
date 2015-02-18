/**
 *       @file  itkDWIGeneratorBase.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-22-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDWIGeneratorBase_h
#define __itkDWIGeneratorBase_h

#include "itkPeakContainerHelper.h"
#include "itkImage.h"
#include "itkImageSource.h"

#include <tr1/memory>
#include "itkSamplingSchemeQSpace.h"
#include "itkCylinderModelGenerator.h"

namespace itk
{
  
/** \class DWIGeneratorBase
 *  \brief Generate DWI data based on provided parameter file.
 *
 *  \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template <class TOutputImage, class TScalarImage=Image<float,3> >
class ITK_EXPORT DWIGeneratorBase : public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DWIGeneratorBase            Self;
  typedef ImageSource<TOutputImage>   Superclass;
  typedef SmartPointer<Self>          Pointer;
  
  // [>* Method for creation through the object factory. <]
  // itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DWIGeneratorBase, ImageSource);

  /** Some convenient typedefs. */
  typedef TOutputImage  OutputImageType;
  typedef typename OutputImageType::Pointer               OutputImagePointer;
  typedef typename OutputImageType::SizeType              OutputImageSizeType;
  typedef typename OutputImageType::SpacingType           OutputImageSpacingType;
  typedef typename OutputImageType::IndexType             OutputImageIndexType;
  typedef typename OutputImageType::PointType             OutputImagePointType;
  typedef typename OutputImageType::DirectionType         OutputImageDirectionType;
  typedef typename OutputImageType::RegionType            OutputImageRegionType;
  typedef typename OutputImageType::PixelType             OutputImagePixelType;
  typedef typename OutputImageType::InternalPixelType     OutputImageInternalPixelType;

  /** Output image dimension */
  itkStaticConstMacro(OutputImageDimension, unsigned int, OutputImageType::ImageDimension);

  /** B0 Image */
//  typedef Image<float, OutputImageDimension>     ScalarImageType;
  typedef TScalarImage                               ScalarImageType;
  typedef typename ScalarImageType::Pointer          ScalarImagePointer;
    
  /** Orientation Matrice Type */
  typedef double                                           PrecisionType;
  typedef utl::NDArray<PrecisionType,2>                    MatrixType;
  typedef utl::NDArray<PrecisionType,1>                    VectorType;
  typedef utl_shared_ptr<MatrixType >                      MatrixPointer;
  typedef std::vector<PrecisionType>                       STDVectorType;
  typedef utl_shared_ptr<STDVectorType>                    STDVectorPointer;
  
  /** Some convenient typedefs for diffusion parameters. */  
  typedef typename std::vector<double>                       DiffusionParameterValuesType;
  typedef typename std::vector<DiffusionParameterValuesType> DiffusionParameterContainerType;
  
  typedef CylinderModelGenerator<double>                     CylinderModelType;
  typedef typename CylinderModelType::Pointer                CylinderModelPointer;
  
  typedef SamplingSchemeQSpace<double>                       SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer         SamplingSchemeQSpacePointer;
  
  typedef SamplingScheme3D<double>                           SamplingSchemeRSpaceType;
  typedef typename SamplingSchemeRSpaceType::Pointer         SamplingSchemeRSpacePointer;
  
  /** \param TENSOR_SYM_CARTESIAN symmetric tensor with the principal direction (x,y,z) in Cartesian coordinates
   *  \param TENSOR_SYM_SPHERICAL symmetric tensor with the principal direction (theta,phi) in spherical coordinates
   *  \param TENSOR general tensor with three Eular angles
   *  \param CYLINDER_SPHERICAL cylinder spherical model with two spherical angles (theta,phi)
   **/
  typedef enum 
    {
    SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS=0,
    SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS,
    TENSOR_IN_EULER_ANGLES,
    CYLINDER_SPHERICAL_MODEL
    } ModelType;
  
  itkSetEnumMacro( ModelType, ModelType );
  itkGetEnumMacro( ModelType, ModelType );
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);

  itkSetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);
  itkGetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);
  
  itkSetObjectMacro(CylinderModel, CylinderModelType);
  itkGetObjectMacro(CylinderModel, CylinderModelType);

  itkSetMacro(PeakType, PeakType);
  itkGetMacro(PeakType, PeakType);
  
  itkSetMacro( B0Scale, double );
  itkGetMacro( B0Scale, double );
  
  /** Image size with the same diffusion parameters in all voxels  */
  itkSetMacro( OutputSize, OutputImageSizeType );
  itkGetMacro( OutputSize, OutputImageSizeType );
    
  unsigned int
  GetNumberOfQSpaceSamples() const
  {
    return m_SamplingSchemeQSpace->GetNumberOfSamples();
  }
  unsigned int
  GetNumberOfRSpaceSamples() const
  {
    return m_SamplingSchemeRSpace->GetNumberOfSamples();
  }
  
  /** Parameters */
  itkSetMacro(NoiseSigma, double);
  itkGetMacro(NoiseSigma, double);
  itkSetMacro(SNR, double);
  itkGetMacro(SNR, double);
  itkSetMacro(ODFOrder, unsigned int);
  itkGetMacro(ODFOrder, unsigned int);
  
  itkSetMacro(MaxNumberOfPeaks, int);
  itkGetMacro(MaxNumberOfPeaks, int);
  
  itkSetGetBooleanMacro(IsOutputDWI);
  itkSetGetBooleanMacro(IsOutputEAP);
  itkSetGetBooleanMacro(IsOutputODF);
  itkSetGetBooleanMacro(IsOutputRTO);
  itkSetGetBooleanMacro(IsOutputMSD);
  
  /** Get B0Image */
  ScalarImageType* GetB0Image()
  {
    ScalarImageType* ptr = static_cast<ScalarImageType*>(this->GetOutputs()[1].GetPointer());
    return ptr;    
  }
  
  OutputImageType* GetDWIImage()
  {
    OutputImageType* ptr = static_cast<OutputImageType*>(this->GetOutputs()[0].GetPointer());
    return ptr;    
  }
  
  OutputImageType* GetODFImage()
  {
    OutputImageType* ptr = static_cast<OutputImageType*>(this->GetOutputs()[2].GetPointer());
    return ptr;    
  }
  
  OutputImageType* GetEAPImage()
  {
    OutputImageType* ptr = static_cast<OutputImageType*>(this->GetOutputs()[3].GetPointer());
    return ptr;    
  }
  
  OutputImageType* GetPeakImage()
  {
    OutputImageType* ptr = static_cast<OutputImageType*>(this->GetOutputs()[4].GetPointer());
    return ptr;    
  }
  
  ScalarImageType* GetRTOImage()
  {
    ScalarImageType* ptr = static_cast<ScalarImageType*>(this->GetOutputs()[5].GetPointer());
    return ptr;    
  }
  ScalarImageType* GetMSDImage()
  {
    ScalarImageType* ptr = static_cast<ScalarImageType*>(this->GetOutputs()[6].GetPointer());
    return ptr;    
  }

protected:
  DWIGeneratorBase();
  ~DWIGeneratorBase();

  typename LightObject::Pointer InternalClone() const;
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData()=0;
  
  /** initialization, test  */
  virtual void Initialization();
  /** allocate all outputs  */
  virtual void AllocateOutputs();
  
  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;

  SamplingSchemeRSpacePointer m_SamplingSchemeRSpace;

  bool m_IsOutputDWI;
  bool m_IsOutputODF;
  bool m_IsOutputEAP;
  bool m_IsOutputRTO;
  bool m_IsOutputMSD;
  
  double m_NoiseSigma;
  double m_SNR;

  unsigned int m_ODFOrder;
  
  ModelType m_ModelType;
  double m_B0Scale;
  
  int m_MaxNumberOfPeaks;
  PeakType m_PeakType;

  CylinderModelPointer m_CylinderModel;  
  
  OutputImageSizeType m_OutputSize;
  OutputImageSpacingType m_OutputSpacing;
  OutputImagePointType m_OutputOrigin;
  OutputImageDirectionType m_OutputDirection;

private:
  DWIGeneratorBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
};


} //namespace ITK


#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDWIGeneratorBase_hxx)
#include "itkDWIGeneratorBase.hxx"
#endif

#endif 
