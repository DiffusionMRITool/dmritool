/**
 *       @file  itkAddNoiseToDWIImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-15-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkAddNoiseToDWIImageFilter_h
#define __itkAddNoiseToDWIImageFilter_h

#include "itkImageToImageFilter.h"
 
namespace itk
{

/**
 *   \class   AddNoiseToDWIImageFilter
 *   \brief   add noise to DWI data
 *
 *   \ingroup DWIProcessing
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template< typename TInputImage, typename TB0Image=Image<double,3>,
  typename TMaskImage=Image<double,3> >
class AddNoiseToDWIImageFilter
  : public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef AddNoiseToDWIImageFilter                        Self;
  typedef ImageToImageFilter< TInputImage, TInputImage >  Superclass;
  typedef SmartPointer< Self >                            Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(AddNoiseToDWIImageFilter, ImageToImageFilter);

  typedef typename Superclass::InputImageType InputImageType;
//  typedef typename Superclass::InputImageType DWIImageType;
  typedef TB0Image  B0ImageType;
  typedef TMaskImage  MaskImageType;

  typedef typename InputImageType::PixelType PixelType;
  
  typedef enum 
    {
    GAUSSIAN=0,
    RICIAN
    } NoiseType;
 
  
  itkSetMacro(Sigma, double);
  itkGetMacro(Sigma, double);
  itkSetMacro(Noisetype, NoiseType);
  itkGetMacro(Noisetype, NoiseType);
 
  itkSetInputMacro(B0Image, B0ImageType);
  itkSetInputMacro(MaskImage, MaskImageType);

  itkGetInputMacro(B0Image, B0ImageType);
  itkGetInputMacro(MaskImage, MaskImageType);
   
protected:
  AddNoiseToDWIImageFilter();
  virtual ~AddNoiseToDWIImageFilter(){}
 
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
 
  /** Does the real work. */
  virtual void GenerateData() ITK_OVERRIDE;
 
private:
  AddNoiseToDWIImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  
  double m_Sigma;
  NoiseType m_Noisetype;
 
};
} //namespace ITK

// Define instantiation macro for this template.
#define ITK_TEMPLATE_AddNoiseToDWIImageFilter(_, EXPORT, TypeX, TypeY)                  \
  namespace itk                                                                   \
  {                                                                               \
  _( 2 ( class EXPORT AddNoiseToDWIImageFilter< ITK_TEMPLATE_2 TypeX > ) )              \
  namespace Templates                                                             \
  {                                                                               \
  typedef AddNoiseToDWIImageFilter< ITK_TEMPLATE_2 TypeX > AddNoiseToDWIImageFilter##TypeY; \
  }                                                                               \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkAddNoiseToDWIImageFilter+-.h"
#endif
 
#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkAddNoiseToDWIImageFilter_hxx)
#include "itkAddNoiseToDWIImageFilter.hxx"
#endif
 
#endif 
