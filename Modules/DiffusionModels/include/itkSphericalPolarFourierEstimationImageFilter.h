/**
 *       @file  itkSphericalPolarFourierEstimationImageFilter.h
 *      @brief  estimate the coeffcients of generalized Spherical Polar Fourier basis 
 *   which can be separated into differe radial basis and sphercial basis.
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-20-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkSphericalPolarFourierEstimationImageFilter_h
#define __itkSphericalPolarFourierEstimationImageFilter_h

#include "itkDiffusionModelEstimationInSphericalCoordinateImageFilter.h"

#include "itkL2RegularizedLeastSquaresSolver.h"
#include "itkL1RegularizedLeastSquaresFISTASolver.h"
#include "itkSpamsWeightedLassoSolver.h"

namespace itk
{

/**
 *   \class   SphericalPolarFourierEstimationImageFilter
 *   \brief   estimate the coeffcients of generalized Spherical Polar Fourier basis 
 *   which can be separated into differe radial basis and sphercial basis.
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT SphericalPolarFourierEstimationImageFilter :
  public DiffusionModelEstimationInSphericalCoordinateImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SphericalPolarFourierEstimationImageFilter         Self;
  typedef DiffusionModelEstimationInSphericalCoordinateImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( SphericalPolarFourierEstimationImageFilter, DiffusionModelEstimationInSphericalCoordinateImageFilter );
  
  /** Convenient Typedefs. */
  itkTypedefMaskedImageToImageMacro(Superclass);

  typedef typename Superclass::MatrixType         MatrixType;
  typedef typename Superclass::VectorType         VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType      STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;

  typedef enum 
    {
    LS=0, 
    L1_2,
    L1_DL
    } EstimationType;
  
  /** Type of L1 solver: 
   * \param FISTA_LS: FISTA with least square initialization
   * \param SPAMS: weighted lasso in spams
   * */
  typedef enum 
    {
    FISTA_LS=0, 
    SPAMS
    } L1SolverType;

  typedef L2RegularizedLeastSquaresSolver<double> L2SolverType;
  typedef L1RegularizedLeastSquaresFISTASolver<double> L1FISTASolverType;
  typedef SpamsWeightedLassoSolver<double> L1SpamsSolverType;

  itkSetMacro(EstimationType, EstimationType);
  itkGetMacro(EstimationType, EstimationType);
  
  itkSetMacro(L1SolverType, L1SolverType);
  itkGetMacro(L1SolverType, L1SolverType);
  
  // itkSetMacro(NeedToConvertToQValue,bool);
  // itkGetMacro(NeedToConvertToQValue,bool);
  
  itkSetMacro(BasisCombinationMatrix, MatrixPointer);
  itkGetMacro(BasisCombinationMatrix, MatrixPointer);
  itkSetMacro(BasisEnergyDL, VectorPointer);
  itkGetMacro(BasisEnergyDL, VectorPointer);
  itkSetMacro(BasisEnergyPowerDL, double);
  itkGetMacro(BasisEnergyPowerDL, double);
  
  // itkSetMacro(RegularizationWeight, VectorPointer);
  itkGetMacro(RegularizationWeight, VectorPointer);
  
  itkSetMacro(IsOriginalBasis,bool);
  itkGetMacro(IsOriginalBasis,bool);
  itkBooleanMacro(IsOriginalBasis);
  
  itkSetMacro(LambdaSpherical, double);
  itkGetMacro(LambdaSpherical, double);
  itkSetMacro(LambdaRadial, double);
  itkGetMacro(LambdaRadial, double);
  itkSetMacro(LambdaL1, double);
  itkGetMacro(LambdaL1, double);


  virtual void SetBasisScale(const double scale);
  itkGetMacro(BasisScale, double);
  
  /** Set/Get the MD image. */
  itkSetObjectMacro(MDImage, ScalarImageType);
  itkGetObjectMacro(MDImage, ScalarImageType);
  itkGetConstObjectMacro(MDImage, ScalarImageType);
  
  /** Set/Get the scale image, which is normally determined by MDImage. */
  itkSetObjectMacro(ScaleImage, ScalarImageType);
  itkGetObjectMacro(ScaleImage, ScalarImageType);
  itkGetConstObjectMacro(ScaleImage, ScalarImageType);

  itkSetGetBooleanMacro(IsAnalyticalB0);
  
  itkSetMacro(B0Weight,double);
  itkGetMacro(B0Weight,double);
  
  itkSetMacro(L1FISTASolver, typename L1FISTASolverType::Pointer);
  itkSetMacro(L1SpamsSolver, typename L1SpamsSolverType::Pointer);

  /** from dimension to rank */
  virtual std::vector<int> DimToRank(const int dimm) const 
    {return std::vector<int>();}
  /** from rank to dimension */
  virtual int RankToDim(const bool is_radial=false, const int radialRank=-1, const int shRank=-1) const  
    {return -1;}

  virtual std::vector<int> GetIndexNLM(const int index) const
    {return std::vector<int>();}; 
  virtual int GetIndexJ(const int n, const int l, const int m) const
    {return -1;};

  /** need to be overidden by subclasses  */
  virtual double ComputeScale(const bool setScale=true);

  bool IsAdaptiveScale() const 
    {
    return !IsImageEmpty(this->m_MDImage);
    }
  
protected:
  SphericalPolarFourierEstimationImageFilter();
  virtual ~SphericalPolarFourierEstimationImageFilter()
    {};
  
  void VerifyInputParameters() const ITK_OVERRIDE;
  
  void InitializeThreadedLibraries() ITK_OVERRIDE;
  
  /** The filter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information accordingly. */
  void GenerateOutputInformation() ITK_OVERRIDE;

  void BeforeThreadedGenerateData () ITK_OVERRIDE;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  

  /** scale for radial basis  */
  double m_BasisScale;
  
  /** regularization parameter in sphercial and radial parts  */
  double m_LambdaSpherical;
  double m_LambdaRadial;
  /** single lambda for L1 norm regularization */
  double m_LambdaL1;
  /** single lambda for L2 norm regularization  */
  double m_LambdaL2;
  
  /** used for learned basis (L1_DL is used)  */
  MatrixPointer  m_BasisCombinationMatrix;
  /** energy of samples in each atom of the learned basis. (L1_DL is used) */
  VectorPointer  m_BasisEnergyDL;
  /**  power of the energy used to calculate m_RegularizationWeight. (L1_DL is used) */
  double m_BasisEnergyPowerDL;

  /** If it is true, use analytical way to ensure E(0)=1  */
  bool m_IsAnalyticalB0;
  /** used for artificial samples for E(0)=1  */
  double m_B0Weight;
  /** basis matrix for artificial b0 shell  */
  MatrixPointer m_BasisMatrixForB0;

  /** used for adaptive scale  */
  typename ScalarImageType::Pointer m_MDImage; 
  typename ScalarImageType::Pointer m_ScaleImage; 

  EstimationType m_EstimationType;

  typename L2SolverType::Pointer   m_L2Solver;
  typename L1FISTASolverType::Pointer   m_L1FISTASolver;
  typename L1SpamsSolverType::Pointer   m_L1SpamsSolver;
  
  L1SolverType   m_L1SolverType;
  
  /** Original SPF basis or dual SPF basis  */
  bool m_IsOriginalBasis;

private:
  SphericalPolarFourierEstimationImageFilter(const Self&);//purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalPolarFourierEstimationImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphericalPolarFourierEstimationImageFilter_hxx)
#include "itkSphericalPolarFourierEstimationImageFilter.hxx"
#endif

#endif 
