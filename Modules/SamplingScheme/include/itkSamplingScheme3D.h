/**
 *       @file  itkSamplingScheme3D.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-20-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingScheme3D_h
#define __itkSamplingScheme3D_h


#include "itkVectorContainer.h"
#include "itkPoint.h"
#include "itkIntTypes.h"

#include "utl.h"


namespace itk
{

/**
 *   \class   SamplingScheme3D
 *   \brief   this class describes sampling in a 3D space (Q space or R space). 
 *
 *   The sampling can be single shell sampling or multiple shell sampling, 
 *   which is determined by m_IndicesInShells. 
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template< class TPixelType=double >
class ITK_EXPORT SamplingScheme3D : public VectorContainer<IdentifierType, Point< TPixelType, 3 > > 
{

public:

  /** Standard class typedefs. */
  typedef SamplingScheme3D             Self;
  typedef VectorContainer<IdentifierType, Point< TPixelType, 3 > >      Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SamplingScheme3D, VectorContainer);

  typedef VectorContainer<IdentifierType, Point< TPixelType, 3 > >    PointsContainer;
  typedef Point<TPixelType, 3>            PointType;
  typedef TPixelType                      ValueType;

  /** Create types that are pointers to each of the container types. */
  typedef typename PointsContainer::Pointer         PointsContainerPointer;
  typedef typename PointsContainer::ConstPointer    PointsContainerConstPointer;

  /** Create types that are iterators for each of the container types. */
  typedef typename PointsContainer::ConstIterator    PointsContainerConstIterator;
  typedef typename PointsContainer::Iterator         PointsContainerIterator;

  typedef utl::NDArray<double,2>                     MatrixType;
  typedef utl_shared_ptr<MatrixType >                MatrixPointer;
  typedef std::vector<double>                        STDVectorType;
  typedef utl_shared_ptr<STDVectorType >             STDVectorPointer;
  typedef std::vector<int>                           IndexVectorType;
  typedef std::vector<IndexVectorType >              Index2DVectorType;
  typedef utl_shared_ptr<Index2DVectorType>          Index2DVectorPointer;
  
  /** Orientations as matrix, return a pointer */
  MatrixPointer GetOrientationsCartesian(const bool alwarysReCalculate=false);
  MatrixPointer GetOrientationsSpherical(const bool alwarysReCalculate=false);
  void SetOrientationsCartesian(const MatrixPointer mat);
  void SetOrientationsSpherical(const MatrixPointer mat);
  MatrixPointer GetOrientationsCartesianInShell(const unsigned int shellIndex) const;
  MatrixPointer GetOrientationsSphericalInShell(const unsigned int shellIndex);
  

  void SetRadiusVector(const STDVectorPointer radiusVec);
  itkGetMacro(RadiusVector, STDVectorPointer);
  STDVectorPointer GetRadiusVectorInShell(unsigned int shellIndex);
  
  itkSetMacro(Tau, double);
  itkGetMacro(Tau, double);
  
  itkSetMacro(DeltaSmall, double);
  itkGetMacro(DeltaSmall, double);
  
  itkSetMacro(DeltaBig, double);
  itkGetMacro(DeltaBig, double);

  itkSetMacro(RadiusThresholdSingleShell, double);
  itkGetMacro(RadiusThresholdSingleShell, double);

  itkSetMacro(IndicesInShells, Index2DVectorPointer);
  itkGetMacro(IndicesInShells, Index2DVectorPointer);
  

  unsigned int GetNumberOfSamples() const
    {
    return this->Size();
    }
  unsigned int GetNumberOfSamplesInShell(const unsigned int shellIndex) const
    {
    if (shellIndex<m_IndicesInShells->size())
      return (*m_IndicesInShells)[shellIndex].size();
    else
      return 0;
    }
  unsigned int GetNumberOfShells() const
    {
    return m_IndicesInShells->size();
    }
  
  IndexVectorType GetNumberOfSamplesAtEachShell() const
    {
    IndexVectorType num;
    if (m_IndicesInShells->size()==0)
      {
      num.push_back(GetNumberOfSamples());
      return num;
      }
    for ( unsigned int i = 0; i < m_IndicesInShells->size(); i += 1 ) 
      num.push_back((*m_IndicesInShells)[i].size());
    return num;
    }
  
  void NormalizeDirections();
  
  void Clear();

  /** remove samples not in m_IndicesInShells  */
  void RemoveSamplesNotIndexed();

  /** generate sampling scheme from random points in sphere  */
  void GenerateFromRandomPoints(const std::vector<int>& numberOfPoints);

  /** Nx3 gradient file, each row is a point in S^2, catesian format  */
  void ReadOrientationFile(const std::string gradFile, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, 
      const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true);
  
  /** Read a list of gradient files. Each file contains gradients for a single shell.   */
  void ReadOrientationFileList(const std::vector<std::string>& gradFileVec, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, 
      const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true);

  std::vector<STDVectorType> GroupRadiusValues();
  
  virtual void CorrectRadiusValues();

  /** push_back a new orientation  */
  void AppendOrientation(const double x, const double y, const double z, const int shell=-1);
  /** push_back a new orientation  */
  void AppendOrientation(const PointType& point, const int shell=-1);
  /** push_back a new orientation and a radius value  */
  void AppendOrientationAndRadiusValue(const double x, const double y, const double z, const double radius, const int shell=-1);
  
  PointType GetOrientation(unsigned int index)
    {
    return this->GetElement(index);
    }
  double GetRadiusValue(unsigned int index)
    {
    return (*m_RadiusVector)[index];
    }
  
  /** return dot products between samples.  */
  MatrixPointer CalculateInnerProductMatrix(const bool isAbsolute=true) const;
  /** return electrostatic energy matrix, the diagonal elements are maximal double value  */
  MatrixPointer CalculateElectrostaticEnergyMatrix(const double order=2.0) const;
  
  /**
   * Packing density: 
   * https://en.wikipedia.org/wiki/Packing_density 
   * https://en.wikipedia.org/wiki/Spherical_cap
   * */
  double CalculatePackingDensity( const bool isSymmetric=true) const;
  double CalculatePackingDensityInShell(const unsigned int shellIndex, const bool isSymmetric=true) const;
  
  /** 
   * Calculate entropy based on spherical cap. 
   * https://en.wikipedia.org/wiki/Spherical_cap  
   * */
  double CalculateSphericalCodeEntropy( const bool isSymmetric=true) const;
  double CalculateSphericalCodeEntropyInShell(const unsigned int shellIndex, const bool isSymmetric=true) const;

  static double CalculateVoronoiEntropy(const MatrixType& grad, const MatrixType& gradTess, const bool isSymmetric=true);
  double CalculateVoronoiEntropy(const int tess=7, const bool isSymmetric=true);
  double CalculateVoronoiEntropyInShell(const unsigned int shellIndex, const int tess=7, const bool isSymmetric=true);
  
  double CalculateMaxDot(const unsigned int index, const bool isSymmetric=true) const;
  double CalculateMaxDotInShell(const unsigned int sampleIndex, const unsigned int shellIndex, const bool isSymmetric=true) const;

  /** calculate the minimum distance for the point with a given index.  */
  double CalculateMinDistance(const unsigned int index, const bool isSymmetric=true) const;
  double CalculateMinDistanceInShell(const unsigned int sampleIndex, const unsigned int shellIndex, const bool isSymmetric=true) const;
  /** calculate the minimum distance for all points.  */
  STDVectorType CalculateMinDistance(const bool isSymmetric=true) const;
  /** calculate the minimum distance for all points in a single shell.  */
  STDVectorType CalculateMinDistanceInShell(const unsigned int shellIndex, const bool isSymmetric=true) const;

  void GetNumbers(int& numberUniqueSamples, int& numberAntipodalSamples, int& numberRepeatedSamples ) const;

  /** calculate electrostatic energy  */
  double CalculateElectrostaticEnergy(const double order=2.0, const bool isNormalize=true, const bool countHalf=true ) const;
  double CalculateElectrostaticEnergyInShell(const unsigned int shellIndex, const double order=2.0, const bool isNormalize=true, const bool countHalf=true ) const;

  /** calculate an upper bound of the minimal distance for a given number of points. 
   *  If isSphericalDistance is true, return spherical distance instead of Euclidean distance. 
   *  http://mathworld.wolfram.com/SphericalCode.html 
   * */
  static double CalculateMinDistanceUpperBound(const unsigned int numberOfPoints, const bool isSphericalDistance=true);

protected:
  SamplingScheme3D();
  virtual ~SamplingScheme3D(){}
  
  typename LightObject::Pointer InternalClone() const;
 
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** tau value  */
  double m_Tau;
  
  /** small delta  */
  double m_DeltaSmall;
  
  /** big delta  */
  double m_DeltaBig;

  /** radius values  */
  STDVectorPointer m_RadiusVector;

  /** 2D vector, each element is a vector containing the indices for orientations in each shell. 
   * If m_IndicesInShells is empty, single shell sampling is used. 
   *
   * \note: There may be some samples whose indices are not in m_IndicesInShells. 
   * These samples which are not indexed can be removed by using RemoveSamplesNotIndexed(). 
   * */
  Index2DVectorPointer m_IndicesInShells;

  /** threshould to separate radius values into different shells.  
   * If it is positive, radius values whose distance is smallter the threshold will be considered in 
   * the same shell, and they will be replaced as their mean b value. 
   * */
  double m_RadiusThresholdSingleShell;
  
  /** Orientation matrix  */
  MatrixPointer m_OrientationsCartesian;
  MatrixPointer m_OrientationsSpherical;
 
private:
  SamplingScheme3D(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  void InsertElement(IdentifierType, PointType);  //purposely not implemented
  void SetElement(IdentifierType, PointType); //purposely not implemented
  PointType & CreateElementAt(IdentifierType); //purposely not implemented
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingScheme3D_hxx)
#include "itkSamplingScheme3D.hxx"
#endif

#endif 
