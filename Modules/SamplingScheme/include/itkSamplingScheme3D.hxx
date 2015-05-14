/**
 *       @file  itkSamplingScheme3D.hxx
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


#ifndef __itkSamplingScheme3D_hxx
#define __itkSamplingScheme3D_hxx


#include "itkSamplingScheme3D.h"
#include "utl.h"

namespace itk
{

template <class TPixelType>
SamplingScheme3D<TPixelType>
::SamplingScheme3D() : 
  m_RadiusVector(new STDVectorType()), 
  m_IndicesInShells(new Index2DVectorType()), 
  m_OrientationsCartesian(new MatrixType()), 
  m_OrientationsSpherical(new MatrixType())
{
  m_Tau = ONE_OVER_4_PI_2; // 4*pi*pi*tau ==1
  m_DeltaSmall = 0;
  m_DeltaBig = m_Tau;  // m_Tau = m_DeltaBig - m_DeltaSmall/3
  m_RadiusThresholdSingleShell = 1e-4;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::SetRadiusVector(const STDVectorPointer radiusVec)
{
  if (m_RadiusVector != radiusVec)
    {
    m_RadiusVector = radiusVec;
    if (m_RadiusVector->size()>0)
      GroupRadiusValues();
    this->Modified();
    }
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::STDVectorPointer
SamplingScheme3D<TPixelType>
::GetRadiusVectorInShell(unsigned int shellIndex)
{
  unsigned int num = GetNumberOfShells();
  STDVectorPointer radiusVec(new STDVectorType());
  if (shellIndex<num && m_RadiusVector->size()>num)
    {
    IndexVectorType vecTmp = (*m_IndicesInShells)[shellIndex];
    for ( unsigned int i = 0; i < vecTmp.size(); i += 1 ) 
      {
      radiusVec->push_back( (*m_RadiusVector)[vecTmp[i] ] );
      }
    }
  return radiusVec;
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::GetOrientationsCartesian()
{
  if (m_OrientationsCartesian->Rows()==GetNumberOfSamples())
    return m_OrientationsCartesian;
  else if (m_OrientationsSpherical->Rows()==GetNumberOfSamples())
    {
    m_OrientationsCartesian = MatrixPointer(new MatrixType());
    *m_OrientationsCartesian = utl::SphericalToCartesian(*m_OrientationsSpherical);
    return m_OrientationsCartesian;
    }
  else if (this->GetNumberOfSamples()>0)
    {
    m_OrientationsCartesian = MatrixPointer(new MatrixType());
    utl::PointsContainerToUtlMatrix<Superclass, double>(*this, *m_OrientationsCartesian);
    return m_OrientationsCartesian;
    }
  else
    utlGlobalException(true, "no orientations");
  return MatrixPointer(new MatrixType());
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::GetOrientationsSpherical()
{
  if (m_OrientationsSpherical->Rows()==GetNumberOfSamples())
    return m_OrientationsSpherical;
  else if (m_OrientationsCartesian->Rows()==GetNumberOfSamples())
    {
    m_OrientationsSpherical = MatrixPointer(new MatrixType());
    *m_OrientationsSpherical = utl::CartesianToSpherical(*m_OrientationsCartesian);
    return m_OrientationsSpherical;
    }
  else if (this->GetNumberOfSamples()>0)
    {
    m_OrientationsCartesian = MatrixPointer(new MatrixType());
    m_OrientationsSpherical = MatrixPointer(new MatrixType());
    utl::PointsContainerToUtlMatrix<Superclass, double>(*this, *m_OrientationsCartesian);
    *m_OrientationsSpherical = utl::CartesianToSpherical(*m_OrientationsCartesian);
    return m_OrientationsSpherical;
    }
  else
    utlGlobalException(true, "no orientations");
  return MatrixPointer(new MatrixType());
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::SetOrientationsCartesian(const MatrixPointer mat)
{
  m_OrientationsCartesian = mat;
  utl::UtlMatrixToPointsContainer<double, Superclass>(*m_OrientationsCartesian, *this);
  m_OrientationsSpherical=MatrixPointer(new MatrixType());
  this->Modified();
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::GetOrientationsSphericalInShell(const unsigned int shellIndex)
{
  unsigned int num = GetNumberOfSamplesInShell(shellIndex);
  utlGlobalException(num==0, "need to set m_IndicesInShells");
  utlException(shellIndex>num, "wrong index");
  MatrixPointer mat (new MatrixType(num, 3));
  m_OrientationsSpherical = GetOrientationsSpherical();
  m_OrientationsSpherical->GetRows((*m_IndicesInShells)[shellIndex], *mat);
  return mat;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::SetOrientationsSpherical(const MatrixPointer mat)
{
  m_OrientationsSpherical = mat;
  m_OrientationsCartesian = MatrixPointer(new MatrixType());
  *m_OrientationsCartesian = utl::SphericalToCartesian(*m_OrientationsSpherical);
  utl::UtlMatrixToPointsContainer<double, Superclass>(*m_OrientationsCartesian, *this);
  this->Modified();
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::GetOrientationsCartesianInShell(const unsigned int shellIndex) const
{
  unsigned int num = GetNumberOfSamplesInShell(shellIndex);
  utlGlobalException(num==0, "need to set m_IndicesInShells");
  utlException(shellIndex>num, "wrong index");
  MatrixPointer mat (new MatrixType(num, 3));
  for ( unsigned int i = 0; i < num; i += 1 ) 
    {
    unsigned int ii = (*m_IndicesInShells)[shellIndex][i];
    (*mat)(i,0) = (*this)[ii][0];
    (*mat)(i,1) = (*this)[ii][1];
    (*mat)(i,2) = (*this)[ii][2];
    }
  return mat;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::NormalizeDirections()
{
  unsigned int num = GetNumberOfSamples();
  double xi,yi,zi, norm2;
  bool isPerformNormalization = false;
  for ( unsigned int i = 0; i < num; i += 1 ) 
    {
    xi = (*this)[i][0];
    yi = (*this)[i][1];
    zi = (*this)[i][2];
    norm2 = xi*xi + yi*yi + zi*zi;
    if (norm2>0 && std::fabs(norm2-1)>1e-8 )
      {
      norm2 = std::sqrt(norm2);
      (*this)[i][0] = xi / norm2;
      (*this)[i][1] = yi / norm2;
      (*this)[i][2] = yi / norm2;
      isPerformNormalization = true;
      }
    }
  if (isPerformNormalization)
    {
    m_OrientationsSpherical=MatrixPointer(new MatrixType());
    m_OrientationsCartesian=MatrixPointer(new MatrixType());
    }
}

template <class TPixelType>
typename LightObject::Pointer
SamplingScheme3D<TPixelType>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  for ( int i = 0; i < this->size(); i += 1 ) 
    rval->push_back((*this)[i]);

  rval->m_DeltaSmall = m_DeltaSmall;
  rval->m_DeltaBig = m_DeltaBig;
  rval->m_Tau = m_Tau;
  rval->m_RadiusThresholdSingleShell = m_RadiusThresholdSingleShell;

  // NOTE: shared_ptr is thread safe, if the data is read only, thus do not need to copy the data block. 
  // NOTE: when changing members in m_SamplingScheme3D, use new to generate another pointer for SmartPointer.
  *rval->m_RadiusVector = *m_RadiusVector;
  *rval->m_IndicesInShells = *m_IndicesInShells;
  rval->m_OrientationsCartesian = m_OrientationsCartesian;
  rval->m_OrientationsSpherical = m_OrientationsSpherical;
  return loPtr;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::Clear()
{
  this->Superclass::Initialize();
  m_RadiusVector=STDVectorPointer(new STDVectorType());
  m_IndicesInShells=Index2DVectorPointer(new Index2DVectorType());
  m_OrientationsSpherical=MatrixPointer(new MatrixType());
  m_OrientationsCartesian=MatrixPointer(new MatrixType());
  this->Modified();
}

template <class TPixelType>
std::vector<typename SamplingScheme3D<TPixelType>::STDVectorType>
SamplingScheme3D<TPixelType>
::GroupRadiusValues()
{
  utlGlobalException(m_RadiusVector->size()==0, "no radius values");
  std::vector<STDVectorType> radiusVectors; 
  if (m_IndicesInShells->size()==0)
    {
    this->m_IndicesInShells = Index2DVectorPointer(new Index2DVectorType()); 
    STDVectorType radiusMax, radiusMin;
    for ( int i = 0; i < m_RadiusVector->size(); i += 1 ) 
      {
      double radius = (*m_RadiusVector)[i];
      int j=0;
      for ( j = 0; j < radiusVectors.size(); j += 1 ) 
        {
        if (radius>=radiusMin[j]-m_RadiusThresholdSingleShell && radius<=radiusMax[j]+m_RadiusThresholdSingleShell)
          {
          radiusVectors[j].push_back(radius);
          (*m_IndicesInShells)[j].push_back(i);
          if (radius<radiusMin[j])
            radiusMin[j] = radius;
          else if (radius>radiusMax[j])
            radiusMax[j] = radius;
          break;
          }
        }
      if (j==radiusVectors.size())
        {
        STDVectorType radiusVecTemp;
        radiusVecTemp.push_back(radius);
        radiusVectors.push_back(radiusVecTemp);
        IndexVectorType radiusIndexTemp;
        radiusIndexTemp.push_back(i);
        m_IndicesInShells->push_back(radiusIndexTemp);
        radiusMin.push_back(radius);
        radiusMax.push_back(radius);
        }
      }

    this->Modified();
    for ( int j = 0; j < radiusVectors.size(); j += 1 ) 
      {
      itkAssertOrThrowMacro(radiusMax[j]-radiusMin[j]<100.0, "the range of radius values is larger than 100, which can not be in the same shell");
      }
    }
  else
    {
    for ( int i = 0; i < m_IndicesInShells->size(); i += 1 ) 
      {
      STDVectorType radiusVecTemp;
      IndexVectorType indexTemp = (*m_IndicesInShells)[i];
      for ( int j = 0; j < indexTemp.size(); j += 1 ) 
        {
        radiusVecTemp.push_back((*m_RadiusVector)[ indexTemp[j] ]);
        }
      radiusVectors.push_back(radiusVecTemp);
      }
    }
  return radiusVectors;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::CorrectRadiusValues()
{
  utlGlobalException(m_RadiusVector->size()==0, "no radius values");
  std::vector<STDVectorType> radiusVectors = GroupRadiusValues(); 

  for ( int j = 0; j < radiusVectors.size(); j += 1 ) 
    {
    double radiusMean=0;
    for ( int k = 0; k < radiusVectors[j].size(); k += 1 ) 
      radiusMean += radiusVectors[j][k];
    radiusMean /= radiusVectors[j].size();
    for ( int k = 0; k < radiusVectors[j].size(); k += 1 ) 
      {
      (*m_RadiusVector)[(*m_IndicesInShells)[j][k] ] = radiusMean;
      }
    }
  this->Modified();
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::ReadOrientationFile(const std::string gradFile, const int NoSymmetricDuple, 
  const int flipx, const int flipy, const int flipz, const bool need_normalize)
{
  Clear();
  m_OrientationsCartesian = utl::ReadGrad<double>(gradFile, NoSymmetricDuple, CARTESIAN_TO_CARTESIAN, flipx, flipy, flipz, need_normalize);
  utl::UtlMatrixToPointsContainer<double, Superclass>(*m_OrientationsCartesian, *this);
  std::vector<int> indicesTmp;
  m_IndicesInShells=Index2DVectorPointer(new Index2DVectorType());
  for ( int j = 0; j < this->size(); j += 1 ) 
    indicesTmp.push_back(j);
  m_IndicesInShells->push_back(indicesTmp);
  this->NormalizeDirections();
  this->Modified();
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::ReadOrientationFileList(const std::vector<std::string>& gradFileVec, const int NoSymmetricDuple, 
  const int flipx, const int flipy, const int flipz, const bool need_normalize)
{
  utlGlobalException(gradFileVec.size()==0, "empty file list");
  Clear();
  MatrixPointer matrix(new MatrixType());

  std::vector<int> indicesTmp;
  m_IndicesInShells=Index2DVectorPointer(new Index2DVectorType());
  int ind=0;
  for ( int i = 0; i < gradFileVec.size(); i += 1 ) 
    {
    matrix = utl::ReadGrad<double>(gradFileVec[i], NoSymmetricDuple, CARTESIAN_TO_CARTESIAN, flipx, flipy, flipz, need_normalize);
    Pointer tmpThis = Self::New(); 
    utl::UtlMatrixToPointsContainer<double, Superclass>(*matrix, *tmpThis);
    indicesTmp.clear();
    for ( int j = 0; j < tmpThis->size(); j += 1 ) 
      {
      this->push_back((*tmpThis)[j] );
      indicesTmp.push_back( ind++ );
      }
    m_IndicesInShells->push_back(indicesTmp);
    }

  this->NormalizeDirections();
  this->Modified();
}

template < class TPixelType >
void
SamplingScheme3D<TPixelType>
::GenerateFromRandomPoints ( const std::vector<int>& numberOfPoints )
{
  Clear();
  int ii=0;
  PointType point;
  for ( int i = 0; i < numberOfPoints.size(); i += 1 ) 
    {
    std::vector<int> indicesTmp;
    for ( int j = 0; j < numberOfPoints[i]; j += 1 ) 
      {
      std::vector<double> pp = utl::RandomPointInSphere(true);
      point[0]=pp[0];  point[1]=pp[1];  point[2]=pp[2];
      this->push_back(point);
      indicesTmp.push_back(ii++);
      }
    m_IndicesInShells->push_back(indicesTmp);
    }
}


template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::AppendOrientation(const double x, const double y, const double z, const int shell)
{
  PointType point;
  point[0]=x;  point[1]=y;  point[2]=z;
  AppendOrientation(point, shell);
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::AppendOrientation(const PointType& point, const int shell)
{
  this->push_back(point);
  m_OrientationsSpherical=MatrixPointer(new MatrixType());
  m_OrientationsCartesian=MatrixPointer(new MatrixType());
  utlException(shell<-1, "shell should be more than -1");
  int num = GetNumberOfSamples();
  if (shell>=0)
    {
    if (m_IndicesInShells->size()>shell)
      (*m_IndicesInShells)[shell].push_back(num);
    if (m_IndicesInShells->size()==shell)
      {
      IndexVectorType indexVec;
      indexVec.push_back(num);
      m_IndicesInShells->push_back(indexVec);
      }
    utlSAException(m_IndicesInShells->size()<shell)(m_IndicesInShells->size())(shell).msg("wrong shell index");
    }
  else
    {
    utlSAException(m_IndicesInShells->size()>1)(m_IndicesInShells->size()).msg("need to set shell index, because there are more than 1 shell");
    if (m_IndicesInShells->size()==1)
      (*m_IndicesInShells)[0].push_back(num);
    }
  this->Modified();
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::AppendOrientationAndRadiusValue(const double x, const double y, const double z, const double radius, const int shell)
{
  PointType point;
  point[0]=x;  point[1]=y;  point[2]=z;
  AppendOrientation(point, shell);
  m_RadiusVector->push_back(radius);
  this->Modified();
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::CalculateInnerProductMatrix(const bool isAbsolute) const
{
  unsigned int num = GetNumberOfSamples();
  MatrixPointer result(new MatrixType(num, num));
  double x,y,z, x1,y1,z1, dot;
  for ( unsigned int i = 0; i < num; i += 1 ) 
    {
    result->operator()(i, i) = 1.0;
    x = this->operator[](i)[0];
    y = this->operator[](i)[1];
    z = this->operator[](i)[2];
    for ( int j = 0; j < i; j += 1 ) 
      {
      x1 = this->operator[](j)[0];
      y1 = this->operator[](j)[1];
      z1 = this->operator[](j)[2];
      dot = x*x1 + y*y1 + z*z1;
      if (dot<0 && isAbsolute)
        result->operator()(i,j) = -dot;
      else
        result->operator()(i,j) = dot;
      result->operator()(j,i) = result->operator()(i,j);
      }
    }
  return result;
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::MatrixPointer
SamplingScheme3D<TPixelType>
::CalculateElectrostaticEnergyMatrix(const double order ) const
{
  double temp, xi, yi, zi, xj, yj, zj;
  unsigned int num = GetNumberOfSamples();
  MatrixPointer result(new MatrixType(num, num));
  bool orderEqual2 = std::abs(order-2)<1e-8 ? true : false;
  for(unsigned int i = 0; i < num; i++) 
    {
    xi = (*this)[i][0];
    yi = (*this)[i][1];
    zi = (*this)[i][2];
    result->operator()(i,i)= std::numeric_limits<double>::max();
    for ( unsigned int j = 0; j< i; j += 1 ) 
      {
      xj = (*this)[j][0];
      yj = (*this)[j][1];
      zj = (*this)[j][2];
      double dot = xi*xj + yi*yj + zi*zj;
      double result1 = 2+2*dot;
      double result2 = 2-2*dot;
      utlWarning(result1<1e-8 || result2<1e-8, "identical directions, d1=("<<xi<<","<<yi<<","<<zi<<") , d2=("<<xj<<","<<yj<<","<<zj<<")");
      temp = (!orderEqual2) ? (1.0/std::pow(result1, order*0.5) + 1.0/std::pow(result2, order*0.5)) : (1.0/result1 + 1.0/result2);  
      result->operator()(i,j) = temp;
      result->operator()(j,i) = temp;
      }
    }
  return result;
}

template <class TPixelType>
double
SamplingScheme3D<TPixelType>
::CalculateMinDistance(const unsigned int index, const bool isSymmetric) const
{
  utlGlobalException(this->GetNumberOfSamples()<2, "No enough points!");
  unsigned int num = GetNumberOfSamples();
  
  double angle, dotMax=-3;  

  double x = (*this)[index][0];
  double y = (*this)[index][1];
  double z = (*this)[index][2];

  for(unsigned int i = 0; i < num; i++) 
    {
    if (i==index)
      continue;

    double xi = (*this)[i][0];
    double yi = (*this)[i][1];
    double zi = (*this)[i][2];
    double dot = x*xi + y*yi + z*zi;
    if (isSymmetric && dot<0 )
      dot = -dot;

    if (dot>dotMax && std::abs(1-dot)>1e-8)
      dotMax = dot;
    }

    if(dotMax >= -1.0 - M_EPS && dotMax <= -1.0 + M_EPS)
      angle = M_PI;
    else if(dotMax <= 1.0 + M_EPS && dotMax >= 1.0 - M_EPS)
      angle = 0;
    else 
      angle = std::acos(dotMax);
  return angle;
}

template <class TPixelType>
double
SamplingScheme3D<TPixelType>
::CalculateMinDistanceInShell(const unsigned int sampleIndex, const unsigned int shellIndex, const bool isSymmetric) const
{
  utlGlobalException(this->GetNumberOfSamples()<2, "No enough points!");
  unsigned int num = GetNumberOfSamplesInShell(shellIndex);
  
  double angle, dotMax=-3;  

  double x = (*this)[sampleIndex][0];
  double y = (*this)[sampleIndex][1];
  double z = (*this)[sampleIndex][2];

  for(unsigned int i = 0; i < num; i++) 
    {
    unsigned int ii = (*m_IndicesInShells)[shellIndex][i];
    if (ii==sampleIndex)
      continue;

    double xi = (*this)[ii][0];
    double yi = (*this)[ii][1];
    double zi = (*this)[ii][2];
    double dot = x*xi + y*yi + z*zi;
    if (isSymmetric && dot<0 )
      dot = -dot;

    if (dot>dotMax && std::abs(1-dot)>1e-8)
      dotMax = dot;
    }

    if(dotMax >= -1.0 - M_EPS && dotMax <= -1.0 + M_EPS)
      angle = M_PI;
    else if(dotMax <= 1.0 + M_EPS && dotMax >= 1.0 - M_EPS)
      angle = 0;
    else 
      angle = std::acos(dotMax);
  return angle;
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::STDVectorType
SamplingScheme3D<TPixelType>
::CalculateMinDistance(const bool isSymmetric) const
{
  STDVectorType angleVec;
  for(unsigned int i = 0; i < GetNumberOfSamples(); i++) 
    {
    double angle = CalculateMinDistance(i, isSymmetric);
    angleVec.push_back(angle);
    }
  return angleVec;
}

template <class TPixelType>
typename SamplingScheme3D<TPixelType>::STDVectorType
SamplingScheme3D<TPixelType>
::CalculateMinDistanceInShell(const unsigned int shellIndex, const bool isSymmetric) const
{
  STDVectorType angleVec;
  unsigned int num = GetNumberOfSamplesInShell(shellIndex);
  for(unsigned int i = 0; i < num; i++) 
    {
    unsigned int ii = (*m_IndicesInShells)[shellIndex][i];
    double angle = CalculateMinDistanceInShell(ii, shellIndex, isSymmetric);
    angleVec.push_back(angle);
    }
  return angleVec;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::GetNumbers(int& numberUniqueSamples, int& numberAntipodalSamples, int& numberRepeatedSamples ) const
{
  double xi, yi, zi, xj, yj, zj;
  unsigned int num = GetNumberOfSamples();
  numberUniqueSamples=0, numberAntipodalSamples=0, numberRepeatedSamples=0;
  for(unsigned int i = 0; i < num; ++i) 
    {
    xi = (*this)[i][0];
    yi = (*this)[i][1];
    zi = (*this)[i][2];
    bool isUnique=true;
    for ( unsigned int j = 0;  j < i; ++j ) 
      {
      xj = (*this)[j][0];
      yj = (*this)[j][1];
      zj = (*this)[j][2];
      double resultSame = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);
      double resultAntipodal = (xi+xj)*(xi+xj) + (yi+yj)*(yi+yj) + (zi+zj)*(zi+zj);
      if (resultSame<1e-8 || resultAntipodal<1e-8)
        isUnique=false;
      if (resultSame<1e-8)
        numberRepeatedSamples++;
      if (resultAntipodal<1e-8)
        numberAntipodalSamples++;
      }
    if (isUnique)
      numberUniqueSamples++;
    }
}

template <class TPixelType>
double
SamplingScheme3D<TPixelType>
::CalculateElectrostaticEnergy(const double order, const bool isNormalize, const bool countHalf ) const
{
  double result = 0;
  double xi, yi, zi, xj, yj, zj;
  int kk=0;
  unsigned int num = GetNumberOfSamples();
  for(unsigned int i = 0; i < num; i++) 
    {
    xi = (*this)[i][0];
    yi = (*this)[i][1];
    zi = (*this)[i][2];
    for ( unsigned int j = 0; (countHalf? (j < i) : (j< num) ); j += 1 ) 
      {
      if (i==j)
        continue;
      xj = (*this)[j][0];
      yj = (*this)[j][1];
      zj = (*this)[j][2];
      double result1 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);
      double result2 = (xi+xj)*(xi+xj) + (yi+yj)*(yi+yj) + (zi+zj)*(zi+zj);
      utlWarning(result1<1e-8 || result2<1e-8, "identical directions, d1=("<<xi<<","<<yi<<","<<zi<<") , d2=("<<xj<<","<<yj<<","<<zj<<")");
      double temp = (order!=2) ? (1.0/std::pow(result1, order*0.5) + 1.0/std::pow(result2, order*0.5)) : (1.0/result1 + 1.0/result2);  
      result += temp; 
      kk++;
      }
    }
  if (isNormalize && kk>0)
    result /= double(kk);
  return result;
}

template <class TPixelType>
double
SamplingScheme3D<TPixelType>
::CalculateElectrostaticEnergyInShell(const unsigned int shellIndex, const double order, const bool isNormalize, const bool countHalf ) const
{
  double result = 0;
  double xi, yi, zi, xj, yj, zj;
  int kk=0;
  unsigned int num = GetNumberOfSamplesInShell(shellIndex);
  for(unsigned int i = 0; i < num; i++) 
    {
    unsigned int ii = (*m_IndicesInShells)[shellIndex][i];
    xi = (*this)[ii][0];
    yi = (*this)[ii][1];
    zi = (*this)[ii][2];
    for ( unsigned int j = 0; (countHalf? (j < i) : (j< num) ); j += 1 ) 
      {
      if (i==j)
        continue;
      unsigned int jj = (*m_IndicesInShells)[shellIndex][j];
      xj = (*this)[jj][0];
      yj = (*this)[jj][1];
      zj = (*this)[jj][2];
      double result1 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);
      double result2 = (xi+xj)*(xi+xj) + (yi+yj)*(yi+yj) + (zi+zj)*(zi+zj);
      utlWarning(result1<1e-8 || result2<1e-8, "identical directions, d1=("<<xi<<","<<yi<<","<<zi<<") , d2=("<<xj<<","<<yj<<","<<zj<<")");
      double temp = (order!=2) ? (1.0/std::pow(result1, order*0.5) + 1.0/std::pow(result2, order*0.5)) : (1.0/result1 + 1.0/result2);  
      result += temp;
      kk++;
      }
    }
  if (isNormalize && kk>0)
    result /= double(kk);
  return result;
}

template <class TPixelType>
double
SamplingScheme3D<TPixelType>
::CalculateMinDistanceUpperBound(const unsigned int numberOfPoints, const bool isSphericalDistance)
{
  double bound = numberOfPoints*M_PI/(6.0*(numberOfPoints-2.0));
  bound = -2.0*std::sin(bound)/(std::cos(2*bound)-1);
  bound = std::sqrt( 4 - bound*bound );
  if (isSphericalDistance)
    bound = std::acos((2.0-bound*bound)/2.0);
  return bound;
}

template <class TPixelType>
void
SamplingScheme3D<TPixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{  
  Superclass::PrintSelf(os, indent);
  PrintVar3(true, m_Tau, m_DeltaSmall, m_DeltaBig, os<<indent);
  MatrixType mat;
  utl::PointsContainerToUtlMatrix<Superclass, double>(*this, mat);
  utl::PrintUtlMatrix(mat, "Orientations = ", " ", os<<indent);
  if (m_OrientationsCartesian->Rows()!=0)
    utl::PrintUtlMatrix(*m_OrientationsCartesian, "m_OrientationsCartesian", " ", os<<indent);
  else if (m_OrientationsSpherical->Rows()!=0)
    utl::PrintUtlMatrix(*m_OrientationsSpherical, "m_OrientationsSpherical", " ", os<<indent);

  utl::PrintVector(*m_RadiusVector, "m_RadiusVector", " ", os<<indent);
  if (GetNumberOfShells()>0)
    {
    os << indent << GetNumberOfShells() << " shells" << std::endl << std::flush;
    for ( int i = 0; i < GetNumberOfShells(); i += 1 ) 
      {
      os << indent << "shell " << i << " : "; 
      utl::PrintVector((*m_IndicesInShells)[i], "indexVector", " ", os<<indent );
      }
    }
  os << indent << "m_RadiusThresholdSingleShell = " <<  m_RadiusThresholdSingleShell << std::endl << std::flush;
}

}


#endif 

