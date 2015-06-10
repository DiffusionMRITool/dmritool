/**
 *       @file  itkPeakContainerHelper.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-26-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkPeakContainerHelper_h
#define __itkPeakContainerHelper_h

#include "utlCore.h"
#include "itkVariableLengthVector.h"


namespace itk
{

typedef enum 
{
  /** peaks peaks stored in x,y,z,value, without the number of peaks  */
  XYZV=0,   
  /** peaks peaks stored in x,y,z, without the number of peaks  */
  XYZ,      
  /** peaks peaks stored in x,y,z,value with the number of peaks in the first value */
  NXYZV,
  /** peaks peaks stored in x,y,z,value with the number of peaks in the first value */
  NXYZ
}PeakType;

/**
 *   \class   PeakContainerHelper
 *   \brief   a VariableLengthVector to contain peaks with different storing types
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
class PeakContainerHelper
{
public:
  
  static std::string GetString(const PeakType peakType)
    {
    if (peakType==XYZV) 
      return std::string("XYZV");
    else if (peakType==XYZ)  
      return std::string("XYZ");
    else if (peakType==NXYZV) 
      return std::string("NXYZV");
    else if (peakType==NXYZ) 
      return std::string("NXYZ");
    else
      utlException(true, "wrong peakType");
    }
  static PeakType GetPeakType(const std::string str)
    {
    if (str=="XYZV") 
      return XYZV;
    else if (str=="XYZ")  
      return XYZ;
    else if (str=="NXYZV") 
      return NXYZV;
    else if (str=="NXYZ") 
      return NXYZ;
    else
      utlException(true, "wrong peakType");
    }
  
  static int GetDimensionPerPeak(const PeakType peakType)
    {
    if (peakType==XYZV || peakType==NXYZV)
      return 4;
    else if (peakType==XYZ || peakType==NXYZ) 
      return 3;
    else
      utlException(true, "wrong peakType");
    }

  static int GetDimension(const PeakType peakType, const int numberOfPeaks)
    {
    int d = GetDimensionPerPeak(peakType);
    if (peakType==XYZ || peakType==XYZV)
      return d*numberOfPeaks;
    else if (peakType==NXYZ || peakType==NXYZV)
      return d*numberOfPeaks+1;
    else
      utlException(true, "wrong peakType");
    }
  
  static int GetNumberOfPeaks(const PeakType peakType, const int dimension)
    {
    int d = GetDimensionPerPeak(peakType);
    int off=0;
    if (peakType==XYZ || peakType==XYZV)
      off = 0;
    else if (peakType==NXYZ || peakType==NXYZV)
      off = 1;
    else
      utlException(true, "wrong peakType");
    utlSAException(!utl::IsInt((1.0*dimension-off)/double(d)))(d)(off)(dimension).msg("wrong size");
    return (dimension-off)/d;
    }
  
  template< class ContainerType >
  static int GetNumberOfPeaks(const PeakType peakType, const int dimension, const ContainerType& vec)
    {
    int nn = GetNumberOfPeaks(peakType, dimension);
    int num=0;
    for ( int i = 0; i < nn; i += 1 ) 
      {
      std::vector<double> p = GetPeak(vec, i, peakType);
      if (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
        num++;
      }
    return num;
    }
  
  template< class ContainerType >
  static std::vector<double> GetPeak(const ContainerType& vec, const int index, const PeakType peakType)
    {
    std::vector<double> peak(3,-1);
    int d = GetDimensionPerPeak(peakType);
    int off=0;
    if (peakType==XYZ || peakType==XYZV)
      off = 0;
    else if (peakType==NXYZ || peakType==NXYZV)
      off = 1;
    else
      utlException(true, "wrong peakType");

    peak[0]=vec[d*index+off+0];
    peak[1]=vec[d*index+off+1];
    peak[2]=vec[d*index+off+2];
    return peak;
    }
  
  template< class VectorType, class ContainerType >
  static void SetPeak(const VectorType& peak, ContainerType& vec, const int index, const PeakType peakType)
    {
    int d = GetDimensionPerPeak(peakType);
    int off=0;
    if (peakType==XYZ || peakType==XYZV)
      off = 0;
    else if (peakType==NXYZ || peakType==NXYZV)
      off = 1;
    else
      utlException(true, "wrong peakType");

    vec[d*index+off+0]=peak[0];
    vec[d*index+off+1]=peak[1];
    vec[d*index+off+2]=peak[2];
    }
  
  template< class ContainerType >
  static double GetPeakValue(const ContainerType& vec, const int index, const PeakType peakType)
    {
    utlException(peakType==XYZ || peakType==NXYZ, "no value for peakType=" << peakType);
    int d = GetDimensionPerPeak(peakType);
    int off=0;
    if (peakType==XYZV)
      off = 0;
    else if (peakType==NXYZV)
      off = 1;
    else
      utlException(true, "wrong peakType");

    return vec[d*index+off+3];
    }
  
  template< class ContainerType >
  static void SetPeakValue(const double peakValue, ContainerType& vec, const int index, const PeakType peakType)
    {
    utlException(peakType==XYZ || peakType==NXYZ, "no value for peakType=" << peakType);
    int d = GetDimensionPerPeak(peakType);
    int off=0;
    if (peakType==XYZV)
      off = 0;
    else if (peakType==NXYZV)
      off = 1;
    else
      utlException(true, "wrong peakType");

    vec[d*index+off+3]=peakValue;
    }


};

}


#endif 

