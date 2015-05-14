/**
 *       @file  utlCore.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "09-20-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlCore_h
#define __utlCore_h


#include "utlCoreMacro.h"
#include "utlSmartAssert.h"

// Include OS-specific headers.
#if UTL_OS==1
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#elif UTL_OS==2
#include <windows.h>
#endif

#include <iostream>
#include <fstream>

#include <algorithm>
#include <utility>
#include <vector>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cassert>
#include <complex>
#include <cmath>
#include <limits>
#include <iomanip>
#include <ctime>
// #include <wordexp.h>

#include "utlCoreMKL.h"

#include "utlCommandLineParser.h"

namespace utl 
{

//! Return the value of a system timer, with a millisecond precision.
/**
\note The timer does not necessarily starts from \c 0.
 **/
inline unsigned long 
Time() 
{
#if UTL_OS==1
  struct timeval st_time;
  gettimeofday(&st_time,0);
  return (unsigned long)(st_time.tv_usec/1000 + st_time.tv_sec*1000);
#elif UTL_OS==2
  SYSTEMTIME st_time;
  GetSystemTime(&st_time);
  return (unsigned long)(st_time.wMilliseconds + 1000*(st_time.wSecond + 60*(st_time.wMinute + 60*st_time.wHour)));
#else
  return 0;
#endif
}

/** Implement a tic/toc mechanism to display elapsed time of algorithms, modified from CImg  */
inline unsigned long 
TicToc(const bool is_tic, std::ostream& os=std::cout) 
{
  static unsigned long t0 = 0;
  const unsigned long t = Time();
  if (is_tic) return (t0 = t);
  const unsigned long dt = t>=t0?(t - t0): std::numeric_limits<unsigned long>::max();
  const unsigned int
    edays = (unsigned int)(dt/86400000.0),
  ehours = (unsigned int)((dt - edays*86400000.0)/3600000.0),
  emin = (unsigned int)((dt - edays*86400000.0 - ehours*3600000.0)/60000.0),
  esec = (unsigned int)((dt - edays*86400000.0 - ehours*3600000.0 - emin*60000.0)/1000.0),
  ems = (unsigned int)(dt - edays*86400000.0 - ehours*3600000.0 - emin*60000.0 - esec*1000.0);
  if (!edays && !ehours && !emin && !esec)
    os << "[UTL] Elapsed time: " << ems << " ms" <<  std::endl << std::flush;
  else {
    if (!edays && !ehours && !emin)
      os << "[UTL] Elapsed time: " << esec << " sec "<< ems << " ms" <<  std::endl << std::flush;
    else {
      if (!edays && !ehours)
        os << "[UTL] Elapsed time: " << emin << " min " << esec << " sec "<< ems << " ms" <<  std::endl << std::flush;
      else{
        if (!edays)
          os << "[UTL] Elapsed time: " << ehours << " hours "<< emin << " min " << esec << " sec "<< ems << "ms" <<  std::endl << std::flush;
        else{
          os << "[UTL] Elapsed time: " << edays << " dats " << ehours << " hours "<< emin << " min " << esec << " sec "<< ems << "ms" <<  std::endl << std::flush;
        }
      }
    }
  }
  return t;
}

//! Start tic/toc timer for time measurement between code instructions.
/**
  \return Current value of the timer (same value as time()).
 **/
inline unsigned long 
Tic(std::ostream& os=std::cout) 
{
  return TicToc(true, os);
}

//! End tic/toc timer and displays elapsed time from last call to tic().
/**
  \return Time elapsed (in ms) since last call to tic().
 **/
inline unsigned long 
Toc(std::ostream& os=std::cout) 
{
  return TicToc(false, os);
}



/** SwapBytes(ptr,sizeof(short),8)  */
inline void 
SwapBytes(void *ptr, const int sizePerElement, const int count) 
{
  char minibuf[32], *buffer=(char *)ptr;
  for (int y=0; y<count; y++) 
    {
    for (int x=0; x<sizePerElement; x++) 
      minibuf[x] = *(buffer++);
    for (int x=0; x<sizePerElement; x++) 
      *(--buffer) = minibuf[x];
    buffer+=sizePerElement;
    }
}


template <typename T>
inline T median(std::vector<T> values) 
{
  // values does not change after median
  const unsigned int N = values.size();
  std::nth_element(values.begin(),values.begin()+N/2,values.end());
  return values[N/2];
}

template <class T>
inline T min(const std::vector<T> v)
{
  return *std::min_element(v.begin(), v.end());
}

template <class T>
inline T max(const std::vector<T> v)
{
  return *std::max_element(v.begin(), v.end());
}


template <class Iterator>
inline unsigned int argmin(Iterator i1, Iterator i2)
{
  Iterator iter = std::min_element(i1, i2);
  int n=0;
  for ( Iterator ii = i1; ii!=iter; ++ii ) 
    n++;
  return n;
}

template <class Iterator>
inline unsigned int argmax(Iterator i1, Iterator i2)
{
  Iterator iter = std::max_element(i1, i2);
  int n=0;
  for ( Iterator ii = i1; ii!=iter;  ++ii ) 
    n++;
  return n;
}

//! Return the minimum between \p a and \p b.
template<typename T> inline const T& min(const T& a,const T& b) { return a<=b?a:b; }
//! Return the minimum between \p a,\p b and \a c.
template<typename T> inline const T& min(const T& a,const T& b,const T& c) { return min(min(a,b),c); }
//! Return the minimum between \p a,\p b,\p c and \p d.
template<typename T> inline const T& min(const T& a,const T& b,const T& c,const T& d) { return min(min(min(a,b),c),d); }
//! Return the maximum between \p a and \p b.
template<typename T> inline const T& max(const T& a,const T& b) { return a>=b?a:b; }
//! Return the maximum between \p a,\p b and \p c.
template<typename T> inline const T& max(const T& a,const T& b,const T& c) { return max(max(a,b),c); }
//! Return the maximum between \p a,\p b,\p c and \p d.
template<typename T> inline const T& max(const T& a,const T& b,const T& c,const T& d) { return max(max(a,b,c),d); }
//! Return the sign of \p x.
template<typename T> inline char sign(const T& x) { return (x<0)?-1:(x==0?0:1); }

inline std::string 
StringToLowerCase(const std::string str)
{
  std::string result(str);
  for (std::string::iterator c = result.begin(); c != result.end(); ++c) 
    *c = tolower(*c);
  return result;
}

inline std::string 
StringToUpperCase(const std::string str)
{
  std::string result(str);
  for (std::string::iterator c = result.begin(); c != result.end(); ++c) 
    *c = toupper(*c);
  return result;
}

/** compare two std::string objects, ignore case */
inline bool 
StringCompareCaseIgnored(const std::string str1, const std::string str2) 
{
  return StringToLowerCase(str1)==StringToLowerCase(str2);
}

// replace with as many times as it shows up in source.
// write the result into source.
inline void 
ReplaceString(std::string& source, const char* replace, const char* with)
{
  const char *src = source.c_str();
  char *searchPos = const_cast<char *>(strstr(src,replace));

  // get out quick if string is not found
  if (!searchPos)
    {
    return;
    }

  // perform replacements until done
  size_t replaceSize = strlen(replace);
  // do while hangs if replaceSize is 0
  if(replaceSize == 0)
    {
    return;
    }
  char *orig = strdup(src);
  char *currentPos = orig;
  searchPos = searchPos - src + orig;

  // initialize the result
  source.erase(source.begin(),source.end());
  do
    {
    *searchPos = '\0';
    source += currentPos;
    currentPos = searchPos + replaceSize;
    // replace
    source += with;
    searchPos = strstr(currentPos,replace);
    }
  while (searchPos);

  // copy any trailing text
  source += currentPos;
  free(orig);
}

/** remove double slashes not at the start  */
inline std::string 
ConvertToWindowsOutputPath(const char* path)
{
  std::string ret;
  // make it big enough for all of path and double quotes
  ret.reserve(strlen(path)+3);
  // put path into the string
  ret.assign(path);
  ret = path;
  std::string::size_type pos = 0;
  // first convert all of the slashes
  while((pos = ret.find('/', pos)) != std::string::npos)
    {
    ret[pos] = '\\';
    pos++;
    }
  // check for really small paths
  if(ret.size() < 2)
    {
    return ret;
    }
  // now clean up a bit and remove double slashes
  // Only if it is not the first position in the path which is a network
  // path on windows
  pos = 1; // start at position 1
  if(ret[0] == '\"')
    {
    pos = 2;  // if the string is already quoted then start at 2
    if(ret.size() < 3)
      {
      return ret;
      }
    }
  while((pos = ret.find("\\\\", pos)) != std::string::npos)
    {
    ret.erase(pos, 1);
    }
  // now double quote the path if it has spaces in it
  // and is not already double quoted
  if(ret.find(' ') != std::string::npos
     && ret[0] != '\"')
    {
    ret.insert(static_cast<std::string::size_type>(0),
               static_cast<std::string::size_type>(1), '\"');
    ret.append(1, '\"');
    }
  return ret;
}

/** change // to /, and escape any spaces in the path  */
inline std::string 
ConvertToUnixOutputPath(const char* path)
{
  std::string ret = path;

  // remove // except at the beginning might be a cygwin drive
  std::string::size_type pos=1;
  while((pos = ret.find("//", pos)) != std::string::npos)
    {
    ret.erase(pos, 1);
    }
  // escape spaces and () in the path
  if(ret.find_first_of(" ") != std::string::npos)
    {
    std::string result = "";
    char lastch = 1;
    for(const char* ch = ret.c_str(); *ch != '\0'; ++ch)
      {
        // if it is already escaped then don't try to escape it again
      if((*ch == ' ') && lastch != '\\')
        {
        result += '\\';
        }
      result += *ch;
      lastch = *ch;
      }
    ret = result;
    }
  return ret;
}

/**  convert windows slashes to unix slashes */ 
inline void ConvertToUnixSlashes(std::string& path)
{
  const char* pathCString = path.c_str();
  bool hasDoubleSlash = false;

  const char* pos0 = pathCString;
  const char* pos1 = pathCString+1;
  for (std::string::size_type pos = 0; *pos0; ++ pos )
    {
    // make sure we don't convert an escaped space to a unix slash
    if ( *pos0 == '\\' && *pos1 != ' ' )
      {
      path[pos] = '/';
      }

    // Also, reuse the loop to check for slash followed by another slash
    if (*pos1 == '/' && *(pos1+1) == '/' && !hasDoubleSlash)
      {
#ifdef _WIN32
      // However, on windows if the first characters are both slashes,
      // then keep them that way, so that network paths can be handled.
      if ( pos > 0)
        {
        hasDoubleSlash = true;
        }
#else
      hasDoubleSlash = true;
#endif
      }

    pos0 ++;
    pos1 ++;
    }

  if ( hasDoubleSlash )
    {
    ReplaceString(path, "//", "/");
    }

  // remove any trailing slash
  if(!path.empty())
    {
    // if there is a tilda ~ then replace it with HOME
    pathCString = path.c_str();
    if(pathCString[0] == '~' && (pathCString[1] == '/' || pathCString[1] == '\0'))
      {
      const char* homeEnv = getenv("HOME");
      if (homeEnv)
        {
        path.replace(0,1,homeEnv);
        }
      }
    // remove trailing slash if the path is more than
    // a single /
    pathCString = path.c_str();
    if(path.size() > 1 && *(pathCString+(path.size()-1)) == '/')
      {
      // if it is c:/ then do not remove the trailing slash
      if(!((path.size() == 3 && pathCString[1] == ':')))
        {
        path = path.substr(0, path.size()-1);
        }
      }
    }
}

inline std::string 
CreateExpandedPath(const std::string & path)
{  
  // std::string result;
  // wordexp_t exp_result;
  // wordexp(path.c_str(), &exp_result, 0);
  // result.assign(exp_result.we_wordv[0]);
  // wordfree(&exp_result);

  std::string result(path);
  ConvertToUnixSlashes(result);

#ifdef _WIN32
  return ConvertToWindowsOutputPath(result.c_str());
#else
  return ConvertToUnixOutputPath(result.c_str());
#endif
}


inline bool 
IsFileExist ( const std::string file )
{
  std::fstream ii;
  ii.open(file.c_str(), std::fstream::in);
  if (ii)
    {
    ii.close();
    return true;
    }
  else
    return false;
}

/** path with the last "/" 
*  http://www.cplusplus.com/reference/string/string/find_last_of/ 
*  GetPath("/home/my.cpp", path, file) will get path=="/home/", file=="my.cpp"
*  */
inline void 
GetPath(const std::string fileNameAbsolute, std::string& path, std::string& file)
{
  size_t found = fileNameAbsolute.find_last_of("/\\");
  path = fileNameAbsolute.substr( 0, found + 1);
  file = fileNameAbsolute.substr( found + 1);
}

/** GetFileExtension("/home/dwi.hdr", ext, file) will get ext=="hdr" and file=="/home/dwi"  */
inline void 
GetFileExtension(const std::string fileNameAbsolute, std::string& ext, std::string& fileNoExt)
{
  size_t found = fileNameAbsolute.find_last_of(".");
  ext = fileNameAbsolute.substr(found + 1); 
  fileNoExt = fileNameAbsolute.substr(0,found); 
  // utlAssert(ext!=fileNoExt, "Logical Error! There is no ext! ext="<<ext <<", fileNoExt="<<fileNoExt);
  if (ext==fileNoExt)
    {
    // no ".", no extension, ext is empty
    ext = "";
    }
  else if (ext=="gz")
    {
    size_t found2 = fileNoExt.find_last_of(".");
    std::string ext2 = fileNoExt.substr(found2 + 1);
    if (ext2=="nii")
      {
      ext = "nii.gz";
      fileNoExt = fileNoExt.substr(0,found2);
      }
    }
}

/** http://www.cplusplus.com/reference/iostream/manipulators/setfill/ 
 * Convert a 'number' into a zero padded string.
 * ZeroPad(5, 4); produces "0005" */
inline std::string 
ZeroPad(const unsigned int number, const unsigned int paddedLength)
{
  std::stringstream padded;
  padded << std::setfill('0') << std::setw(paddedLength) << number;
  return padded.str();
}

/** getSequentialFileName("dwi", 3, "hdr") will generate dwi_000003.hdr */
inline std::string 
GetSequentialFileName(const std::string filePrefix, const unsigned int iteration, const std::string fileExtension, const unsigned int paddedLength=6)
{
  std::stringstream padded;
  padded << filePrefix << "_" << ZeroPad(iteration, paddedLength) << "." << fileExtension;
  return padded.str();
}

/** Normalize values into [0,1] using the minimal value and the maximal value  */
template<class VectorType>
inline VectorType
NormalizeMinMax(const VectorType & v, const int nSize) 
{
  VectorType norm(v);
  double max = -std::numeric_limits<double>::max();
  double min = std::numeric_limits<double>::max();

  for(unsigned int i = 0; i < nSize; i++) 
    {
    if(v[i] > max)
      max = v[i];
    if(v[i] < min)
      min = v[i];
    }

  for(unsigned int i = 0; i < v.size(); i++) 
    {
    if((max-min) != 0)
      norm[i] = ((v[i]-min)/(max-min));
    else 
      norm[i] = 0;
    }

  return norm;
}

template<typename T>
inline std::vector<T> 
NormalizeMinMax(std::vector<T> & v) 
{
  return NormalizeMinMax<std::vector<T> >(v, v.size());
}

/** normalize values using L2 norm   */
template<class VectorType>
inline VectorType
NormalizeUnitNorm(const VectorType & v, const int nSize) 
{
  VectorType v1(v);
  double  norm = 0.0;
  for(unsigned int i = 0; i < nSize; i++) 
    norm += v1[i]*v1[i];
  if (norm>1e-20)
    {
    double normsqrt = std::sqrt(norm);
    for ( int i = 0; i < v.size(); i += 1 ) 
      v1[i] /= normsqrt;
    }
  return v1;
}

template<typename T>
inline std::vector<T> 
NormalizeUnitNorm(const std::vector<T> & v) 
{
  return NormalizeUnitNorm<std::vector<T> >(v, v.size());
}

/** normalize values using the maximal value  */
template<class VectorType>
inline VectorType
NormalizeMax(const VectorType & v, const int nSize) 
{
  VectorType norm(v);
  double max = -std::numeric_limits<double>::max();

  for(unsigned int i = 0; i < nSize; i++) 
    {
    if(v[i] > max)
      max = v[i];
    }

  for(unsigned int i = 0; i < nSize; i++) 
    {
    if( max != 0)
      norm[i] = v[i]/max;
    }

  return norm;
}

template<typename T>
inline std::vector<T> 
NormalizeMax(const std::vector<T> & v) 
{
  return NormalizeMax<std::vector<T> >(v, v.size());
}

/** test if a string means a number   */
inline bool 
IsNumber( const std::string input)
{
  if (input.size()==0) 
    return false;
  char *pend;
  double number = strtod(input.c_str(),&pend);
  if (*pend!='\0') 
    return false;
  return true; 
} 

/** convert number to string  */
template <class T>
inline std::string
ConvertNumberToString ( const T value )
{
  std::ostringstream strs;
  strs << value;
  return strs.str();
}

/** convert string to number  */
template <class T>
inline T
ConvertStringToNumber ( const std::string input )
{
  // // The first way also works
  // utlException(input.size()==0, input << " is null");
  // char *pend;
  // T result = strtod(input.c_str(),&pend);
  // utlException(*pend!='\0', input << " is not a number");
  // return result;

  // the second way
  std::istringstream ss(input);
  T result;
  ss >> result;
  if (!ss)
    utlException(true, "input="<<input);
  return result;
}

template <class T>
inline int 
RoundNumber ( const T x )
{
  return x>=0 ? std::floor(x+0.5) : std::ceil(x-0.5);
}


inline bool 
IsInt( const std::string input, const double epss=1e-10)
{
  if (input.size()==0) 
    return false;
  char *pend;
  double dd = strtod(input.c_str(),&pend);
  if (*pend!='\0') 
    return false;
  double diff = dd - (long)(dd);
  return diff<epss && diff>-epss;
} 

inline bool 
IsInt( const double dd, const double epss=1e-10)
{
  double diff = dd - (long)(dd);
  return diff<epss && diff>-epss;
} 

inline bool 
IsInt( const float dd, const double epss=1e-10)
{
  double diff = dd - (long)(dd);
  return diff<epss && diff>-epss;
} 

inline bool 
IsEven(const int value)
{
  return (value % 2 == 0);
}

inline
bool IsOdd(const int value)
{
  // NOTE: use value%2==1 is wrong, because it can be -1.
  return !IsEven(value);
}

template <class VectorType>
inline bool 
IsContainsNaN(const VectorType& a, const int size)
{
  for(unsigned int i = 0; i < size; ++i)
  {
    if(std::isnan(a[i]))
      return true;
  }
  return false;
}

/** separate a string into a string vector based on given delimit  */
inline void
SplitString (const std::string str, std::vector<std::string>& strVec, const char* cc=" ")
{
  strVec.clear();
  char str2[str.size()+1];
  for ( int i = 0; i < str.size(); i += 1 ) 
    str2[i] = str[i];
  str2[str.size()]='\0';
  char *pch = strtok(str2, cc);
  while (pch)
    {
    strVec.push_back(std::string(pch));
    pch = strtok(NULL, cc);
    }
}

/** Read a file, put all strings into a string 2D vector, i.e. vector<vector<string> >  */
inline void
ReadLines ( const std::string filename, std::vector < std::vector<std::string> >& strVec, const char* cc=" ")
{
  strVec.clear();
  std::string line;

  std::ifstream infile(filename.c_str(), std::ios_base::in);
  utlGlobalAssert(infile, "failed to open input file " << filename);
  while (std::getline(infile, line, '\n'))
    {
    std::vector<std::string> tmp;
    SplitString(line, tmp, cc);
    // char *pch = strtok((char *)line.c_str(), cc);
    // while (pch)
    //   {
    //   tmp.push_back(std::string(pch));
    //   pch = strtok(NULL, cc);
    //   }
    // do not add blank lines
    if (tmp.size()>0)
      strVec.push_back(tmp);
    }

  // std::cout << "read " << strVec.size() << " lines from " << filename << ".\n";
  infile.close();
}

/** same with ReadLines(), ignore the first line if it is shows the number of rows.  */
inline void 
ReadLinesFirstlineCheck ( const std::string filename, std::vector < std::vector<std::string> >& strVec, const char* cc=" " )
{
  ReadLines(filename, strVec, cc);
  if (strVec.size()>0 && strVec[0].size()==1 && IsInt(strVec[0][0]))
    {
    int nSize;
    std::istringstream ( strVec[0][0] ) >> nSize;
    if (nSize+1==strVec.size())
      {
      std::vector<std::vector<std::string> > strVec_back(strVec);
      strVec.clear();
      for ( int i = 1; i < strVec_back.size(); i += 1 ) 
        strVec.push_back(strVec_back[i]);
      }
    }
}

template <class TVectorType>
inline double
GetSumOfVector( const TVectorType& vec, const int NSize)
{
  double sumTemp=0;
  for ( int i = 0; i < NSize; i += 1 ) 
    sumTemp += vec[i];
  return sumTemp;
}

template <typename T>
std::ostream& 
operator<<(std::ostream& out, const std::vector<T>& vec )
{
  if (vec.size()>0)
    {
    for ( int i = 0; i < vec.size()-1; i += 1 ) 
      {
      out << vec[i] << " ";
      }
    out << vec[vec.size()-1];
    }
  else
    out << "(empty vector)";
  return out;
}

/** Get statistics (min, max, mean, std) from a given range of values.  */
template <class IteratorType>
inline std::vector<double> 
GetContainerStats ( IteratorType v1, IteratorType v2 )
{
  if (v1==v2)
    return std::vector<double>(4,0);

  double m = *v1, M = m;
  double S = 0, S2 = 0;
  int n=0;
  for ( IteratorType ptr=v1; ptr!=v2;  ++ptr ) 
    {
    const double val = *ptr;
    if (val<m) 
      m = val; 
    if (val>M) 
      M = val;
    S += val;
    S2 += val*val;
    n++;
    }
  double mean = S/(double)n;
  double variance = n>1?(S2 - S*S/n)/(n-1):0;
  variance = std::sqrt(variance);
  
  std::vector<double> result(4);
  result[0]=m, result[1]=M, result[2]=mean, result[3]=variance;
  return result;
}

template <class IteratorType>
int
GetNumberOfNonZeroValues ( IteratorType v, IteratorType v2, const double threshold=1e-6 )
{
  if (v==v2)
    return 0;

  int nnz = 0;
  for ( IteratorType ptr=v; ptr != v2;  ++ptr ) 
    {
    if (std::abs(*ptr)>threshold)
      nnz++;
    }
  return nnz;
}

template <typename T>
inline void 
PrintVector ( const std::vector<T>& vec, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  // utlShowPosition(true);
  std::vector<double> st = GetContainerStats(vec.begin(), vec.end());
  char tmp[1024];
  sprintf(tmp, "%-8s(%p):  size = %lu,  stat = { %g, %g [%g], %g } : ", str==""?"vector":str.c_str(), &vec, vec.size(), st[0], st[2], st[3], st[1] );
  std::string strr(tmp);
  if (vec.size()>0)
    {
    os << strr << " = [ ";
    for ( int i = 0; i < vec.size(); i += 1 ) 
      os << vec[i] << separate;
    os << "];\n" << std::flush;
    }
  else
    os << strr << " is an empty vector" << std::endl;
}

template <>
inline void
PrintVector<std::string> ( const std::vector<std::string>& vec, const std::string str, const char* separate, std::ostream& os)
{
  char tmp[1024];
  sprintf(tmp, "%-8s(%p):  size = %lu, : ", str==""?"vector":str.c_str(), &vec, vec.size());
  std::string strr=tmp;
  if (vec.size()>0)
    {
    os << strr << " = [ ";
    for ( int i = 0; i < vec.size(); i += 1 ) 
      os << vec[i] << separate;
    os << "];\n" << std::flush;
    }
  else
    os << strr << " is an empty vector" << std::endl;
}

template <class VectorType>
inline void 
PrintVector ( const VectorType& vec, const int NSize, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  // utlShowPosition(true);
  char tmp[1024];
  std::vector<double> st = GetContainerStats(&vec[0], &vec[NSize]);
  sprintf(tmp, "%-8s(%p):  size = %lu,  stat = { %g, %g [%g], %g } : ", str==""?"vector":str.c_str(), &vec, (long unsigned int)NSize, st[0], st[2], st[3], st[1] );
  std::string strr(tmp);
  if (NSize>0)
    {
    os << strr << " = [ ";
    for ( int i = 0; i < NSize; i += 1 ) 
      os << vec[i] << separate;
    os << "];\n" << std::flush;
    }
  else
    os << strr << " is an empty vector" << std::endl;
}

template <class IteratorType>
inline void 
PrintContainer ( IteratorType v1, IteratorType v2, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  // utlShowPosition(true);
  char tmp[1024];
  std::vector<double> st = GetContainerStats(v1, v2);
  long unsigned int NSize=0;
  IteratorType q = v1;
  for ( ; q!=v2; ++q ) 
    NSize++;
  sprintf(tmp, "%-8s(%p):  size = %lu,  stat = { %g, %g [%g], %g } : ", str==""?"vector":str.c_str(), v1, NSize, st[0], st[2], st[3], st[1] );

  std::string strr(tmp);
  if (v1!=v2)
    {
    os << strr << " = [ ";
    IteratorType q = v1;
    for ( ; q!=v2; ++q ) 
      os << (*q) << separate;
    os << "];\n" << std::flush;
    }
  else
    os << strr << " is an empty container" << std::endl;
}

inline  int 
strfind(const char *s,const char c) 
{
  if (s) 
    { 
    int l; 
    for (l=0; l<=strlen(s) && s[l]!=c; l++) 
      ;
    if (l==strlen(s)+1)
      l = -1;
    return l; 
    }
  return -1; 
}

// to set std::vector<T> with variable length
template <typename T>
  int
SetVector(const char *s, std::vector<T>& vec, const int least_num=0, const char& c=',') 
{
  vec.clear();
  if (strlen(s)==0)
    return vec.size();
  if (s) 
    { 
    // std::cout << "strlen(s) = " << strlen(s) << std::endl;
    unsigned int l; 
    char temp[1024];
    int k = 0;
    for (l=0; l<=strlen(s); l++)
      {
      if (s[l]==c || s[l]=='\0')
        {
        temp[k]='\0';
        k = 0;
        T ftmp = (T)atof(temp); 
        vec.push_back(ftmp);
        }
      else
        {
        temp[k] = s[l];
        k++;
        }
      }
    if (l==strlen(s)+1)
      l = -1;
    utlAssert(least_num==0 || vec.size()>=least_num, "no enough numbers were imported in vector");
    return l; 
    }
  utlAssert(least_num==0 || vec.size()>=least_num, "no enough numbers were imported in vector");
  return -1; 
}

template <typename T>
  int
SetVector(const std::string s, std::vector<T>& vec, const int least_num=0, const char& c=',') 
{
  SetVector(s.c_str(), vec, least_num, c);
}

template < class T >
inline void
ReadVector ( const std::string vectorStr, std::vector<T>& vec)
{
  utlException(vectorStr=="", "need to set vectorStr");
  vec.clear();

  std::vector<std::vector<std::string> > strVec;
  ReadLinesFirstlineCheck(vectorStr, strVec, " ");
  utlException(strVec.size()==0, "wrong file");

  bool isRowVector = strVec.size()==1;
  if (isRowVector)
    {
    for ( int i = 0; i < strVec[0].size(); i += 1 ) 
      vec.push_back(ConvertStringToNumber<T>(strVec[0][i]));
    }
  else
    {
    for ( int i = 0; i < strVec.size(); i += 1 ) 
      {
      utlException(strVec[i].size()!=1, "should keep one float number in one line");
      vec.push_back(ConvertStringToNumber<T>(strVec[i][0]));
      }
    }
}

template <>
inline void
ReadVector<std::string> ( const std::string vectorStr, std::vector<std::string>& vec)
{
  utlException(vectorStr=="", "need to set vectorStr");
  vec.clear();

  std::vector<std::vector<std::string> > strVec;
  ReadLinesFirstlineCheck(vectorStr, strVec, " ");
  utlException(strVec.size()==0, "wrong file");
  
  bool isRowVector = strVec.size()==1;
  if (isRowVector)
    {
    for ( int i = 0; i < strVec[0].size(); i += 1 ) 
      vec.push_back(strVec[0][i]);
    }
  else
    {
    for ( int i = 0; i < strVec.size(); i += 1 ) 
      {
      utlException(strVec[i].size()!=1, "should keep one float number in one line");
      vec.push_back(strVec[i][0]);
      }
    }
}

template <typename VectorType>
void 
SaveVector ( const VectorType& vv, const int NSize, const std::string vectorStr, const bool is_save_number=false)
{
  std::ofstream  out;                                // create ofstream object
  out.open ( vectorStr.c_str() );           // open ofstream
  utlException(!out, "\nERROR : failed to open output file " << vectorStr );
  
  if (is_save_number)
    out << NSize << "\n";
  for ( int i = 0; i < NSize; i += 1 ) 
    out << vv[i] << "\n";
  out.close();
}

template <typename T>
void 
SaveVector ( const std::vector<T>& vv, const std::string vectorStr, const bool is_save_number=false)
{
  int NSize = vv.size();
  SaveVector<std::vector<T> >(vv, NSize, vectorStr, is_save_number);
}

template <typename T>
inline std::vector<int> 
FindVector ( const std::vector<T>& vec, const T elem, const double gap=M_EPS )
{
  std::vector<int> index;
  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    if (std::fabs(elem-vec[i])<gap)
      index.push_back(i);
    }
  return index;
}

template <>
inline std::vector<int> 
FindVector<std::string> ( const std::vector<std::string>& vec, const std::string elem, const double )
{
  std::vector<int> index;
  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    if (elem==vec[i])
      index.push_back(i);
    }
  return index;
}

template <typename T>
std::vector<std::vector<int> >
SeparateVector ( const std::vector<T>& vec, std::vector<T>& vec_sep, const double gap=M_EPS )
{
  std::vector<std::vector<int> > index;
  vec_sep.clear();

  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    bool is_test = false;
    for ( int i_ind = 0; i_ind < index.size(); i_ind += 1 ) 
      {
      for ( int j_ind = 0; j_ind < index[i_ind].size(); j_ind += 1 ) 
        {
        if (i==index[i_ind][j_ind])
          {
          is_test = true;
          break;
          }
        }
      if (is_test)
        break;
      }

    if (!is_test)
      {
      vec_sep.push_back(vec[i]);
      std::vector<int> index_tmp = FindVector(vec, vec[i], gap);
      index.push_back(index_tmp);
      }
    }

  for ( int m = 0; m < index.size(); m += 1 ) 
    {
    T sum = 0.0;
    for ( int n = 0; n < index[m].size(); n += 1 ) 
      sum += vec[index[m][n]];
    vec_sep[m] = sum/(double)index[m].size();
    }
  return index;
}

template <typename T>
void
ConnectVector ( std::vector<T>& vec1, const std::vector<T>& vec2)
{
  for ( int i = 0; i < vec2.size(); i += 1 ) 
    vec1.push_back(vec2[i]);
}

template <typename T>
std::vector<T> 
SelectVector ( const std::vector<T>& vec, const std::vector<int>& index )
{
  std::vector<T> select;
  for ( int i = 0; i < index.size(); i += 1 ) 
    {
    select.push_back( vec[index[i]] );
    }
  return select;
}

template <typename T>
std::vector<T> 
SelectVector ( const std::vector<T>& vec, const int startIndex, const int numberOfElement )
{
  std::vector<T> select;
  for ( int i = 0; i < numberOfElement; i += 1 ) 
    {
    if (startIndex+i<vec.size())
      select.push_back( vec[startIndex+i] );
    }
  return select;
}

template <class VectorType, class T>
inline bool 
IsInVector ( const VectorType& vec, const int size, const T& num, const double eps=1e-10 )
{
  for ( int i = 0; i < size; i += 1 ) 
    {
    if (std::abs(vec[i]-num) < eps)
      return true;
    }
  return false;
}

template <typename T>
inline bool 
IsInVector ( const std::vector<T>& vec, const T& num, const double eps=1e-10 )
{
  return IsInVector<std::vector<T>, T>(vec, vec.size(), num, eps);
}

template <>
inline bool 
IsInVector<int> ( const std::vector<int>& vec, const int& num, const double )
{
  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    if (vec[i]==num)
      return true;
    }
  return false;
}

template <>
inline bool 
IsInVector<std::string> ( const std::vector<std::string>& vec, const std::string& num, const double )
{
  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    if (vec[i]==num)
      return true;
    }
  return false;
}

template <>
inline bool 
IsInVector<char> ( const std::vector<char>& vec, const char& num, const double )
{
  for ( int i = 0; i < vec.size(); i += 1 ) 
    {
    if (vec[i]==num)
      return true;
    }
  return false;
}

/** generate random float value in [d1,d2]  */
template < typename T >
inline T 
Random( const T d1 = T(0.0), const T d2 = T(1.0) )
{
  utlAssert(d2>d1, "wrong setting");
  static bool first_time = true;
  if (first_time) 
    {
    std::srand(std::time(NULL)); 
    first_time = false; 
    }
  T ra = 1e-10 +  (1-2e-10) * std::rand() / RAND_MAX;
  return (d2 - d1)*ra + d1;
}  // -----  end random  -----

/** generate random int value in [d1,d2]  */
inline long 
RandomInt( const long d1 = 0, const long d2 = 1 )
{
  if (d2==d1)
    return d1;
  utlAssert(d2>d1, "wrong setting, d1="<<d1 << ", d2="<<d2);
  static bool first_time = true;
  if (first_time) 
    {
    std::srand(std::time(NULL)); 
    first_time = false; 
    }
  int tmp = d2-d1;
  int rand_tmp = std::rand() % (tmp+1);
  return rand_tmp + d1;
}  // -----  end random  -----

template < typename T >
inline std::vector<T> 
RandomVec( const std::vector<T>& vec, const int num, std::vector<int>& result_index)
{
  utlAssert(vec.size()>num && num>0, "wrong size, vec.size()="<<vec.size() << ", num="<<num);
  result_index.clear();
  std::vector<T> result;
  static bool first_time = true;
  if (first_time) 
    {
    std::srand(std::time(NULL)); 
    first_time = false; 
    }
  for ( int nn = 0; nn < num; nn += 1 ) 
    {
    int i = RandomInt(0, vec.size()-num+nn);
    bool is_choose = false;
    for ( int n_i = 0; n_i < result_index.size(); n_i += 1 ) 
      {
      if (i==result_index[n_i])
        {
        is_choose=true;
        break;
        }
      }
    if (!is_choose)
      {
      result_index.push_back(i);
      result.push_back(vec[i]);
      }
    else
      {
      result_index.push_back(vec.size()-num+nn);
      result.push_back(vec[vec.size()-num+nn]);
      }
    }
  return result;
}  // -----  end random  -----

/** generate random point in sphere. 
 * http://mathworld.wolfram.com/SpherePointPicking.html 
 * */
inline std::vector<double> 
RandomPointInSphere ( const bool hemis )
{
  std::vector<double> result(3);
  double u = hemis ? Random(0.0,1.0) : Random(-1.0,1.0); // hemisphere
  double theta = Random(0.0,2*M_PI); // hemisphere
  result[0] = std::sqrt(1-u*u)*std::cos(theta);
  result[1] = std::sqrt(1-u*u)*std::sin(theta);
  result[2] = u;
  return result;
}

/** Generate random value from Gaussian distribution.  
 * http://en.wikipedia.org/wiki/Gaussian_random_variable#Generating_values_from_normal_distribution  */
template < typename T >
inline T 
GaussRand ( const T value, const double sigma )
{
  if (sigma<=0)
    return value;
  static bool first_time = true;
  if (first_time) 
    {
    std::srand(std::time(NULL)); 
    first_time = false; 
    }
  T A;
  T U1,U2;
  U1 = 1e-10 +  (1-2e-10) * std::rand() / RAND_MAX;
  U2 = 1e-10 +  (1-2e-10) * std::rand() / RAND_MAX;  	  
  A = std::sqrt(-2*std::log(U1))*std::cos(2*M_PI*U2);
  return (value + sigma*A);
}  // -----  end rician_rand  -----

/**  Generate random value from Rician distribution.  
 * */
template < typename T >
inline T 
RicianRand ( const T value, const double sigma )
{
  if (sigma<=0)
    return value;
  static bool first_time = true;
  if (first_time) 
    {
    std::srand(std::time(NULL)); 
    first_time = false; 
    }

  T real_part = value;
  T imag_part = 0;

  T n1 = GaussRand(T(0.0),sigma);
  T n2 = GaussRand(T(0.0),sigma);

  T noise_real = real_part + n1;
  T noise_imag = imag_part + n2;

  T noise = std::sqrt(noise_real*noise_real + noise_imag*noise_imag);
  return noise;
}  // -----  end rician_rand  -----

template < typename VectorType>
inline VectorType
AddNoise( const VectorType& signal, const int size, const double sigma, const bool is_rician=true)
{
  if (sigma<=0)
    return signal;

  VectorType result(signal);
  for ( int i = 0; i < size; i += 1 ) 
    {
    if (is_rician)
      result[i] = RicianRand(signal[i],sigma);
    else
      result[i] = GaussRand(signal[i],sigma);
    }
  return result;
}

template < typename T >
inline std::vector<T>
AddNoise( const std::vector<T>& signal, const double sigma, const bool is_rician=true)
{
  return AddNoise(signal, signal.size(), sigma, is_rician);
}


/** http://en.wikipedia.org/wiki/Spherical_coordinate_system  */
template < typename T >
inline void
cartesian2Spherical( const T x, const T y, const T z, T& r, T & theta, T & phi)
{
  r = sqrt(x*x + y*y + z*z);
  if ( r >= M_EPS )
    {
    theta   = std::acos( z / r );
    phi = std::atan2( y, x );
    }
  else
    theta = phi = 0.0;
}

template < typename T >
inline void
cartesian2Spherical( T& x, T& y, T& z)
{
  T x1=x, y1=y, z1=z;
  cartesian2Spherical(x1,y1,z1, x,y,z);
}

/** http://en.wikipedia.org/wiki/Spherical_coordinate_system  */
template < typename T >
inline void
spherical2Cartesian( const T r, const T theta, const T phi, T& x, T & y, T & z)
{
  x = r * std::sin(theta) * std::cos(phi);
  y = r * std::sin(theta) * std::sin(phi);
  z = r * std::cos(theta);
}

template < typename T >
inline void
spherical2Cartesian( T& x, T& y, T& z)
{
  T x1=x, y1=y, z1=z;
  spherical2Cartesian(x1,y1,z1, x,y,z);
}

template <class IteratorType>
inline void 
PowerVector (  IteratorType v1, IteratorType v2, const double poww )
{
  for ( IteratorType ptr=v1; ptr != v2;  ++ptr ) 
    *ptr = std::pow(*ptr, poww);
}

template <class VectorType>
inline VectorType
GetVectorLinspace( const double valueMin, const double valueMax, const int num )
{
  VectorType result(num);
  if (num==1)
    result[0] = valueMin;
  else
    {
    double step = (valueMax-valueMin)/(double)(num-1);
    for ( int i = 0; i < num; i += 1 ) 
      result[i] = valueMin+step*i;
    }
  return result;
}

inline std::vector<int> 
GetRange ( const int start, const int end, const int space=1)
{
  std::vector<int> result;
  int val = start;
  do 
    {
    result.push_back(val);
    val += space;
    } while ( val < end );
  return result;
}

template <class VectorType>
inline void
VectorShrinkage( VectorType& vec, const int N, const double kappa )
{
  for ( int i = 0; i < N; i += 1 ) 
    vec[i] = (vec[i]>=kappa) ? (vec[i]-kappa) : ( (vec[i]<=-kappa) ? (vec[i]+kappa) : 0.0 );
}

template <class VectorType>
inline VectorType
GetVectorShrinkage( const VectorType& vec, const int N, const double kappa )
{
  VectorType result(vec);
  VectorShrinkage(result, N, kappa);
  return result;
}

template <class VectorType>
inline void
AbsoluteVector( VectorType& vec, const int N)
{
  for ( int i = 0; i < N; i += 1 ) 
    vec[i] = std::abs(vec[i]);
}

template <class VectorType>
inline VectorType
GetAbsoluteVector( const VectorType& vec, const int N )
{
  VectorType result(vec);
  AbsoluteVector(result, N);
  return result;
}

/** copy a vector to another vector with different type.  */
template <class T1, class T2>
void 
VectorToVector ( const T1& v1,  T2& v2, const int N)
{
  for ( int i = 0; i < N; i += 1 ) 
    v2[i] = v1[i];
}

/** copy a matrix to another matrix with different type.  */
template <class T1, class T2>
void 
MatrixToMatrix ( const T1& mat1,  T2& mat2, const int NRows, const int NColumns)
{
  for ( int i = 0; i < NRows; i += 1 ) 
    for ( int j = 0; j < NColumns; j += 1 ) 
    mat2(i,j) = mat1(i,j);
}

/** print statistics from a matrix  */
template <class TMatrixType>
void
PrintMatrixStats ( const TMatrixType& matrix, const int NumberRows, const int NumberColumns, const std::string str="", const char* separate=" ", std::ostream& os=std::cout )
{
  char tmp[1024];
  std::vector<double> st = NumberRows*NumberColumns>0 ? GetContainerStats(&matrix(0,0), &matrix(0,0)+NumberRows*NumberColumns) : std::vector<double>(4,0);
  sprintf(tmp, "%-8s(%p):  size = (%d, %d),  stat = { %g, %g [%g], %g } : ", str==""?"matrix":str.c_str(), &matrix, NumberRows, NumberColumns, st[0], st[2], st[3], st[1] );
  std::string strr(tmp);
  if (NumberRows>0 && NumberColumns>0)
    os << strr << std::endl;
}

template <class TMatrixType>
void
PrintMatrix ( const TMatrixType& matrix, const int NumberRows, const int NumberColumns, const std::string str="", const char* separate=" ", std::ostream& os=std::cout )
{
  char tmp[1024];
  std::vector<double> st = NumberRows*NumberColumns>0 ? GetContainerStats(&matrix(0,0), &matrix(0,0)+NumberRows*NumberColumns) : std::vector<double>(4,0);
  sprintf(tmp, "%-8s(%p):  size = (%d, %d),  stat = { %g, %g [%g], %g } : ", str==""?"matrix":str.c_str(), &matrix, NumberRows, NumberColumns, st[0], st[2], st[3], st[1] );
  std::string strr(tmp);
  
  if (NumberRows>0 && NumberColumns>0)
    {
    os << strr << " = [ ";
    for ( int i = 0; i < NumberRows-1; i += 1 ) 
      {
      for ( int j = 0; j < NumberColumns-1; j += 1 ) 
        os << matrix(i,j)<<  separate;
      os << matrix(i,NumberColumns-1) << "; ";
      }
    for ( int j = 0; j < NumberColumns-1; j += 1 ) 
      os << matrix(NumberRows-1,j)<<  separate;
    os << matrix(NumberRows-1,NumberColumns-1) << "];\n" << std::flush;
    }
  else
    os << strr << " is an empty matrix" << std::endl;
  return;
}

template <class TMatrixType>
inline double
ArgmaxSymmetricMatrix(const TMatrixType matrix, const int size, int& row, int& colomn, const bool includeDiagonalElements)
{
  double maxValue = -std::numeric_limits<double>::max();
  double value;
  for ( unsigned int i = 0; i < size; i += 1 ) 
    {
    for ( unsigned int j = 0; (includeDiagonalElements?(j<=i):(j<i)) ; j += 1 ) 
      {
      value = matrix(i,j);
      if (value > maxValue)
        {
        maxValue = value;
        row = i;
        colomn = j;
        }
      }
    }
  return maxValue;
}

template <class TMatrixType>
inline double
ArgminSymmetricMatrix(const TMatrixType& matrix, const int size, int& row, int& colomn, const bool includeDiagonalElements)
{
  double minValue = std::numeric_limits<double>::max();
  double value;
  for ( unsigned int i = 0; i < size; i += 1 ) 
    {
    for ( unsigned int j = 0; (includeDiagonalElements?(j<=i):(j<i)) ; j += 1 ) 
      {
      value = matrix(i,j);
      if (value < minValue)
        {
        minValue = value;
        row = i;
        colomn = j;
        }
      }
    }
  return minValue;
}

template <class TMatrixType>
inline double
ArgmaxMatrix(const TMatrixType& matrix, const int rows, const int columns, int& row, int& colomn)
{
  double maxValue = -std::numeric_limits<double>::max();
  double value;
  for ( unsigned int i = 0; i < rows; i += 1 ) 
    {
    for ( unsigned int j = 0; j < columns; j += 1 ) 
      {
      value = matrix(i,j);
      if (value > maxValue)
        {
        maxValue = value;
        row = i;
        colomn = j;
        }
      }
    }
  return maxValue;
}

template <class TMatrixType>
inline double
ArgminMatrix(const TMatrixType matrix, const int rows, const int columns, int& row, int& colomn)
{
  double minValue = std::numeric_limits<double>::max();
  double value;
  for ( unsigned int i = 0; i < rows; i += 1 ) 
    {
    for ( unsigned int j = 0; j < columns; j += 1 ) 
      {
      value = matrix(i,j);
      if (value < minValue)
        {
        minValue = value;
        row = i;
        colomn = j;
        }
      }
    }
  return minValue;
}

template <typename Vector2D>
void 
Save2DVector ( const Vector2D& vv, std::ostream& out=std::cout)
{
  for ( int i = 0; i < vv.size(); i += 1 ) 
    {
    for ( int j = 0; j < vv[i].size()-1; j += 1 ) 
      out << vv[i][j]<<  " ";
    out << vv[i][vv[i].size()-1] << "\n";
    }
  return;
}

template <typename Vector2D>
void 
Save2DVector ( const Vector2D& vv, const std::string file)
{
  std::ofstream  out;                                // create ofstream object
  out.open ( file.c_str() );           // open ofstream
  utlException (!out, "\nERROR : failed to open output file " << file );
  Save2DVector(vv, out);
  out.close();
  return;
}

template <class TMatrixType>
void
SaveMatrix ( const TMatrixType& matrix, const int NumberRows, const int NumberColumns, const std::string file )
{
  std::ofstream  out;                                // create ofstream object
  out.open ( file.c_str() );           // open ofstream
  utlException (!out, "\nERROR : failed to open output file " << file );
  for ( int i = 0; i < NumberRows; i += 1 ) 
    {
    for ( int j = 0; j < NumberColumns-1; j += 1 ) 
      out << matrix(i,j)<<  " ";
    out << matrix(i,NumberColumns-1) << "\n";
    }
  out.close();
  return;
}

template <class TMatrixType>
void
ReadMatrix ( const std::string file, TMatrixType& matrix )
{
  std::vector<std::vector<std::string> > matrixStr;
  utl::ReadLines(file, matrixStr);
  int NumberRows = matrixStr.size();
  utlException(NumberRows==0, "wrong file");
  int NumberColumns = matrixStr[0].size();
  matrix = TMatrixType(NumberRows, NumberColumns);
  for ( int i = 0; i < NumberRows; i += 1 ) 
    {
    utlException(NumberColumns!=matrixStr[i].size(), "wrong size");
    for ( int j = 0; j < NumberColumns; j += 1 ) 
      std::istringstream ( matrixStr[i][j] ) >> matrix(i,j);
    }
  return;
}

}

#endif 
