/**
 *       @file  utlDMRI.h
 *      @brief  
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlDMRI_h
#define __utlDMRI_h

#include "utlCore.h"
#include "utlMatrix.h"
#include "utlSTDHeaders.h"


enum {
  CARTESIAN_TO_CARTESIAN=0,
  CARTESIAN_TO_SPHERICAL=1,
  SPHERICAL_TO_CARTESIAN=3,
  SPHERICAL_TO_SPHERICAL=4
};

enum {
  DIRECTION_NODUPLICATE=0,
  DIRECTION_DUPLICATE=1
};

enum {
  DIRECTION_NOFLIP=0, 
  DIRECTION_FLIP=1
};

namespace utl
{


/** get the maximal rank of spherical harmonic coefficient vector with a given dimension  */
inline int
DimToRankSH ( const int dimm ) 
{
  double result = -1.5 + std::sqrt(0.25 + 2.0*dimm);
  
  utlAssert(IsInt(result), "Logical ERROR! Wrong dimension! dimm = " << dimm);
  
  return (int) result;
}

/** get the dimension of Spherical Harmonic coefficients from given rank  */
inline int
RankToDimSH (const int shRank)
{
  if (shRank<0)
    return 0;
  else
    return (shRank + 1)*(shRank + 2)/2;
}

/** get SH index j from given rank (l,m). [0, -2,-1,0,1,2,  -4,-3,-2,-1,0,1,2,3,4, ...]   */
inline int 
GetIndexSHj(const int l, const int m) 
{
  return (l*l + l)/2 + m;
}

/** get rank index (l,m) from given linear index j. [0, -2,-1,0,1,2,  -4,-3,-2,-1,0,1,2,3,4, ...]  */
inline std::vector<int> 
GetIndexSHlm(const int j) 
{
  std::vector<int> lm(2,0);
  
  if (j==0)
    return lm;

  int l = std::ceil(0.5*(std::sqrt(9.+8.*j)-3.0));
  if (utl::IsOdd(l))
    l+=1;
  int residual = j+1 - (l-1)*l/2;
  lm[0] = l;
  lm[1] = residual-1 - l;

  return lm;
}

/** for symmetric tensor with eigenvalues (e1,e2,e2), get (e1,e2) from the fa and (e1+e2+e3)/3  */
template < typename T >
std::vector<T>
GetE1E2FromFAMD ( const T fa, const T meanEigenValue )
{
  std::vector<T> e1e2(2);
  e1e2[1] = meanEigenValue*(1.0 - fa/std::sqrt(3.0-2.0*fa*fa));
  e1e2[0] = 3.0*meanEigenValue - 2.0*e1e2[1];
  return e1e2;
}

template <class T>
inline utl_shared_ptr<NDArray<T,2> >
ReadGrad(const std::string grad_str, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode=CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
{
  std::vector < std::vector<std::string> > strVec;

  // use three delimiters ' ', ',', '\t'  
  ReadLinesFirstlineCheck(grad_str, strVec, " \t,");

  for ( int i = 0; i < strVec.size(); i += 1 ) 
    utlException(strVec[i].size()!=3, "wrong dimension in gradient file! strVec[i].size()="<< strVec[i].size() <<", grad_str="<<grad_str);

  utl_shared_ptr<NDArray<T,2> > grad (new NDArray<T,2>() );
  if (NoSymmetricDuple==DIRECTION_DUPLICATE)
    grad->ReSize(2*strVec.size(),3);
  else
    grad->ReSize(strVec.size(),3);

  std::vector<double> xyz(3);
  int jj=0;
  for (int i=0; i<strVec.size(); i++) 
    {
    utlException(strVec[i].size()!=3, "wrong size for grad");
    std::istringstream ( strVec[i][0] ) >> xyz[0];
    std::istringstream ( strVec[i][1] ) >> xyz[1];
    std::istringstream ( strVec[i][2] ) >> xyz[2];

    if (mode==CARTESIAN_TO_SPHERICAL || mode==CARTESIAN_TO_CARTESIAN)
      {
      if (need_normalize)
        xyz = NormalizeUnitNorm(xyz);
      xyz[0] = (flipx==DIRECTION_FLIP?-1:1)*xyz[0]; xyz[1] = (flipy==DIRECTION_FLIP?-1:1)*xyz[1]; xyz[2] = (flipz==DIRECTION_FLIP?-1:1)*xyz[2];
      if (mode==CARTESIAN_TO_SPHERICAL)
        cartesian2Spherical(xyz[0],xyz[1],xyz[2]);
      }
    else if(mode==SPHERICAL_TO_CARTESIAN || mode==SPHERICAL_TO_SPHERICAL)
      {
      spherical2Cartesian(xyz[0],xyz[1],xyz[2]);
      if (need_normalize)
        xyz = NormalizeUnitNorm(xyz);
      xyz[0] = (flipx==DIRECTION_FLIP?-1:1)*xyz[0]; xyz[1] = (flipy==DIRECTION_FLIP?-1:1)*xyz[1]; xyz[2] = (flipz==DIRECTION_FLIP?-1:1)*xyz[2];
      if (mode==SPHERICAL_TO_SPHERICAL)
        cartesian2Spherical(xyz[0],xyz[1],xyz[2]);
      }
    else
      utlGlobalException(true, "wrong mode");

    (*grad)(jj,0)=xyz[0], (*grad)(jj,1)=xyz[1], (*grad)(jj,2)=xyz[2];
    jj++;

    if(NoSymmetricDuple==DIRECTION_DUPLICATE) 
      {
      if (mode==SPHERICAL_TO_CARTESIAN || mode==CARTESIAN_TO_CARTESIAN)
        { xyz[0] = -xyz[0]; xyz[1] = -xyz[1]; xyz[2] = -xyz[2]; }
      else
        { xyz[1]=M_PI-xyz[1]; xyz[2]=M_PI+xyz[2]; xyz[2]=(xyz[2]>=2*M_PI)?(xyz[2]-2*M_PI):xyz[2]; }
      (*grad)(jj,0)=xyz[0], (*grad)(jj,1)=xyz[1], (*grad)(jj,2)=xyz[2];
      jj++;
      }
    }

  return grad;
}

template < class T >
void
MatchBVectorAndGradientMatrix (const T& br, std::vector<T>& vec, const NDArray<T,2>& grad )
{
  utlException(grad.Rows()==0, "grad.size()==0");
  vec = std::vector<T>(grad.Rows(), br);
}

template < class T >
void
MatchBVectorAndGradientMatrix ( std::vector<T>& vec, NDArray<T,2>& grad )
{
  utlException(vec.size()==0 || grad.Rows()==0, "wrong size! vec.size()="<<vec.size()<<", grad.Rows()="<<grad.Rows());
  
  if (grad.Rows()==vec.size())
    return;

  std::vector<T> vec_result;
  if (vec.size()==1)
    {
    MatchBVectorAndGradientMatrix(vec[0], vec_result, grad);
    vec = vec_result;
    return;
    }

  vec_result.resize(vec.size()* grad.Rows());
  NDArray<T,2> grad_result(vec.size()* grad.Rows(), grad.Columns());
  for (int i = 0; i<vec.size(); i++)
    {
    for ( int j = 0; j < grad.Rows(); j += 1 ) 
      {
      vec_result[i*grad.Rows()+j] = vec[i];
      for ( int kk = 0; kk < 3; kk += 1 ) 
        grad_result(i*grad.Rows()+j,kk) = grad(j,kk);
      }
    }
  vec = vec_result;
  grad = grad_result;
  return;
}

/* Function taken from 3D Slicer, SuperquadricTensorGlyph
 *
 * This is sort of the inverse of code from Gordon Kindlmann for mapping
 * the mode (index value) to RGB. See vtkTensorMathematics for that code.
 * There may be a simpler way to do this but this works.
 * Note this expects a "0 1" Hue Range in the vtkLookupTable used to
 * display the glyphs.
 */
inline void 
RGBToIndex(double R, double G, double B, double &index) 
{

  // remove the gray part of the color.
  // this is so we can use the model where either R,G, or B is 0.
  // then we scale so that the max of the other two is one.
  double min = R;
  int minIdx = 0;
  
  if (G < min)
    {
      min = G;
      minIdx = 1;
    }
  if (B < min)
    {
      min = B;
      minIdx = 2;
    }

  // make the smallest of R,G,B equal 0
  R = R - min;
  G = G - min;
  B = B - min;

  // now take the max, and scale it to be 1.
  double max = R;
  int maxIdx = 0;
  if (G > max)
    {
      max = G;
      maxIdx = 1;
    }
  if (B > max)
    {
      max = B;
      maxIdx = 2;
    }

  R = R/max;
  G = G/max;
  B = B/max;


  // now using the inverse sextants, map this into an index.
  // switch (sextant) {
  //   case 0: { R = 1;      G = frac;   B = 0;      break; }
  //   case 1: { R = 1-frac; G = 1;      B = 0;      break; }
  //   case 2: { R = 0;      G = 1;      B = frac;   break; }
  //   case 3: { R = 0;      G = 1-frac; B = 1;      break; }
  //   case 4: { R = frac;   G = 0;      B = 1;      break; }
  //   case 5: { R = 1;      G = 0;      B = 1-frac; break; }
  // }
  int sextant;
  if (maxIdx == 0 && minIdx == 2) sextant = 0;
  if (maxIdx == 1 && minIdx == 2) sextant = 1;
  if (maxIdx == 1 && minIdx == 0) sextant = 2;
  if (maxIdx == 2 && minIdx == 0) sextant = 3;
  if (maxIdx == 2 && minIdx == 1) sextant = 4;
  if (maxIdx == 0 && minIdx == 1) sextant = 5;

  double offset;
  offset = 256/6;

  switch (sextant) 
    {
    case 0: { index =  G*offset;     break; }
    case 1: { index = offset + (1-R)*offset;      break; }
    case 2: { index = offset*2 + B*offset;   break; }
    case 3: { index = offset*3 + (1-G)*offset;      break; }
    case 4: { index = offset*4 + R*offset;      break; }
    case 5: { index = offset*5 + (1-B)*offset; break; }
    }  
}

}


#endif 
