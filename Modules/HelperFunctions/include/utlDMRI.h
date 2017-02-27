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
#include "utlNDArray.h"
#include "utlSTDHeaders.h"

/** @addtogroup utlHelperFunctions
@{ */

enum {
  CARTESIAN_TO_CARTESIAN=0,
  CARTESIAN_TO_SPHERICAL=1,
  SPHERICAL_TO_CARTESIAN=2,
  SPHERICAL_TO_SPHERICAL=3
};

enum {
  RADIAN_TO_DEGREE=0,
  RADIAN_TO_RADIAN=1,
  DEGREE_TO_RADIAN=2,
  DEGREE_TO_DEGREE=3
};

enum {
  DIRECTION_NODUPLICATE=0,
  DIRECTION_DUPLICATE=1
};

enum {
  DIRECTION_NOFLIP=0, 
  DIRECTION_FLIP=1
};

// enum {
//   CARTESIAN=0, 
//   SPHERICAL=1
// };

enum{
  // [xx, xy, xz, yy, yz, zz], fsl
  TENSOR_UPPER_TRIANGULAR=0,
  // [xx, yx, yy, zx, zy, zz]
  TENSOR_LOWER_TRIANGULAR, 
  // [xx, yy, zz, sqrt(2)*xy, sqrt(2)*xz, sqrt(2)*yz], embed 3x3 matrix into 6x1 vector
  TENSOR_EMBED6D 
};

namespace utl
{


/** get the maximal rank of spherical harmonic coefficient vector with a given dimension  */
inline int
DimToRankSH ( const int dimm ) 
{
  double result = -1.5 + std::sqrt(0.25 + 2.0*dimm);
  
  utlGlobalAssert(IsInt(result), "Logical ERROR! Wrong dimension! dimm = " << dimm);
  
  return (int) result;
}

/** get the dimension of Spherical Harmonic coefficients from given rank  */
inline int
RankToDimSH (const int shRank)
{
  utlSAGlobalException(!utl::IsEven(shRank))(shRank).msg("shRank should be even");
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

/** 
 * If isE2E3Equal==true, for symmetric tensor with eigenvalues (e1,e2,e2), get (e1,e2) from the fa and (e1+e2+e2)/3.
 * If isE2E3Equal==false, for symmetric tensor with eigenvalues (e1,e1,e2), get (e1,e2) from the fa and (e1+e1+e2)/3.
 * \note: when isE2E3Equal==false, fa cannot be larger than 1/sqrt(2), otherwise e3 is negative. 
 * */
template < typename T >
std::vector<T>
GetE1E2FromFAMD ( const T fa, const T meanEigenValue, const bool isE2E3Equal=true )
{
  std::vector<T> eigenValues(2);
  double tmp = fa/std::sqrt(3.0-2.0*fa*fa);
  if (isE2E3Equal)
    {
    eigenValues[0] = meanEigenValue*(1.0 + 2.0*tmp);
    eigenValues[1] = meanEigenValue*(1.0 - tmp);
    }
  else
    {
    utlSAException(fa>1.0/utl::SQRT2)(fa).msg("fa is too large so that e3 is negative");
    eigenValues[0] = meanEigenValue*(1.0 + tmp);
    eigenValues[1] = meanEigenValue*(1.0 - 2.0*tmp);
    }
  return eigenValues;
}

/** Covert 6D tensor format to 9D format (3x3 symmetric matrix).  
 * v9d is in ROW_MAJOR.*/
template<typename V1Type, typename V2Type>
void 
ConvertTensor6DTo9D(const V1Type& v6d, V2Type& v9d, int v6dStoreWay)
{
  if (v6dStoreWay==TENSOR_LOWER_TRIANGULAR)
    {
    v9d[0] = v6d[0]; // xx
    v9d[1] = v6d[1]; // xy
    v9d[2] = v6d[3]; // xz
    v9d[3] = v6d[1]; // yx
    v9d[4] = v6d[2]; // yy
    v9d[5] = v6d[4]; // yz
    v9d[6] = v6d[3]; // zx
    v9d[7] = v6d[4]; // zy
    v9d[8] = v6d[5]; // zz
    }
  else if (v6dStoreWay==TENSOR_UPPER_TRIANGULAR)
    {
    v9d[0] = v6d[0]; // xx
    v9d[1] = v6d[1]; // xy
    v9d[2] = v6d[2]; // xz
    v9d[3] = v6d[1]; // yx
    v9d[4] = v6d[3]; // yy
    v9d[5] = v6d[4]; // yz
    v9d[6] = v6d[2]; // zx
    v9d[7] = v6d[4]; // zy
    v9d[8] = v6d[5]; // zz
    }
  else if (v6dStoreWay==TENSOR_EMBED6D)
    {
    v9d[0] = v6d[0]; // xx
    v9d[1] = v6d[3]*utl::SQRT1_2; // xy
    v9d[2] = v6d[4]*utl::SQRT1_2; // xz
    v9d[3] = v6d[3]*utl::SQRT1_2; // yx
    v9d[4] = v6d[1]; // yy
    v9d[5] = v6d[5]*utl::SQRT1_2; // yz
    v9d[6] = v6d[4]*utl::SQRT1_2; // zx
    v9d[7] = v6d[5]*utl::SQRT1_2; // zy
    v9d[8] = v6d[2]; // zz
    }
  else
    utlException(true, "wrong v1StoreWay");
}

/** Covert 9D tensor format to 6D format (3x3 symmetric matrix).  
 * v9d is in ROW_MAJOR.*/
template<typename V1Type, typename V2Type>
void 
ConvertTensor9DTo6D(const V1Type& v9d, V2Type& v6d, int v6dStoreWay)
{
  if (v6dStoreWay==TENSOR_LOWER_TRIANGULAR)
    {
    v6d[0] = v9d[0]; // xx
    v6d[1] = v9d[3]; // yx
    v6d[2] = v9d[4]; // yy
    v6d[3] = v9d[6]; // zx
    v6d[4] = v9d[7]; // zy
    v6d[5] = v9d[8]; // zz
    }
  else if (v6dStoreWay==TENSOR_UPPER_TRIANGULAR)
    {
    v6d[0] = v9d[0]; // xx
    v6d[1] = v9d[1]; // xy
    v6d[2] = v9d[2]; // xz
    v6d[3] = v9d[4]; // yy
    v6d[4] = v9d[5]; // yz
    v6d[5] = v9d[8]; // zz
    }
  else if (v6dStoreWay==TENSOR_EMBED6D)
    {
    v6d[0] = v9d[0]; // xx
    v6d[1] = v9d[4]; // yy
    v6d[2] = v9d[8]; // zz
    v6d[3] = v9d[1]*utl::SQRT2; // sqrt(2)*xy
    v6d[4] = v9d[2]*utl::SQRT2; // sqrt(2)*xz
    v6d[5] = v9d[5]*utl::SQRT2; // sqrt(2)*yz
    }
  else
    utlException(true, "wrong v1StoreWay");
}

/** The function only works for 2th order tensors which are in 3D and are symmetric.  */
template <class T>
inline void
Convert2To1Tensor(const utl::NDArray<T,2>& mat, utl::NDArray<T,1>& vec)
{
  utlException(mat.Rows()!=3 || mat.Cols()!=3, "It should be in 3D space");
  vec.ReSize(6);
  utl::ConvertTensor9DTo6D(mat, vec, TENSOR_EMBED6D);
}

template <class T>
inline void
Convert1To2Tensor(const utl::NDArray<T,1>& vec, utl::NDArray<T,2>& mat)
{
  utlException(vec.Size()!=6, "It should be in 3D space");
  mat.ReSize(3,3);
  utl::ConvertTensor6DTo9D(vec, mat, TENSOR_EMBED6D);
}

template <class T>
inline void
NormalizeGrad(const utl::NDArray<T,2>& grad)
{
  if (grad.Size()==0)
    return;
  utlSAException(grad.Cols()!=3)(grad.Rows())(grad.Cols())("grad should be Nx3 matrix");
  utl::NDArray<T,1> vec;
  for ( int i = 0; i < grad.Rows(); ++i ) 
    {
    vec = grad.GetRow(i);
    double norm = vec.GetTwoNorm();
    if (norm>1e-10)
      {
      vec /= norm;
      grad.SetRow(i, vec);
      }
    }
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

    /** @} */

#endif 
