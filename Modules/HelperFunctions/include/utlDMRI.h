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
  // 9d [xx, xy, xz, yx, yy, yz, zx, zy, zz]
  TENSOR_9D=0,
  // [xx, xy, xz, yy, yz, zz], fsl
  TENSOR_UPPER_TRIANGULAR,
  // [xx, yx, yy, zx, zy, zz]
  TENSOR_LOWER_TRIANGULAR, 
  // [xx, yy, zz, xy, xz yz],  mrtrix, internal in vistasoft, AFQ
  TENSOR_DIAGONAL_FIRST, 
  // [xx, yy, zz, sqrt(2)*xy, sqrt(2)*xz, sqrt(2)*yz], embed 3x3 matrix into 6x1 vector
  TENSOR_EMBED6D 
};

enum{
  TRACTS_UNKNOWN=0,
  TRACTS_TRK=1,
  TRACTS_TCK=2,
  TRACTS_VTK=3
};

enum{
  DWI_NORMALIZE=0,
  DWI_NONORMALIZE
};

namespace utl
{


// [>* The name and dimention map for vectors.  <]
// typedef utl_unordered_map<std::string, int>  VecDimType;
// static VecDimType vectorNameDimMap{
//     {"scalar",1},
//     {"frame", 9}, 
//     {"point",3}
// };  

/** 
* We define a name of scalar vector with its dimension, if the dimension is not 1. 
* For example, the name "frame_9" means a frame with 9 dimensions. 
* The name "orientational order" means a scalar value with 1 dimension. 
* It is used to determine the dimension of scalars in a fiber track.
* */
inline int
GetScalarsDimentionByName(const std::string& name)
{
  // auto ii = vectorNameDimMap.find(name);
  // if (ii!=vectorNameDimMap.end())
  //   return ii->second;
  // else
  //   return 1;

  size_t found = name.find_last_of("_");
  if (found==std::string::npos)
    return 1;

  std::string dim = name.substr(found+1);
  if (utl::IsInt(dim))
    return utl::ConvertStringToNumber<int>(dim);
  else
    return 1;
}

/** Get scalar vector by its name. The name-dimention map is in vectorNameDimMap */
template <class T>
inline std::vector<T>
GetScalarsByName(const std::vector<T>& vec, const std::vector<std::string>& nameVec, const std::string& name)
{
  int start=0;
  std::vector<T> result;
  for ( int i = 0; i < nameVec.size(); ++i ) 
    {
    int len = GetScalarsDimentionByName(nameVec[i]);
    if (nameVec[i]==name)
      {
      for ( int j = 0; j < len; ++j ) 
        result.push_back(vec[start+j]);
      return result;
      }
    start += len;
    }
  return result;
}

template <class T>
inline void
SetScalarsByName(const std::vector<T>& vec, const std::vector<std::string>& nameVec, const std::vector<T>& scalars, const std::string& name)
{
  int start=0;
  for ( int i = 0; i < nameVec.size(); ++i ) 
    {
    int len = GetScalarsDimentionByName(nameVec[i]);
    utlSAException(len!=scalars.size())(nameVec[i])(len)(scalars.size()).msg("the size of scalars is not consistent");
    if (nameVec[i]==name)
      {
      for ( int j = 0; j < len; ++j ) 
        vec[start+j] = scalars[j];
      }
    start += len;
    }
}

template <class T>
inline void
RemoveScalarsByName(std::vector<T>& vec, const std::vector<std::string>& nameVec, const std::string& name)
{
  int start=0, len=0, i;
  for ( i = 0; i < nameVec.size(); ++i ) 
    {
    len = GetScalarsDimentionByName(nameVec[i]);
    if (nameVec[i]==name)
      break;
    start += len;
    }

  // no name found
  if (i==nameVec.size())
    return;

  // // remove name
  // int jj = i;
  // nameVec.erase(nameVec.begin()+jj);

  // remove scalars
  for ( i = start; i+len < vec.size(); ++i ) 
    vec[i] = vec[i+len];

  vec.resize(vec.size()-len);
}
 

/** get the maximal rank of spherical harmonic coefficient vector with a given dimension  */
inline int
DimToRankSH ( const int dimm ) 
{
  double result = -1.5 + std::sqrt(0.25 + 2.0*dimm);
  
  utlGlobalAssert(IsInt(result) && IsEven(result), "Logical ERROR! Rank of SH here must be an even integer. Wrong dimension! rank=" << result <<", dimm = " << dimm);
  
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
  else if (v6dStoreWay==TENSOR_DIAGONAL_FIRST)
    {
    v9d[0] = v6d[0]; // xx
    v9d[1] = v6d[3]; // xy
    v9d[2] = v6d[4]; // xz
    v9d[3] = v6d[3]; // yx
    v9d[4] = v6d[1]; // yy
    v9d[5] = v6d[5]; // yz
    v9d[6] = v6d[4]; // zx
    v9d[7] = v6d[5]; // zy
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
  else if (v6dStoreWay==TENSOR_DIAGONAL_FIRST)
    {
    v6d[0] = v9d[0]; // xx
    v6d[1] = v9d[4]; // yy
    v6d[2] = v9d[8]; // zz
    v6d[3] = v9d[1]; // xy
    v6d[4] = v9d[2]; // xz
    v6d[5] = v9d[5]; // yz
    }
  else
    utlException(true, "wrong v1StoreWay");
}

template<typename V1Type, typename V2Type>
void 
ConvertTensor6DTo6D(const V1Type& v6d1, V2Type& v6d2, int s1, int s2)
{
  if (s1==s2)
    {
    for ( int i = 0; i < 6; ++i ) 
      v6d2[i] = v6d1[i];
    }
  else
    {
    std::vector<double> v9d(9);
    ConvertTensor6DTo9D(v6d1, v9d, s1);
    ConvertTensor9DTo6D(v9d, v6d2, s2);
    }
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
ReadGrad(const std::string& grad_str, const int NoSymmetricDuple=DIRECTION_NODUPLICATE, const int mode=CARTESIAN_TO_SPHERICAL, 
  const int flipx=DIRECTION_NOFLIP, const int flipy=DIRECTION_NOFLIP, const int flipz=DIRECTION_NOFLIP, const bool need_normalize=true) 
{
  std::vector < std::vector<std::string> > strVec;

  // use three delimiters ' ', ',', '\t'  
  ReadLinesFirstlineCheck(grad_str, strVec, " \t,");

  for ( int i = 0; i < strVec.size(); i += 1 ) 
    utlGlobalException(strVec[i].size()!=3, "wrong dimension in gradient file! strVec[i].size()="<< strVec[i].size() <<", grad_str="<<grad_str);

  utl_shared_ptr<NDArray<T,2> > grad (new NDArray<T,2>() );
  if (NoSymmetricDuple==DIRECTION_DUPLICATE)
    grad->ReSize(2*strVec.size(),3);
  else
    grad->ReSize(strVec.size(),3);

  std::vector<double> xyz(3);
  int jj=0;
  for (int i=0; i<strVec.size(); i++) 
    {
    utlGlobalException(strVec[i].size()!=3, "wrong size for grad");
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

inline int 
GetFiberTractsFormatFromFileExtension(const std::string& filename)
{
  std::string ext, fileNoExt;
  utl::GetFileExtension(filename, ext, fileNoExt);
  int format = TRACTS_UNKNOWN;
  if (ext=="trk")
    format = TRACTS_TRK;
  else if (ext=="tck")
    format = TRACTS_TCK;
  else if (ext=="vtk")
    format = TRACTS_VTK;
  else
    format = TRACTS_UNKNOWN;
  return format;
}

inline int 
GetFiberTractsFormat(const std::string& filename)
{
  std::string ext, fileNoExt;
  utl::GetFileExtension(filename, ext, fileNoExt);

  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (!file)
    utlGlobalException(true, "Unable to open file " + filename);

  int format = TRACTS_UNKNOWN;
  if (ext=="trk")
    {
    format = TRACTS_TRK;
    char code[6];
    fseek (file, 0, SEEK_SET);  
    size_t rsize = fread(code, 1, 6, file);
    utlGlobalException(strcmp(code,"TRACK")!=0, "Wrong data format TRACTS_TRK. Magic code is wrong. It is " + std::string(code) +". Should be 'TRACK'");
    utlGlobalException(rsize!=6, "wrong size read");
    }
  else if (ext=="tck")
    {
    format = TRACTS_TCK;
    }
  else if (ext=="vtk")
    {
    format = TRACTS_VTK;
    }
  else
    format = TRACTS_UNKNOWN;
  fclose(file);
  return format;
}

/** 
 * Get design matrix for DTI from gradient matrix and b values.  
 *
 * \param gradMat Nx3 gradient matrix in cartesian format. 
 * \param dwi_normalize whether DWI samples are normalized using the b0 image. 
 *                      For normalized DWI samples, b0 is not needed to re-estimated. 
 * 
 * \return Nx6 or Nx7 design matrix for DTI. First 6 values are for 
 * \f$ -b*[g_x g_x, 2.0*g_x g_y, 2.0*g_x g_z, g_y g_y, 2.0*g_y g_z, g_z g_z] \f$.
 * */
template <class T>
utl::Matrix<T> 
GetDTIDesignMatrix ( const utl::Matrix<T>& gradMat, const std::vector<T>& bVec, int dwi_normalize )
{
  utlSAException(gradMat.Rows()!=bVec.size())(gradMat.Rows())(bVec.size()).msg("wrong size of gradMat and bVec");
  utl::Matrix<T> mat;
  if (dwi_normalize==DWI_NORMALIZE)
    mat.ReSize(gradMat.Rows(), 6);
  else if (dwi_normalize==DWI_NONORMALIZE)
    mat.ReSize(gradMat.Rows(), 7);
  else
    utlException(true, "wrong dwi_normalize");

  for ( int i = 0; i < gradMat.Rows(); ++i ) 
    {
    mat(i,0) = -1.0 * bVec[i] * gradMat(i,0) * gradMat(i,0);  // xx
    mat(i,1) = -2.0 * bVec[i] * gradMat(i,0) * gradMat(i,1);  // xy
    mat(i,2) = -2.0 * bVec[i] * gradMat(i,0) * gradMat(i,2);  // xz
    mat(i,3) = -1.0 * bVec[i] * gradMat(i,1) * gradMat(i,1);  // yy
    mat(i,4) = -2.0 * bVec[i] * gradMat(i,1) * gradMat(i,2);  // yz
    mat(i,5) = -1.0 * bVec[i] * gradMat(i,2) * gradMat(i,2);  // zz
    if (dwi_normalize==DWI_NONORMALIZE)
      mat(i,6) = 1;                                           // S0
    }
  return mat;
}

/* Function taken from 3D Slicer, vtkDiffusionTensorMathematics
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
