/**
 *       @file  utl4thOrderTensor.h
 *      @brief  
 *     Created  "06-21-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */
#ifndef __utl4thOrderTensor_h
#define __utl4thOrderTensor_h

#include "utlNDArray.h"

namespace utl
{


template <class T>
inline void 
Convert1To2Tensor(const utl::NDArray<T,1>& vec, utl::NDArray<T,2>& mat);


/**
 *   \class   NDArray<T,4>
 *   \brief   NDArray<T,4> is a 4th order tensor
 *
 *   \author  Jian Cheng
 *   \date    05-22-2016
 *   \ingroup utlNDArray Math
 */
template < class T >
class NDArray<T, 4> : public NDArrayBase<T,4>
{
public:
  typedef NDArray                  Self;
  typedef NDArrayBase<T,4>            Superclass;

  typedef typename Superclass::ValueType         ValueType;
  typedef typename Superclass::ScalarValueType    ScalarValueType;

  typedef typename Superclass::SizeType          SizeType;
  typedef typename Superclass::ShapeType         ShapeType;
  typedef typename Superclass::Pointer           Pointer;
  typedef typename Superclass::ConstPointer      ConstPointer;
  typedef typename Superclass::Reference         Reference;
  typedef typename Superclass::ConstReference    ConstReference;
  typedef typename Superclass::Iterator          Iterator;
  typedef typename Superclass::ConstIterator     ConstIterator;

  using Superclass::Dimension;
  using Superclass::SubDimension;
  
  using Superclass::SetData;
  using Superclass::CopyData;
  using Superclass::ReSize;
  using Superclass::operator();
  using Superclass::operator=;
  using Superclass::operator+=;
  using Superclass::operator-=;
  using Superclass::operator%=;
  using Superclass::operator/=;

  NDArray() : Superclass() 
    { }
  NDArray(const NDArray<T,4>& mat) : Superclass(mat)
    { }
  template<typename EType>
  NDArray(const Expr<EType, typename EType::ValueType>& expr) : Superclass(expr)
    { }

  template< typename TMatrixValueType >
  NDArray(const NDArray< TMatrixValueType,4> & r) : Superclass(r)
    {  }
  
  NDArray(NDArray<T,4>&& r) 
  { 
    operator=(std::move(r)); 
  }

  explicit NDArray(const ShapeType& shape) : Superclass(shape) { }
  
  NDArray(const T* vec, const SizeType s0, const SizeType s1, const SizeType s2, const SizeType s3) : Superclass() 
    { 
    SizeType shape[4];
    shape[0]=s0, shape[1]=s1, shape[2]=s2, shape[3]=s3;
    __utl_ndarray_alloc_blah(shape);
    utl::cblas_copy<T>(this->Size(), vec, 1, this->m_Data, 1);
    }
  
  NDArray(const SizeType s0, const SizeType s1, const SizeType s2, const SizeType s3) : Superclass()
    { 
    SizeType shape[4];
    shape[0]=s0, shape[1]=s1, shape[2]=s2, shape[3]=s3;
    __utl_ndarray_alloc_blah(shape);
    }
  
  NDArray(const SizeType s0, const SizeType s1, const SizeType s2, const SizeType s3, const T r) : Superclass()
    { 
    SizeType shape[4];
    shape[0]=s0, shape[1]=s1, shape[2]=s2, shape[3]=s3;
    __utl_ndarray_alloc_blah(shape);
    this->Fill(r);
    }

/**
 * Constructor assumes input points to array of correct size.
 * Values are copied individually instead of with a binary copy.  This
 * allows the T's assignment operator to be executed.
 */
  NDArray(const T* vec, const ShapeType& shape) : Superclass(vec,shape)   {  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArray(const ShapeType& shape, const T r) : Superclass(shape, r)   {  }

  NDArray<T,4>& operator=(NDArray<T,4> & r)
    {
    Superclass::operator=(r);
    return *this;
    }

  NDArray<T,4>& operator=(NDArray<T,4> && r)
    {
    if ( this != &r ) 
      {
      this->Clear();
      this->Swap(r);
      }
    return *this;
    }

  // 4th order tensor related
  
  UTL_ALWAYS_INLINE Reference operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
  { return this->m_Data[i*this->m_OffSetTable[1] + j*this->m_OffSetTable[2] + k*this->m_OffSetTable[3]+l];}
  UTL_ALWAYS_INLINE ConstReference operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
  { return this->m_Data[i*this->m_OffSetTable[1] + j*this->m_OffSetTable[2] + k*this->m_OffSetTable[3]+l];}
  
  inline bool ReSize(const SizeType s0, const SizeType s1, const SizeType s2, const SizeType s3)
    {
    SizeType shape[4];
    shape[0]=s0, shape[1]=s1, shape[2]=s2, shape[3]=s3;
    return Superclass::ReSize(shape);
    }

  NDArray<T,2> GetRefSubMatrix(const int i, const int j) const
    {
    NDArray<T,2> mat;
    unsigned index[4];
    index[0]=i, index[1]=j, index[2]=0, index[3]=0;
    unsigned shape[2];
    shape[0]=this->m_Shape[2], shape[1]=this->m_Shape[3];
    mat.SetData(this->m_Data+this->GetOffset(index), shape);
    return mat;
    }

  /** Make it major symmetric  */
  void MajorSymmetrize()
    {
    int size = this->m_Shape[0];
    utlException(this->m_Shape[1]!=size || this->m_Shape[2]!=size || this->m_Shape[3]!=size, "should have the same size in every dimension");
    SizeType index1[Dimension], index2[Dimension];

    NDArray<int, Dimension> isFilled(this->m_Shape, 0);
    for ( int i = 0; i < size; ++i ) 
      for ( int j = 0; j < size; ++j ) 
        for ( int k = 0; k < size; ++k ) 
          for ( int l = 0; l < size; ++l ) 
            {
            // utlPrintVar4(true, i,j,k,l);
            index1[0]=i, index1[1]=j, index1[2]=k, index1[3]=l;
            index2[0]=k, index2[1]=l, index2[2]=i, index2[3]=j;
            if (isFilled(index1)==0 || isFilled(index2)==0)
              {
              T val = ((*this)(index1) + (*this)(index2))/2.0;
              (*this)(index1)=val;
              (*this)(index2)=val;

              isFilled(index1)=1;
              isFilled(index2)=1;
              }
            }
    }

  /** Test whether it is major symmetric  */
  bool IsMajorSymmetric() const
    {
    int size = this->m_Shape[0];
    utlException(this->m_Shape[1]!=size || this->m_Shape[2]!=size || this->m_Shape[3]!=size, "should have the same size in every dimension");
    SizeType index1[Dimension], index2[Dimension];
    NDArray<int, Dimension> isFilled(this->m_Shape, 0);
    for ( int i = 0; i < size; ++i ) 
      for ( int j = 0; j < size; ++j ) 
        for ( int k = 0; k < size; ++k ) 
          for ( int l = 0; l < size; ++l ) 
            {
            index1[0]=i, index1[1]=j, index1[2]=k, index1[3]=l;
            index2[0]=k, index2[1]=l, index2[2]=i, index2[3]=j;
            if (isFilled(index1)==0 || isFilled(index2)==0)
              {
              T val1 = (*this)(index1);
              T val2 = (*this)(index2);
              if (!utl::IsSame(val1, val2, 1e-10) )
                {
                // utl::PrintVector(index1, 4, "index1");
                // utl::PrintVector(index2, 4, "index2");
                // utlPrintVar2(true, (*this)(index1), (*this)(index2));
                return false;
                }
              isFilled(index1)=1;
              isFilled(index2)=1;
              }
            }
    return true;
    }
  
  void MinorSymmetrize()
    {
    int size = this->m_Shape[0];
    utlException(this->m_Shape[1]!=size || this->m_Shape[2]!=size || this->m_Shape[3]!=size, "should have the same size in every dimension");
    SizeType index1[Dimension], index2[Dimension], index3[Dimension], index4[Dimension];

    NDArray<int, Dimension> isFilled(this->m_Shape, 0);
    for ( int i = 0; i < size; ++i ) 
      for ( int j = 0; j < size; ++j ) 
        for ( int k = 0; k < size; ++k ) 
          for ( int l = 0; l < size; ++l ) 
            {
            // utlPrintVar4(true, i,j,k,l);
            index1[0]=i, index1[1]=j, index1[2]=k, index1[3]=l;
            index2[0]=j, index2[1]=i, index2[2]=k, index2[3]=l;
            index3[0]=i, index3[1]=j, index3[2]=l, index3[3]=k;
            index4[0]=j, index4[1]=i, index4[2]=l, index4[3]=k;

            if (isFilled(index1)==0 || isFilled(index2)==0 || isFilled(index3)==0 || isFilled(index4)==0)
              {
              T val1 = (*this)(index1),
                val2 = (*this)(index2),
                val3 = (*this)(index3),
                val4 = (*this)(index4);
              (*this)(index1) = (val1+val2+val3+val4)/4.0;
              (*this)(index2) = (val1+val2+val3+val4)/4.0;
              (*this)(index3) = (val1+val2+val3+val4)/4.0;
              (*this)(index4) = (val1+val2+val3+val4)/4.0;

              isFilled(index1)=1;
              isFilled(index2)=1;
              isFilled(index3)=1;
              isFilled(index4)=1;
              }
            }
    }

  bool IsMinorSymmetric() const
    {
    int size = this->m_Shape[0];
    utlException(this->m_Shape[1]!=size || this->m_Shape[2]!=size || this->m_Shape[3]!=size, "should have the same size in every dimension");

    NDArray<int, Dimension> isFilled(this->m_Shape, 0);
    SizeType index1[Dimension], index2[Dimension], index3[Dimension], index4[Dimension];
    for ( int i = 0; i < size; ++i ) 
      for ( int j = 0; j < size; ++j ) 
        for ( int k = 0; k < size; ++k ) 
          for ( int l = 0; l < size; ++l ) 
            {
            index1[0]=i, index1[1]=j, index1[2]=k, index1[3]=l;
            index2[0]=j, index2[1]=i, index2[2]=k, index2[3]=l;
            index3[0]=i, index3[1]=j, index3[2]=l, index3[3]=k;
            index4[0]=j, index4[1]=i, index4[2]=l, index4[3]=k;

            if (isFilled(index1)==0 || isFilled(index2)==0 || isFilled(index3)==0 || isFilled(index4)==0)
              {
              T val1 = (*this)(index1),
                val2 = (*this)(index2),
                val3 = (*this)(index3),
                val4 = (*this)(index4);
              if (!utl::IsSame(val1, val2, 1e-10) )
                return false;
              if (!utl::IsSame(val1, val3, 1e-10) )
                return false;
              if (!utl::IsSame(val1, val4, 1e-10) )
                return false;
              isFilled(index1)=1;
              isFilled(index2)=1;
              isFilled(index3)=1;
              isFilled(index4)=1;
              }
            }
    return true;
    }
  
  void EigenDecompositionWithMinorSymmetry (NDArray<T,1>& valReal, NDArray<T,1>& valImg) const
    {
    utlException(this->m_Shape[0]!=3 || this->m_Shape[1]!=3 || this->m_Shape[2]!=3 || this->m_Shape[3]!=3, "should have size (3,3,3,3)");
    NDArray<T, 2> mat; 
    utl::Convert4To2Tensor(*this, mat);
    mat.EigenDecompositionNonSymmetricMatrix(valReal, valImg);
    }

  void EigenDecompositionWithMinorSymmetry (NDArray<T,1>& valReal, NDArray<T,1>& valImg, std::vector<NDArray<T,2> >& matRealR, std::vector<NDArray<T,2> >& matImgR) const
    {
    utlException(this->m_Shape[0]!=3 || this->m_Shape[1]!=3 || this->m_Shape[2]!=3 || this->m_Shape[3]!=3, "should have size (3,3,3,3)");
    matRealR.resize(6); matImgR.resize(6);
    NDArray<T, 2> mat, vecRealR, vecImgR; 
    utl::Convert4To2Tensor(*this, mat);
    mat.EigenDecompositionNonSymmetricMatrix(valReal, valImg, vecRealR, vecImgR);
    NDArray<T,1> vecTmp(6);
    for ( int i = 0; i < valReal.Size(); ++i ) 
      {
      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecRealR(j,i);
      utl::Convert1To2Tensor(vecTmp, matRealR[i]);
      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecImgR(j,i);
      utl::Convert1To2Tensor(vecTmp, matImgR[i]);
      }
    }

  void EigenDecompositionWithMinorSymmetry (NDArray<T,1>& valReal, NDArray<T,1>& valImg, std::vector<NDArray<T,2> >& matRealR, std::vector<NDArray<T,2> >& matImgR, std::vector<NDArray<T,2> >& matRealL, std::vector<NDArray<T,2> >& matImgL) const
    {
    utlException(this->m_Shape[0]!=3 || this->m_Shape[1]!=3 || this->m_Shape[2]!=3 || this->m_Shape[3]!=3, "should have size (3,3,3,3)");
    matRealR.resize(6); matImgR.resize(6);
    matRealL.resize(6); matImgL.resize(6);
    NDArray<T, 2> mat, vecRealR, vecImgR, vecRealL, vecImgL; 
    utl::Convert4To2Tensor(*this, mat);
    mat.EigenDecompositionNonSymmetricMatrix(valReal, valImg, vecRealR, vecImgR, vecRealL, vecImgL);
    NDArray<T,1> vecTmp(6);
    for ( int i = 0; i < valReal.Size(); ++i ) 
      {
      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecRealR(j,i);
      utl::Convert1To2Tensor(vecTmp, matRealR[i]);
      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecImgR(j,i);
      utl::Convert1To2Tensor(vecTmp, matImgR[i]);

      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecRealL(j,i);
      utl::Convert1To2Tensor(vecTmp, matRealL[i]);
      for ( int j = 0; j < 6; ++j ) 
        vecTmp[j] = vecImgL(j,i);
      utl::Convert1To2Tensor(vecTmp, matImgL[i]);
      }
    }

};


}


#endif 
