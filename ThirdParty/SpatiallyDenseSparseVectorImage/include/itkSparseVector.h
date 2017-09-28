/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkSparseVector_h
#define __itkSparseVector_h

#include "itkIntTypes.h"

#include "utlSTDHeaders.h"

namespace itk
{
/** \class SparseVector
 * \brief Represents a sparse array.
 *
 * \author Pew-Thian Yap (ptyap@med.unc.edu)
 */
template< typename TValueType, typename TKeyType = SizeValueType >
class SparseVector
{
public:

  /** The element type stored at each location in the Array. */
  typedef TValueType                               ValueType;
  typedef TKeyType                                 KeyType;
  typedef TKeyType                                 ElementIdentifier;
  typedef TValueType                               ComponentType;
  typedef utl_unordered_map<TKeyType, TValueType> InternalDataType;
  typedef SparseVector                             Self;

  /** Default constructor */
  SparseVector()
  {
    m_FillBufferValue = 0;
    m_Data.clear();
  }

  /** Destructor */
  ~SparseVector()
  {
  }

  /** Set the all the elements of the array to the specified value */
  void Fill(TValueType const & v)
  {
    m_FillBufferValue = v;
    m_Data.clear();
  }

  /** Clear data */
  void Clear()
  {
    m_Data.clear();
  }

  /** Return modifiable lvalue reference to the element at specified index. No range checking. */
  TValueType & operator[](unsigned int i)
  {
    return m_Data[i];
  }

  /** Return constant rvalue reference to the element at specified index. No range checking. */
  TValueType operator[](unsigned int i) const
  {
    TValueType value;
    typename InternalDataType::const_iterator it = m_Data.find( i );

    if ( it == m_Data.end() )
      {
      value = m_FillBufferValue;
      }
    else
      {
      value = it->second;
      }

    return value;
  }

  inline unsigned int GetSize(void) const
  {
    return m_Data.size();
  }

  inline unsigned int Size(void) const
  {
    return m_Data.size();
  }

  /** Pointer to internal data */
  InternalDataType * GetDataPointer()
  {
    return &m_Data;
  }

  const InternalDataType * GetDataPointer() const
  {
    return &m_Data;
  }

private:

  InternalDataType  m_Data;                 // data
  ValueType         m_FillBufferValue;
};

} // namespace itk

#endif
