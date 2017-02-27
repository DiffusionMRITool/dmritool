/**
 *       @file  utlITKConceptChecking.h
 *      @brief  
 *     Created  "06-10-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlITKConceptChecking_h
#define __utlITKConceptChecking_h

#include <itkConceptChecking.h>

namespace itk
{

namespace Detail
{
}

namespace Concept
{

template< int D1, int D2 >
struct SameInteger 
  {
  struct Constraints 
    {
    typedef Detail::UniqueType_int< D1 > DT1;
    typedef Detail::UniqueType_int< D2 > DT2;
    void constraints()
      {
      DT1 a = DT2();

      Detail::IgnoreUnusedVariable(a);
      }
    };
  itkConceptConstraintsMacro();
  };


}

}


#endif 
