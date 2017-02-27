/**
 *       @file  utlDMRIGTest.cxx
 *      @brief  
 *     Created  "03-01-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "utlDMRI.h"
#include "itkSphericalHarmonicsGenerator.h"
#include "utlDMRIStoredTables.h"

TEST(utlDMRI, GetIndexSHlm)
{
  std::vector<int> lmVec(2);
  int j=0;
  for ( int l = 0; l <= 20; l += 2 ) 
    {
    for ( int m = -l; m <= l; m += 1 ) 
      {
      lmVec = utl::GetIndexSHlm(j);
      EXPECT_EQ(lmVec[0], l);
      EXPECT_EQ(lmVec[1], m);
      int jj = utl::GetIndexSHj(l,m);
      EXPECT_EQ(jj, j);
      j++;
      }
    }
}

TEST(utlDMRIStoredTables, ReadGrad)
{
  utl::ReadGrad<double>(3);
}

TEST(utlDMRIStoredTables_DeathTest, ReadGrad)
{
  EXPECT_DEATH(utl::ReadGrad<double>(8), "");
}
