/**
 *       @file  gtest_test.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-19-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include <algorithm>
#include "gtest/gtest.h"
#include <cmath>

TEST(gtest_test, build)
{
  double a=0.1+0.2, b=0.3;
  std::cout << "(0.1+0.2)-0.3 = " << a-b << std::endl << std::flush;
  std::cout << "(0.1+0.2) == 0.3 = " << (a==b) << std::endl << std::flush;
  EXPECT_LE(std::fabs(0.1+0.2-0.3), 1e-10);
  EXPECT_EQ(5*5, 25+1e-15);
}

TEST(gtest_test, build2)
{
  EXPECT_EQ(5+5, 10);
}

