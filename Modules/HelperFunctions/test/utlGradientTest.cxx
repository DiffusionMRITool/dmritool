/**
 *       @file  utlGradientTest.cxx
 *      @brief  
 *     Created  "10-03-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "utlCore.h"
#include "itkSphericalHarmonicsGenerator.h"
#include "utlDMRIStoredTables.h"
#include "utlCommandLineParser.h"

/**
 * \brief  test gradient table read
 */
int 
main (int argc, char const* argv[])
{
  utl_usage("test gradient table read");
  int tess = utl_option("-t", 3, "tess order");
  bool duple = utl_option("-d", false, "DUPLICATE for directions");
  bool flipx = utl_option("-flipx", false, "flipx");
  bool flipy = utl_option("-flipy", false, "flipy");
  bool flipz = utl_option("-flipz", false, "flipz");
  int mode = utl_option("-mode", 0, "0: CARTESIAN_TO_CARTESIAN; 1: CARTESIAN_TO_SPHERICAL; 2: SPHERICAL_TO_CARTESIAN; 3: SPHERICAL_TO_SPHERICAL");
  
  if (utl_option("-h",(const char *)NULL,0)) return 0;

    {
    utl::GradientTable<double>::Initialize(tess);
    utl_shared_ptr< utl::Matrix<double> > grad = utl::GradientTable<double>::GetGrad(tess, duple, mode, flipx, flipy, flipz);
    utl::PrintUtlMatrix(*grad, "grad");
    }
    
    {
    utl_shared_ptr< utl::Matrix<float> > grad1 = utl::ReadGrad<float>(tess, duple, mode, flipx, flipy, flipz);
    utl::PrintUtlMatrix(*grad1, "grad");
    }
  
  return 0;
}

