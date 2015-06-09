/**
 *       @file  ComputeSHCoefficientsOfDWIFromSymmetricTensor.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-25-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "ComputeSHCoefficientsOfDWIFromSymmetricTensorCLP.h"
#include "itkSpecialFunctionGenerator.h"
#include "utlCore.h"

void 
DoMain ( std::vector<float>& e1e2, std::vector<float>& bvec, const int shOrder, std::ostream& os, const bool outputAll, const bool nob )
{
  for ( int i = 0; i < bvec.size(); i += 1 ) 
    {
    double b = bvec[i];
    std::vector<double> coef = utl::GetSymmetricTensorSHCoef<double>(b, e1e2[0], e1e2[1], shOrder, 0.0, 0.0 );
    if (!nob)
      os << b << " ";
    if (outputAll)
      {
      for ( int i = 0; i < coef.size(); i += 1 ) 
        os << coef[i] << " ";
      }
    else
      {
      for ( int l = 0; l <= shOrder; l += 2 ) 
        {
        int jj = utl::GetIndexSHj(l, 0);
        os << coef[jj] << " ";
        }
      }
    os << std::endl;
    }
}

/**
 * \brief  Compute SH coefficients of DWI signal in DTI with symmetric tensor  (along z-axis).
 *
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  utlGlobalException(_EigenValuesArg.isSet() && _FAMDArg.isSet(), "only one of --eigenvalues and --famd can be set");
  utlGlobalException(!_EigenValuesArg.isSet() && !_FAMDArg.isSet(), "need to set one of --eigenvalues and --famd");
  utlGlobalException(_BValues.size()==0, "need to set --bvalues");

  if (_EigenValuesArg.isSet())
    {
    utlGlobalException(_EigenValues.size()!=2, "wrong size in --eigenvalues");
    }
  else
    {
    utlGlobalException(_FAMD.size()!=2, "wrong size in --eigenvalues");
    _EigenValues = utl::GetE1E2FromFAMD(_FAMD[0], _FAMD[1]);
    }

  if (_OutputFileArg.isSet())
    {
    std::ofstream  out; 
    out.open ( _OutputFile.c_str() );
    DoMain(_EigenValues, _BValues, _SHOrder, out, _OutputAllCoefficientsArg.isSet(), _NoBArg.isSet());
    out.close();
    }
  else
    DoMain(_EigenValues, _BValues, _SHOrder, std::cout, _OutputAllCoefficientsArg.isSet(), _NoBArg.isSet());

  return 0;
}
