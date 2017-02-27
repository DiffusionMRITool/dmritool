/**
 *       @file  TextFileOperator.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "02-05-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlCore.h"
#include "TextFileOperatorCLP.h"


enum
{
  OP_NULL,
  OP_INFO, 
  OP_TRANSPOSE, 
  OP_EXTRACT,
  OP_CONNECT_ROW,
  OP_CONNECT_COLUMN,
  OP_INDEX_ROW,
  OP_SCALE
};

void
SetOperationWithChecking( int &operation, int value )
{
  if ( operation == OP_NULL )
    {
    operation = value;
    }
  else
    {
    std::cerr << "Only one type of operation is allowed!" << std::endl;
    exit( EXIT_FAILURE );
    }
}


/**
 * \brief  operators for text file
 */
int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;
  
  int operation = OP_NULL;

  if ( _InfoArg.isSet() )
    SetOperationWithChecking( operation, OP_INFO );
  if ( _TransposeArg.isSet() )
    SetOperationWithChecking( operation, OP_TRANSPOSE );
  if (_FileConnectRowArg.isSet())
    SetOperationWithChecking( operation, OP_CONNECT_ROW );
  if (_FileConnectColumnArg.isSet())
    SetOperationWithChecking( operation, OP_CONNECT_COLUMN );
  if ( _FileIndexRowsArg.isSet() )
    SetOperationWithChecking( operation, OP_INDEX_ROW );
  if ( _ExtractSubMatrixArg.isSet() )
    SetOperationWithChecking( operation, OP_EXTRACT );
  if ( _ScaleArg.isSet() )
    SetOperationWithChecking( operation, OP_SCALE );
  // if ( _ScaleColumnArg.isSet() )
  //   SetOperationWithChecking( operation, OP_SCALE_COLUMN );
  // if ( _ScaleRowArg.isSet() )
  //   SetOperationWithChecking( operation, OP_SCALE_ROW );

  std::vector<std::vector<std::string> > input, input2, output;
  utl::ReadLines(_InputFile, input, " \t,");

  int NumberRows = input.size();
  utlException(NumberRows==0, "no data");

  int NumberColumns=input[0].size();
  bool isSameColumn=true;
  for ( int i = 1; i < input.size(); i += 1 ) 
    {
    if (NumberColumns!=input[i].size())
      isSameColumn = false;
    if (NumberColumns<input[i].size())
      NumberColumns = input[i].size();
    }
  for ( int i = 0; i < input.size(); i += 1 ) 
    {
    int numberPush=NumberColumns-input[i].size();
    for ( int j = 0; j < numberPush; j += 1 ) 
      input[i].push_back(" ");
    }
  

  std::ofstream  out;   
  if (_OutputfileArg.isSet())
    {
    out.open ( _Outputfile.c_str() );           
    utlException (!out,  "ERROR : failed to open output file " << _Outputfile);
    }

  switch ( operation )
    {
  case OP_INFO :
      {
      std::cout << "File information" << std::endl << std::flush;
      std::cout << "Row: " << NumberRows << ", Column: " << NumberColumns << std::endl << std::flush;
      if (isSameColumn)
        std::cout << "all rows have the same number of columns" << std::endl << std::flush;
      else
        std::cout << "not all rows have the same number of columns" << std::endl << std::flush;
      break;
      }
  case OP_TRANSPOSE :
      {
      std::cout << "File transpose" << std::endl << std::flush;
      utlException (!_OutputfileArg.isSet(), "no output");
      for ( int i = 0; i < NumberColumns; i += 1 ) 
        {
        int NumberColumns = input[i].size();
        for ( int j = 0; j < NumberRows-1; j += 1 ) 
          out << input[j][i] <<  " ";
        out << input[NumberRows-1][i] << "\n";
        }
      break;
      }
  case OP_CONNECT_COLUMN :
      {
      std::cout << "Connect columns" << std::endl << std::flush;
      utlException (!_OutputfileArg.isSet(), "no output");
      utl::ReadLines(_FileConnectColumn, input2);
      for ( int i = 0; i < input.size(); i += 1 ) 
        {
        std::vector<std::string> element1 = input[i];
        if (i<input2.size())
          {
          std::vector<std::string> element2 = input2[i];
          for ( int j = 0; j < element2.size(); j += 1 ) 
            element1.push_back(element2[j]);
          }
        input[i] = element1;
        }
      utl::Save2DVector(input, out);
      break;
      }
  case OP_CONNECT_ROW :
      {
      std::cout << "Connect rows" << std::endl << std::flush;
      utlException (!_OutputfileArg.isSet(), "no output");
      utl::ReadLines(_FileConnectRow, input2);
      utl::Save2DVector(input, out);
      utl::Save2DVector(input2, out);
      break;
      }
  case OP_INDEX_ROW :
      {
      std::cout << "select rows" << std::endl << std::flush;
      utlException (!_OutputfileArg.isSet(), "no output");
      std::vector<int> index;
      utl::ReadVector(_FileIndexRows, index);
      input2 = utl::SelectVector(input, index);
      utl::Save2DVector(input2, out);
      break;
      }
  case OP_EXTRACT:
      {
      std::cout << "Extract sub-matrix" << std::endl << std::flush;
      _ExtractSubMatrix[0] = utl::max(_ExtractSubMatrix[0], 0);
      _ExtractSubMatrix[1] = utl::min(_ExtractSubMatrix[1], NumberRows-1);
      _ExtractSubMatrix[2] = utl::max(_ExtractSubMatrix[2], 0);
      _ExtractSubMatrix[3] = utl::min(_ExtractSubMatrix[3], NumberColumns-1);
      // utlPrintVar4(true, _ExtractSubMatrix[0], _ExtractSubMatrix[1], _ExtractSubMatrix[2], _ExtractSubMatrix[3]);
      utlException (!_OutputfileArg.isSet(), "no output");
      for ( int i = _ExtractSubMatrix[0]; i <= _ExtractSubMatrix[1]; i += 1 ) 
        {
        std::vector<std::string> tmp;
        for ( int j = _ExtractSubMatrix[2]; j <= _ExtractSubMatrix[3]; j += 1 ) 
          tmp.push_back(input[i][j]);
        output.push_back(tmp);
        }
      utl::Save2DVector(output, out);
      break;
      }
  case OP_SCALE:
      {
      utlException (!_OutputfileArg.isSet(), "no output");
      for ( int i = 0; i < NumberRows; i += 1 ) 
        {
        std::vector<std::string> tmp;
        for ( int j = 0; j < NumberColumns; j += 1 ) 
          {
          std::string strTemp = input[i][j];
          if (utl::IsNumber(strTemp))
            {
            double floatTemp = utl::ConvertStringToNumber<double>(strTemp);
            floatTemp *= _Scale;
            strTemp = utl::ConvertNumberToString(floatTemp);
            }
          tmp.push_back(strTemp);
          }
        output.push_back(tmp);
        }
      utl::Save2DVector(output, out);
      break;
      }
  default :
    utlException(true, "wrong operation");
    break;
    }

  if (_OutputfileArg.isSet())
    {
    std::cout << "save to " << _Outputfile << std::endl << std::flush;
    out.close();
    }

  return 0;
}
