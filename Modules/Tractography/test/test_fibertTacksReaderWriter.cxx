/**
 *       @file  test_fibertTacksReaderWriter.cxx
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "itkFiberTractsReader.h"
#include "itkFiberTractsWriter.h"

#include "utlCommandLineParser.h"


int 
main (int argc, char const* argv[])
{
  std::string _fileIn  = utl_option("-i", "", "input file of tracks");
  std::string _fileOut  = utl_option("-o", "", "output file of tracks");

  utl_showdoc("-h"); 
  utl_showdoc("--help"); 

  itk::FiberTractsReader::Pointer reader = itk::FiberTractsReader::New();
  reader->SetFileName(_fileIn);

  reader->Update();

  itk::FiberTracts<>::Pointer fibers = reader->GetOutput();

  fibers->Print(std::cout<<"fibers=\n");

  if (_fileOut!="")
    {
    itk::FiberTractsWriter::Pointer writer = itk::FiberTractsWriter::New();
    writer->SetFileName(_fileOut);
    writer->SetFiberTracts(fibers);
    writer->Update();
    }
  
  return 0;
}
