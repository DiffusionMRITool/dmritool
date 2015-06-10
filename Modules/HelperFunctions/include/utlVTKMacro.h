/**
 *       @file  utlVTKMacro.h
 *      @brief  
 *     Created  "02-18-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlVTKMacro_h
#define __utlVTKMacro_h

/** @addtogroup utlHelperFunctions
@{ */

#if VTK_MAJOR_VERSION <= 5
  #define vtkSetInputData(x,y) do { (x)->SetInput(y); } while (0) 
  #define vtkAddInputData(x,y) do { (x)->AddInput(y); } while (0)
#else
  #define vtkSetInputData(x,y) do { (x)->SetInputData(y); } while(0)
  #define vtkAddInputData(x,y) do { (x)->AddInputData(y); } while(0)
#endif

    /** @} */

#endif 
