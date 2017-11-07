/*=========================================================================

 Program:   Command Progress Update

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkCommandProgressUpdate_h
#define __itkCommandProgressUpdate_h

#include "itkCommand.h"
#include "itkObject.h"
#include "itkProcessObject.h"
#include "itkTimeProbe.h"
#include "itkRealTimeClock.h"
#include <iostream>
#include <iomanip>

namespace itk
{

/** \brief Observer for progress notification.
 *
 *  \author Pew-Thian Yap, UNC Chapel Hill, ptyap@med.unc.edu
 */

class CommandProgressUpdate : public Command
{
public:
  typedef CommandProgressUpdate   Self;
  typedef Command                 Superclass;
  typedef SmartPointer<Self>      Pointer;
  itkNewMacro(Self);
  
  typedef RealTimeClock ClockType;
  typedef ClockType::TimeStampType TimeStampType;
  
  itkSetMacro(UseColor, bool);
  itkGetConstMacro(UseColor, bool);
  itkBooleanMacro(UseColor);

  typedef enum {LIST=0, BAR} StyleType;
  const static unsigned int NumberOfBarSegments = 20;
  
  itkSetMacro(Style, StyleType);
  itkGetConstMacro(Style, StyleType);
  
protected:
  CommandProgressUpdate()
    {
    this->UseColorOn();
    m_Style = BAR;
    m_SuppressOutput = false; // avoid repeating output for 100%
    m_Clock = ClockType::New();
    m_StartTimeStamp = m_Clock->GetTimeInSeconds();
    m_LastProgress = 0;
    }

private:
  TimeProbe clock;
  bool m_UseColor;
  StyleType m_Style;
  bool m_SuppressOutput;
  ClockType::Pointer m_Clock;
  TimeStampType m_StartTimeStamp;
  float m_LastProgress;
  
public:
  void Execute(Object *caller, const EventObject & event) ITK_OVERRIDE
    {
    Execute( (const Object *)caller, event );
    }
  
  void Execute(const Object * object, const EventObject & event) ITK_OVERRIDE
    {
    const ProcessObject * filter =
      dynamic_cast< const ProcessObject * >(object);
    if (!ProgressEvent().CheckEvent( &event ))
      {
      return;
      }
    // In case filter needs to be executed multiple times
    if (m_SuppressOutput)
      {
      m_SuppressOutput = false;
      return;
      }
    float progress = filter->GetProgress();
    unsigned int progressInRoundedPercent =
      static_cast<unsigned int>( 100 * progress );
    
    std::string redPrefix = "";
    std::string greenPrefix = "";
    std::string yellowPrefix = "";
    std::string bluePrefix = "";
    std::string magentaPrefix = "";
    std::string cyanPrefix = "";
    std::string suffix = "";
    
    if ( m_UseColor )
      {
      redPrefix = "\033[1;31m";
      greenPrefix = "\033[32m";
      yellowPrefix = "\033[33m";
      bluePrefix = "\033[1;34m";
      magentaPrefix = "\033[35m";
      cyanPrefix = "\033[36m";
      suffix = "\033[0m";
      }
    
    switch ( m_Style )
      {
      case LIST:
        
        std::cout << "  >>>> "
        //      << std::setiosflags(std::ios::fixed)
        //      << std::setprecision(2)
          << yellowPrefix
          << std::setw(3)
          << std::setiosflags(std::ios::right)
          << progressInRoundedPercent << "%"
          << suffix << " completed" << std::flush;
        
        if ( progress == 0 )
          {
          m_SuppressOutput = false;
          std::cout << std::endl;
          clock.Start();
          }
        else
          {
          clock.Stop();
          std::cout
            << ", "
            << magentaPrefix
            << std::setiosflags(std::ios::fixed)
            << std::setprecision(3)
          //      << std::setw(5)
          //      << setiosflags(std::ios::right)
           << clock.GetMean() << "s" << suffix << " elapsed" << std::endl;
          std::cout.precision(6);
          clock.Start();
          }
        if (progress == 1) m_SuppressOutput = true;
        break;
        
      case BAR:
        if ( progress == 0 )
          {
          m_LastProgress = progress;
          m_SuppressOutput = false;
          m_StartTimeStamp = m_Clock->GetTimeInSeconds();
          clock.Start();
          }

        clock.Stop();

        std::cout
//          << "\e[?25l"
          << "\r"
          << yellowPrefix
          << std::setw(3)
          << std::setiosflags(std::ios::right)
          << progressInRoundedPercent << "%" << suffix << std::flush;
        
        std::cout << " [" << std::flush;
        for (unsigned int k=0; k<NumberOfBarSegments+1; k++)
          {
          if ( k == (unsigned int)(progress * NumberOfBarSegments) )
            {
            std::cout << greenPrefix << ">" << suffix << std::flush;
            }
          else if ( k < (unsigned int)(progress * NumberOfBarSegments) )
            {
            std::cout << greenPrefix << "=" << suffix << std::flush;
            }
          else
            {
            std::cout << " " << std::flush;
            }
          }
        
        if ( progress < 1 )
          {
          unsigned int secondsLeft = (unsigned int)
            (clock.GetMean() * (1 - progress) * 100 + 0.5);
          unsigned int hoursLeft = secondsLeft/3600;
          secondsLeft = secondsLeft % 3600;
          unsigned int minutesLeft = secondsLeft/60;
          secondsLeft = secondsLeft % 60;

          unsigned int secondsElapsed = (unsigned int)
            (m_Clock->GetTimeInSeconds() - m_StartTimeStamp + 0.5);
          unsigned int hoursElapsed = secondsElapsed/3600;
            secondsElapsed = secondsElapsed % 3600;
          unsigned int minutesElapsed = secondsElapsed/60;
            secondsElapsed = secondsElapsed % 60;
          
          float secondsPerPercent;
          
          if ( progress > 0.01 )
            {
            secondsPerPercent = clock.GetMean()/(progress - m_LastProgress)/100;
            }
          else
            {
            secondsPerPercent = 0;
            }
          
          m_LastProgress = progress;
          
          std::cout << "]  " << std::flush;
          std::cout
            << magentaPrefix
            << std::setiosflags(std::ios::fixed)
            << std::setprecision(3)
            << secondsPerPercent << " s" << "/%"
            << suffix
            << cyanPrefix
            << "  ETA "
            << std::setfill('0')
            << std::setw(2)
            << hoursLeft
            << ":"
            << std::setw(2)
            << minutesLeft
            << ":"
            << std::setw(2)
            << secondsLeft
            << std::setfill(' ')
            << suffix
            << bluePrefix
            << "  ET "
            << std::setfill('0')
            << std::setw(2)
            << hoursElapsed
            << ":"
            << std::setw(2)
            << minutesElapsed
            << ":"
            << std::setw(2)
            << secondsElapsed
            << std::setfill(' ')
            << suffix
            << " " << std::flush;
          std::cout.precision(6);
          clock.Start();
          }
        else
          {
          std::cout << "]"
//            << "\e[?25h"
            << std::endl;
          m_SuppressOutput = true;
          }
        
        break;
      }

    }
};

} // end namespace itk


#endif



