/**
 *       @file  utlCommandLineParser.h
 *      @brief  small but powerful CMD parser, borrowed from CImg, http://cimg.sourceforge.net/ 
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-28-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlCommandLineParser_h
#define __utlCommandLineParser_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include "utlCoreMacro.h"

#if UTL_OS==1
#include <unistd.h>
#endif


/** 
 * small but powerful CMD parser, borrowed from CImg, http://cimg.sourceforge.net/ 
 * */

/**
 * \code
 * utl_usage("bla bla ...");  
 * \endcode
 * */
#define utl_usage(usage) utl::utlOption((char*)NULL,(unsigned int)argc,(char**)argv,(char*)NULL,(char*)usage)
/**
 * \code
 * utl_help("bla bla ...");  
 * \endcode
 * */
#define utl_help(str) utl::utlOption((char*)NULL,(unsigned int)argc,(char**)argv,str,(char*)NULL)

/** \code 
 * double p1 = utl_option("-i1", 2.0, "bla bla...");
 * int p2   = utl_option("-i2", 3, "bla bla...");
 * std::string p3  = utl_option("-i3", "default", "bla bla...");
 * bool p4  = utl_option("-i4", false, "bla bla...");
 * \endcode
 * */
#define utl_option(name,defaut,usage) utl::utlOption((char*)name,(unsigned int)argc,(char**)argv,defaut,(char*)usage)

/** \code
 * utl_showdoc("-h"); 
 * utl_showdoc("--help"); 
 * \endcode
 * */
#define utl_showdoc(argu) if (utl_option(argu,(const char *)NULL,0)) return 0; 

namespace utl
{


inline const char* 
utlOption(const char *const name,const unsigned int argc,char **argv,const char *const defaut, const char *const usage=NULL) 
{
  static bool first=true, visu=false;
  const char *res = NULL;
  if (first) 
    { 
    first=false; 
    visu = utlOption("-h",argc,argv,(const char*)NULL)!=NULL;  
    visu |= utlOption("-help",argc,argv,(char*)NULL)!=NULL;
    visu |= utlOption("--help",argc,argv,(char*)NULL)!=NULL;
    if (argc==1) 
      visu=true;
    }
  if (!name && visu) 
    {
    if (usage)
      {
      std::string argv0Str(argv[0]);
      unsigned found = argv0Str.find_last_of("/\\");
      std::fprintf(stderr,"\n %s",utl::GetColoredString(argv0Str.substr(0,found), COLOR_RED).c_str());
      std::fprintf(stderr," : %s",usage);
      std::fprintf(stderr," (%s, %s)\n\n",__DATE__,__TIME__);
      }
    if (defaut)
      std::fprintf(stderr,"%s\n",defaut);
    }
  if (name) 
    {
    if (argc>0) 
      {
      unsigned int k=0;
      while (k<argc && strcmp(argv[k],name)) 
        k++;
      res=(k++==argc?defaut:(k==argc?argv[--k]:argv[k]));
      } 
    else 
      res = defaut;
    if (visu && usage) 
      std::fprintf(stderr,"    %-8s = %-12s : %s\n", utl::GetColoredString(std::string(name),COLOR_BOLD).c_str(),res?res:"NULL",
              utl::GetColoredString(std::string(usage),COLOR_PURPLE).c_str());
    }
  return res;
}

inline bool 
utlOption(const char *const name,const unsigned int argc,char **argv, const bool defaut,const char *const usage=NULL) 
{
  const char *s = utlOption(name,argc,argv,(const char*)NULL);
  const bool res = s?(strcmp(s,"false") && strcmp(s,"FALSE") && strcmp(s,"off") && strcmp(s,"OFF") && strcmp(s,"0")):defaut;
  utlOption(name,0,NULL,res?"true":"false",usage);
  return res;
}

inline int 
utlOption(const char *const name,const unsigned int argc,char **argv, const int defaut,const char *const usage=NULL) 
{
  const char *s = utlOption(name,argc,argv,(const char*)NULL);
  const int res = s?atoi(s):defaut;
  char tmp[256];
  std::sprintf(tmp,"%d",res);
  utlOption(name,0,NULL,tmp,usage);
  return res;
}

inline char 
utlOption(const char *const name,const unsigned int argc,char **argv, const char defaut,const char *const usage=NULL) 
{
  const char *s = utlOption(name,argc,argv,(const char*)NULL);
  const char res = s?s[0]:defaut;
  char tmp[8];
  tmp[0] = res;
  tmp[1] ='\0';
  utlOption(name,0,NULL,tmp,usage);
  return res;
}

inline double 
utlOption(const char *const name,const unsigned int argc,char **argv, const double defaut,const char *const usage=NULL) 
{
  const char *s = utlOption(name,argc,argv,(const char*)NULL);
  const double res = s?atof(s):defaut;
  char tmp[256];
  std::sprintf(tmp,"%g",res);
  utlOption(name,0,NULL,tmp,usage);
  return res;
}

}

#endif 


