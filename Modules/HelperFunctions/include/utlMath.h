/**
 *       @file  utlMath.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-29-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlMath_h
#define __utlMath_h

#include "utlCore.h"

namespace utl
{

/**
 * \brief  BesselJPrimeZeros in Mathematica, the k-th solution of J'_m(x)=0
 */
const static double BesselJPrimeZerosTable[] = 
  {
   // 10 in order 0
   0.0,3.8317059702075125,7.0155866698156055,10.173468135062722,13.323691936314223,16.47063005087763,19.615858510468243,22.760084380592772,25.903672087618382,29.046828534916855,           
   // 10 in order 1
   1.841183781340659,5.3314427735250325,8.536316366346284,11.706004902592063,14.863588633909034,18.015527862681804,21.16436985918879,24.311326857210776,27.457050571059245,30.601922972669094,
   // 10 in order 2
   3.0542369282271404,6.706133194158457,9.969467823087596,13.170370856016122,16.347522318321783,19.512912782488204,22.671581772477424,25.826037141785264,28.977672772993678,32.127327020443474,
   // 10 in order 3
   4.201188941210528,8.015236598375951,11.345924310742964,14.585848286167028,17.78874786606647,20.9724769365377,24.144897432909264,27.310057930204348,30.470268806290424,33.62694918279668,
   // 10 in order 4
   5.317553126083994,9.282396285241617,12.68190844263889,15.964107037731551,19.196028800048904,22.401032267689004,25.589759681386735,28.767836217666503,31.938539340972785,35.10391667734677,
   // 10 in order 5
   6.41561637570024,10.519860873772291,13.9871886301403,17.312842487884627,20.57551452138689,23.803581476593862,27.01030789777772,30.20284907898166,33.38544390101012,36.56077768688036,
   // 10 in order 6
   7.5012661446841475,11.734935953042752,15.268181461097873,18.637443009666203,21.931715017802237,25.183925599499627,28.409776362510083,31.617875716105036,34.81339298429743,37.9996408977153,
   // 10 in order 7
   8.57783648971416,12.93238623709143,16.529365884366943,19.941853366527344,23.26805292645757,26.545032061823576,29.790748583196613,33.015178641375144,36.22438054878716,39.42227457893926,
   // 10 in order 8
   9.647421651997242,14.11551890789478,17.774012366915255,21.229062622853125,24.58719748631768,27.889269427955092,31.155326556188324,34.396628554272176,37.620078044197086,40.83017868182204,
   // 10 in order 9
   10.711433970699948,15.286737667333524,19.004593537946054,22.501398726777285,25.891277276839137,29.21856349993608,32.50524735237553,35.7637929288088,39.00190281151422,42.22463843075328,
   // 10 in order 10
   11.770876674955586,16.447852748486536,20.223031412681703,23.760715860327448,27.182021527190532,30.534504754007074,33.84196577513571,37.118000423665606,40.37106890533389,43.60676490137951,
  };

/** pre-computed table for Gamma[1/2], Gamma[3/2], etc. Table[N[Gamma[n/2], 35], {n, 1, 60, 2}] in mathematica. */
const static double GammaHalfIntegerTable[30] =
  {
1.7724538509055160272981674833411452, 
0.88622692545275801364908374167057259, 
1.3293403881791370204736256125058589, 
3.3233509704478425511840640312646472, 
11.631728396567448929144224109426265, 
52.342777784553520181149008492418194, 
287.88527781504436099631954670830007, 
1871.2543057977883464760770536039504, 
14034.407293483412598570577902029628, 
119292.46199460900708784991216725184, 
1.1332783889487855673345741655888925*1e6, 
1.1899423083962248457013028738683371*1e7, 
1.3684336546556585725564983049485877*1e8, 
1.7105420683195732156956228811857346*1e9, 
2.3092317922314238411890908896007417*1e10, 
3.3483860987355645697241817899210754*1e11, 
5.1899984530401250830724817743776669*1e12, 
8.5634974475162063870695949277231504*1e13, 
1.4986120533153361177371791123515513*1e15, 
2.7724322986333718178137813578503700*1e16, 
5.4062429823350750447368736478082214*1e17, 
1.1082798113786903841710590978006854*1e19, 
2.3828015944641843259677770602714736*1e20, 
5.3613035875444147334274983856108156*1e21, 
1.2599063430729374623554621206185417*1e23, 
3.0867705405286967827708821955154271*1e24, 
7.8712648783481767960657495985643390*1e25, 
2.0858851927622668509574236436195498*1e27, 
5.7361842800962338401329150199537621*1e28, 
1.6348125198274266444378807806868222*1e30
  };

  
/** pre-computed table for factorial of integers. Table[n!, {n, 0, 30}] in mathematica  */
const static unsigned long FactorialTable[21]=
  {
  1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 
  39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 
  355687428096000, 6402373705728000, 121645100408832000,2432902008176640000
  // , 
  // 51090942171709440000, 1124000727777607680000, 
  // 25852016738884976640000, 620448401733239439360000, 
  // 15511210043330985984000000, 403291461126605635584000000, 
  // 10888869450418352160768000000, 304888344611713860501504000000, 
  // 8841761993739701954543616000000, 265252859812191058636308480000000
  };

/** calculate exp(-dist) using LUT */
inline double
ExpNegtiveLUT(const double dist, const double distMax=30.0, const int precision=1000 )
{
  if (dist >= distMax-1) 
    return 0.0;

  static int LUT_LENGTH = (int) rintf( distMax * (double) precision);
  static std::vector<double> EXP_LUT(LUT_LENGTH,-1.0);
  
  static bool is_firstTime = true;
  if (is_firstTime)
    {
    for (int i=0; i< LUT_LENGTH; i++)   
      EXP_LUT[i] = std::exp( - (double) i / (double)precision);
    is_firstTime = false;
    }

  double distPrecision = dist* (double)precision;
  int x = (int) std::floor(distPrecision);
  return EXP_LUT[x] + (EXP_LUT[x+1]-EXP_LUT[x])*(distPrecision-x);
}

/** efficient way to calculate std::pow(a,b) when b is integer or half integer */
inline double
PowHalfInteger(const double a, const double b)
{
  if (b==0)
    return 1;

  if (b<0)
    return 1.0/PowHalfInteger(a, -b);

  utlException(!utl::IsInt(2*b), "b need to be an integer or a half integer.");

  if (utl::IsInt(b))
    {
    long b_int = (long)b;
    double result = a;
    for ( long i = 1; i < b_int; i += 1 ) 
      result *= a;
    return result;
    }
  else
    {
    long b_int = long(b-0.5);
    double result = std::sqrt(a);
    for ( long i = 0; i < b_int; i += 1 ) 
      result *= a;
    return result;
    }
}

/** factorial of non-negative value n using a static table. */
inline unsigned long 
Factorial(const int n)
{
  utlException(n<0, "n should be non-negative");
  if(n <= 1) // n=0,1 return 1
    return 1;

  utlException(n>20, "n is too big for the table");
  
  return FactorialTable[n];
}

template <typename T>
T 
Factorial ( const T v1, const int times )
{
  T result = v1;
  for ( int i = 1; i < times; i += 1 ) 
    {
    result *= v1 - i;
    }
  return result;
}

/**
 * \brief  generalized binomial coefficients
 * \note   if times==0, return 1 whatever v1 is.
 */
template <typename T>
T 
Binomial ( const T v1, const int times )
{
  utlException(times<0, "negative value is invalid");
  if (times==0)
    return 1;
  T fac_1 = Factorial(v1,times);
  unsigned long fac_2 = Factorial(times);
  return fac_1 / T(fac_2);
}

double 
LegendrePolynomialAt0(const int order)
{
  if(order == 0)
    return 1.0;
  if(order % 2 != 0)
    return 0.0;
  else 
    {
    double odd = 1;
    double even = 2;

    for(int i = 3; i <= order - 1; i+=2) 
      odd = odd * i;
    for(int i = 4; i <= order; i+=2) 
      even = even * i;

    if((order / 2) % 2 == 0)
      return odd / even;
    else
      return -1 * odd / even;
    }
}

/** get the coefficient vector of nth order Lagurre polynomial L_n^{alpha}(x). 
 *      The default value of alpha is 0.5 for \c DiffusionMRI::EAP class. 
 * \note  the ith coefficient of n order Lagurre polynomial has a closed form 
 *     \f$ (-1)^i\binom{n+0.5}{n-i}\frac{1}{i!} \f$. 
 *     But here we use a recursive way that is more accurate when n is big.
 *   \f$L_n^a(x) =(2+(a-1-x)/n)*L_{n-1}^a(x) - (1+(a-1)/n)*L_{n-2}^a(x)\f$
 *   Also use pre-computed table Table[CoefficientList[LaguerreL[n, a, x], x], {n, 0, 10}] in mathematica.
 *     */ 
inline std::vector<double>
GetCoefLaguerre (const int n, const double a=0.5)
{
  utlException(n<0, "n should > 0");
  std::vector<double> coef(n+1);
  if (n==0)
    {
    coef[0] = 1;
    }
  else if (n==1)
    {
    coef[0]=1+a; coef[1]=-1.0;
    }
  else if (n==2)
    {
    double a2= a*a;
    coef[0]=1+(3*a)/2.0+a2/2.0; coef[1]=-2-a; coef[2]=1.0/2.0;
    }
  else if (n==3)
    {
    double a2= a*a;
    double a3= a2*a;
    coef[0]=1+(11*a)/6.0+a2+a3/6.0; coef[1]=-3-(5*a)/2.0-a2/2.0; 
    coef[2]=3.0/2.0+a/2.0; coef[3]=-1.0/6.0; 
    }
  else if (n==4)
    {
    double a2= a*a;
    double a3= a2*a;
    double a4= a3*a;
    coef[0]=1+ (25*a)/12.0+(35*a2)/24.0+(5*a3)/12.0+a4/24.0; 
    coef[1]=-4-(13*a)/3.0-(3*a2)/2.0-a3/6.0;
    coef[2]=3+(7*a)/4.0+a2/4.0; 
    coef[3]=-(2.0/3.0)-a/6.0; 
    coef[4]=1/24.0;
    }
  else if (n==5)
    {
    double a2= a*a;
    double a3= a2*a;
    double a4= a3*a;
    double a5= a4*a;
    coef[0]=1+(137*a)/60.0+(15*a2)/8.0+(17*a3)/24.0+a4/8.0+a5/120.0; 
    coef[1]=-5-(77*a)/12.0-(71*a2)/24.0-(7*a3)/12.0-a4/24.0; 
    coef[2]=5+(47*a)/12.0+a2+a3/12.0;
    coef[3]=-(5.0/3.0)-(3*a)/4.0-a2/12.0; 
    coef[4]=5.0/24.0+a/24.0;
    coef[5]=-(1/120.0);
    }
  else if (n==6)
    {
    double a2= a*a;
    double a3= a2*a;
    double a4= a3*a;
    double a5= a4*a;
    double a6= a5*a;
    coef[0]=1+(49*a)/20.0+(203*a2)/90.0+(49*a3)/48.0+(35*a4)/144.0+(7*a5)/240.0+a6/720.0;
    coef[1]=-6-(87*a)/10.0-(29*a2)/6.0-(31*a3)/24.0-a4/6.0-a5/120.0;
    coef[2]=15/2.0+(57*a)/8.0+(119*a2)/48.0+(3*a3)/8.0+a4/48.0;
    coef[3]=-(10/3.0)-(37*a)/18.0-(5*a2)/12.0-a3/36.0;
    coef[4]=5/8.0+(11*a)/48.0+a2/48.0;
    coef[5]=-(1/20.0)-a/120.0;
    coef[6]=1/720.0;
    }
  else if (n==7)
    {
    double a2= a*a;
    double a3= a2*a;
    double a4= a3*a;
    double a5= a4*a;
    double a6= a5*a;
    double a7= a6*a;
    coef[0]=1+(363*a)/140.0+(469*a2)/180.0+(967*a3)/720.0+(7*a4)/18.0+(23*a5)/360.0+a6/180.0+a7/5040.0;
    coef[1]=-7-(223*a)/20.0-(319*a2)/45.0-(37*a3)/16.0-(59*a4)/144.0-(3*a5)/80.0-a6/720.0;
    coef[2]=21/2.0+(459*a)/40.0+(235*a2)/48.0+(49*a3)/48.0+(5*a4)/48.0+a5/240.0;
    coef[3]=-(35/6.0)-(319*a)/72.0-(179*a2)/144.0-(11*a3)/72.0-a4/144.0;
    coef[4]=35/24.0+(107*a)/144.0+a2/8.0+a3/144.0;
    coef[5]=-(7/40.0)-(13*a)/240.0-a2/240.0;
    coef[6]=7/720.0+a/720.0;
    coef[7]-(1/5040);
    }
  else if (n==8)
    {
    double a2= a*a;
    double a3= a2*a;
    double a4= a3*a;
    double a5= a4*a;
    double a6= a5*a;
    double a7= a6*a;
    double a8= a7*a;
    coef[0]=1+(761*a)/280.0+(29531*a2)/10080.0+(267*a3)/160.0+(1069*a4)/1920.0+(9*a5)/80.0+(13*a6)/960.0+a7/1120.0+a8/40320.0;
    coef[1]=-8-(481*a)/35.0-(349*a2)/36.0-(329*a3)/90.0-(115*a4)/144.0-(73*a5)/720.0-a6/144.0-a7/5040.0;
    coef[2]=14+(341*a)/20.0+(6077*a2)/720.0+(209*a3)/96.0+(89*a4)/288.0+(11*a5)/480.0+a6/1440.0;
    coef[3]=-(28/3.0)-(743*a)/90.0-(23*a2)/8.0-(71*a3)/144.0-a4/24.0-a5/720.0;
    coef[4]=35/12.0+(533*a)/288.0+(251*a2)/576.0+(13*a3)/288.0+a4/576.0;
    coef[5]=-(7/15.0)-(73*a)/360.0-(7*a2)/240.0-a3/720.0;
    coef[6]=7/180.0+a/96.0+a2/1440.0;
    coef[7]=-(1/630.0)-a/5040.0;
    coef[8]=1/40320.0;
    }
  else
    {
    std::vector<double> coef_1 = GetCoefLaguerre(n-1,a);
    std::vector<double> coef_2 = GetCoefLaguerre(n-2,a);
    for ( int i = 1; i < n-1; i += 1 ) 
      coef[i] = coef_1[i]* (2.0+(a-1)/double(n)) + coef_1[i-1]*(-1.0/double(n)) - coef_2[i]*(1+(a-1)/n);
    coef[0] = coef_1[0]*(2.0+(a-1)/double(n)) - coef_2[0]*(1+(a-1)/n);
    coef[n-1] = coef_1[n-1]* (2.0+(a-1)/double(n)) + coef_1[n-2]*(-1.0/double(n));
    coef[n] = coef_1[n-1]*(-1.0/double(n));
    }
  return coef;
}    // -----  end of method getCoefLaguerre  -----

/** get the coefficient vector of the product of \f$L_{n1}^{a1}(x)\f$ and \f$L_{n2}^{a2}(x)\f$   */
inline std::vector<double>
GetCoefLaguerreProduct ( const int n1, const double a1, const int n2, const double a2)
{
  utlException(n1<0 || n2<0, "n1 and n2 should > 0");
  std::vector<double> coef(n1+n2+1,0);
  std::vector<double> coef1 = GetCoefLaguerre(n1,a1);
  std::vector<double> coef2 = GetCoefLaguerre(n2,a2);
  for ( int nn = 0; nn < n1+n2+1; nn += 1 ) 
    {
    for ( int j = 0; j<n1+1 && j<=nn; j += 1 ) 
      {
      int i = nn-j;
      if (i<n2+1)
        coef[nn] += coef1[j]*coef2[i];
      }
    }
  return coef;
}    // -----  end of method getCoefLaguerreProduct  -----

inline double 
GammaHalfInteger(const double x)
{
  utlException(x<=0, "x should be half of a positive integer");
  if (utl::IsInt(x))
    return Factorial(x-1);

  utlException(!utl::IsInt(2*x), "x should be half of a positive integer");
  
  const static int TABLE_LENGTH = 30; 
  
  if (x>TABLE_LENGTH)
    {
    double result=1, tmp=x-1;
    for ( long x_floor = std::floor(x); x_floor>=TABLE_LENGTH; x_floor-- ) 
      {
      result *= tmp;
      tmp -= 1;
      }
    return result*GammaHalfIntegerTable[TABLE_LENGTH-1];
    }
  else 
    {
    int x_floor = std::floor(x);
    return GammaHalfIntegerTable[x_floor];
    }
}

}


#endif 
