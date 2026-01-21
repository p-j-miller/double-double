/* Basic test program for double-doubles
   Written by Peter Miller 1-12-2025
   Tests check fma (which is assumed present and working by the code), then check double-double multiply and subtract, finally divide.
   These check the "lower level" double-double functions as well.
   Note these tests assume "#define USE_ACCURATE_ALG" has been used in double-doube.c, some (4) will fail if this is not present.
   double-double "power" functions have their own test program in "power10" and so are not checked here.
To compile under Windows using winlibs gcc 15.2.0 (the changes for Linux should be obvious):
  C:\winlibs\winlibs-x86_64-posix-seh-gcc-15.2.0-mingw-w64ucrt-13.0.0-r2\mingw64\bin\gcc -Wall -m64 -fexcess-precision=standard -Ofast  -std=gnu99 -I. main.c ../double-double/double-double.c ../u2_64-128bits-with-two-u64/u2_64.c ../my_printf/my_printf.c  ../fma/fmaq.c -lquadmath -static -o test.exe
  
 expected output:
sizeof float=4 double=8 long double=16
__SIZEOF_FLOAT__ is defined and set to 4
__SIZEOF_DOUBLE__ is defined and set to 8
__SIZEOF_LONG_DOUBLE__ is defined and is set to 16
__SIZEOF_FLOAT128__ is defined and set to 16
__SIZEOF_INT128__ is defined and set to 16
number of bits in mantissa of float is FLT_MANT_DIG=24  and max exponent is 2^FLT_MAX_EXP=128 max exponent (decimal)=FLT_MAX_10_EXP=38
number of bits in mantissa of double is DBL_MANT_DIG=53  and max exponent is 2^DBL_MAX_EXP=1024 max exponent (decimal)=DBL_MAX_10_EXP=308
number of bits in mantissa of long double is LDBL_MANT_DIG=64  and max exponent is 2^LDBL_MAX_EXP=16384 max exponent (decimal)=LDBL_MAX_10_EXP=4932
 max difference between 1 and next larger representable number for doubles is  DBL_EPSILON=((double)2.22044604925031308084726333618164062e-16L)
number of bits in mantissa of float128 is __FLT128_MANT_DIG__=113  and max exponent is 2^__FLT128_MAX_EXP__=16384 max exponent (decimal)=__FLT128_MAX_10_EXP__=4932
__GNUC__ defined, __GNUC__=15 __GNUC_MINOR__=2
__MINGW32__ defined and is set to 1
__USE_MINGW_ANSI_STDIO is defined as 0
Using UCRT for sprintf() &  snprintf()
__MINGW64__ defined and is set to 1
_UCRT defined and is set to
__MSVCRT__ defined and is set to 1
_WIN32 defined and is set to 1
_WIN64 defined and is set to 1
__x86_64 defined and is set to 1

Checking long double functions:
 long double functions: 0 errors

Checking fma functions:
double (DBL_EPSILON = 2.22044604925031308084726333618164062e-16 so e^2=4.9303806576313237838233035330172e-32 :
1+epsilon=1.00000000000000022e+00 1+2*epsilon=1.00000000000000044e+00
   direct=0.00000000000000000e+00
      fma=4.93038065763132378e-32
     f128=4.930380657631323783823303533017414e-32
 all f128=4.930380657631323783823303533017414e-32
  - double check passed

long double (LDBL_EPSILON = 1.08420217248550443400745280086994171e-19 so e^2=1.1754943508222875079687365372222e-38 :
1+epsilon=1.000000000000000000108e+00 1+2*epsilon=1.000000000000000000217e+00
   direct=0.000000000000000000000e+00
     fmal=1.175494350822287507969e-38
     f128=0.000000000000000000000000000000000e+00
 all f128=1.175494350822287507968736537222246e-38
  - long double check passed

float128 (FLT128_EPSILON = 1.92592994438723585305597794258492732e-34 so e^2=3.7092061506874213857317352615475e-68 :
1+epsilon=1.000000000000000000000000000000000193e+00 1+2*epsilon=1.000000000000000000000000000000000385e+00
   direct=0.000000000000000000000000000000000000e+00
     fmaq=3.709206150687421385731735261547639513e-68
  - float128 check passed

double-double (DBL_EPSILON = 2.22044604925031308084726333618164062e-16 so e^2=4.9303806576313237838233035330172e-32 :
 (1+e)*(1+e) - (1+2e) gives:
   double-double=4.93038065763132378e-32 0.00000000000000000e+00
            f128=4.930380657631323783823303533017414e-32
 (1+e)*(1-e) gives (should be 1-e^2):
   double-double=1.00000000000000000e+00 -4.93038065763132378e-32
            f128=9.999999999999999999999999999999507e-01
          f128-1=-4.930380657631323783823303533017414e-32
  - double-double check passed

long double-double (LDBL_EPSILON = 1.08420217248550443400745280086994171e-19 so e^2=1.1754943508222875079687365372222e-38 :
 (1+e)*(1+e) - (1+2e) gives:
   double-double=1.175494350822287507969e-38 0.000000000000000000000e+00
            f128=1.175494350822287507968736537222246e-38
 (1+e)*(1-e) gives (should be 1-e^2):
   double-double=1.000000000000000000000e+00 -1.175494350822287507969e-38
          f128-1=-1.175494350822287507968736537222246e-38
  - long double-double check passed

float128-double-double: (FLT128_EPSILON = 1.92592994438723585305597794258492732e-34 so e^2=3.7092061506874213857317352615475e-68 :
1+epsilon=1.000000000000000000000000000000000193e+00 1+2*epsilon=1.000000000000000000000000000000000385e+00
 (1+e)*(1+e) - (1+2e) gives:
                 direct= 0.000000000000000000000000000000000000e+00
                   fmaq= 3.709206150687421385731735261547639513e-68
  as f128 double-double= 3.709206150687421385731735261547639513e-68 , 0.000000000000000000000000000000000000e+00
 (1+e)*(1-e) gives (should be 1-e^2):
  as f128 double-double= 1.000000000000000000000000000000000000e+00 , -3.709206150687421385731735261547639513e-68
  - double-double-float128 check passed


Checking double double multiply and divide for accuracy
 355/113 (approximation to pi):
 as dd =3.14159292035398252e+00,-2.20079608421269966e-16
 dd's combined to f128's=3.141592920353982300884955752212392e+00
               as f128's=3.141592920353982300884955752212389e+00
              Difference=3.081487911019577364889564708135884e-33
      Test passed
 (355+e)/(113) (e=double epsilon):
 as dd =3.14159292035398252e+00,-2.18114611917508624e-16
 dd's combined to f128's=3.141592920353982302849952255973734e+00
               as f128's=3.141592920353982302849952255973728e+00
              Difference=6.162975822039154729779129416271767e-33
      Test passed
 (355+e)/(113+e) (e=double epsilon):
 as dd =3.14159292035398208e+00,2.19801378827817027e-16
 dd's combined to f128's=3.141592920353982296676733151236769e+00
               as f128's=3.141592920353982296676733151236778e+00
              Difference=-9.629649721936179265279889712924637e-33
      Test passed
 (355+e/1e+14)/(113+e/1e+14) (e=double epsilon):
 as dd =3.14159292035398252e+00,-2.20079608421270015e-16
 dd's combined to f128's=3.141592920353982300884955752212343e+00
               as f128's=3.141592920353982300884955752212347e+00
              Difference=-4.237045877651918876723151473686840e-33
      Test passed

 Long-double-double:
 (355+e/1e+15)/(113+e/1e+15) (e=double epsilon):
 as ldd =3.141592920353982300872e+00,1.343259328742682016264e-20
 ldd's combined to f128's=3.141592920353982300884955752212385e+00
                as f128's=3.141592920353982300884955752212386e+00
               Difference=-3.851859888774471706111955885169855e-34
      Test passed
  above result calculated as double-double f128
                         =3.141592920353982300884955752212385e+00,5.609308119860504523507473129697090e-35
      Test passed


All tests passed

  
*/
/*----------------------------------------------------------------------------
 * MIT License:
 *
 * Copyright (c) 2025,2026 Peter Miller
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHOR OR COPYRIGHT HOLDER BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *--------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h> /* see https://gcc.gnu.org/onlinedocs/libquadmath/quadmath_005fsnprintf.html#quadmath_005fsnprintf - also needs quadmath library linking in */
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <inttypes.h>
#include <ctype.h>
#include <stdarg.h>
#include "../my_printf/my_printf.h"
#include "../double-double/double-double.h"


#define _mkstr(s) # s
#define mkstr(s) _mkstr(s)      /* creates "s" */

/* the line below defines GCC_OPTIMIZE_AWARE when we can use # pragma GCC optimize ("-O2") */
#define GCC_OPTIMIZE_AWARE (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) || defined(__clang__)
/* code below cannot be compiled with -Ofast as this makes the compiler break some C rules that we need, so make sure of this here */
#if GCC_OPTIMIZE_AWARE
 #pragma GCC push_options
 #pragma GCC optimize ("-O3") /* cannot use Ofast, normally -O3 is OK. Note macro expansion does not work here ! */
 #if defined(_WIN32) && !defined(_WIN64)
  #pragma GCC target("sse2")
 #endif
#endif

unsigned int errs=0;

#if 1 /* check a sample of long double maths functions - tests atanl(), sinl(),cosl(),sinhl,coshl(),sqrtl(),powl(),fabsl() as well as +,-,*,/  */

/* this is nothing to do with double-doubles directly, here just in case these functions need to be checked  */
static void check_mldbl(long double ld)
{long double lr,lr1,lpi=atanl(1)+atanl(2)+atanl(3);
 double d=ld,r,r1,pi=atan(1)+atan(2)+atan(3);
 /* from   https://en.wikipedia.org/wiki/List_of_trigonometric_identities pi=arctan(1)+arctan(2)+arctan(3) - 1,2,3 can all b expressed exactly as fp numbers so this should be very accurate */
 if(lpi<3.1415926 || lpi> 3.1415927) // PI=3.14159265....
    {++errs;
     printf("  arctan(1)+arctan(2)+arctan(3)=%Lg should = PI [ as doubles = %g]\n",lpi,pi);
    }
 if((lr=fabsl(lpi-pi))>1e-15) // >1e-15 double should give at least 16 sig figs re 3.xxxxx
    {++errs;
     printf("  arctan(1)+arctan(2)+arctan(3)=%Lg should = PI as doubles = %g, difference=%Lg\n",lpi,pi,lr);
    }
 /* sin(x)^2+cos(x)^2=1 */
 lr=sinl(ld)*sinl(ld)+cosl(ld)*cosl(ld);
 r=sin(d)*sin(d)+cos(d)*cos(d);
 if(lr>1.00000001 || lr<0.99999999)
    {++errs;
     printf("  sin(%Lg)^2+cos(%Lg)^2 = %Lg should = 1 [ as doubles = %g]\n",ld,ld,lr,r);
    }
 lr1=fabsl(lr-1.0L);
 r1=fabs(r-1.0);
 if(r1>DBL_EPSILON && lr1>r1 ) // long double error should be no worse than double error (assuming double error is not 0)
    {++errs;
     printf("  sin(%Lg)^2+cos(%Lg)^2 = %Lg (error %Lg) should = 1 [ as doubles = %g (error=%g)]\n",ld,ld,lr,lr1,r,r1);
    }

 /* cosh(x)^2-sinh(x)^2=1 note sinh(x)=(e^x-e^-x)/2 and cosh(x)=(e^x+e^-x)/2 so we need to be careful of overflow  and loss of resolution due to the subtraction */
 if(fabsl(ld)>3.0L) {ld=3.0L/ld; d=ld;} // try and avoid overflow, |x|<=3 keeps sinh(x) and cosh(x)< ~ +/-10
 lr=coshl(ld)*coshl(ld)-sinhl(ld)*sinhl(ld);
 r=cosh(d)*cosh(d)-sinh(d)*sinh(d);
 if(lr>1.00000001 || lr<0.99999999)
    {++errs;
     printf("  cosh(%Lg)^2-sinh(%Lg)^2 = %.10Lg should = 1 [ as doubles = %.10g]\n",ld,ld,lr,r);
    }
 lr1=fabsl(lr-1.0L);
 r1=fabs(r-1.0);
 if(r1>DBL_EPSILON && lr1>r1) // long double error should be no worse than double error, unless double error is 0
    {++errs;
     printf("  cosh(%Lg)^2-sinh(%Lg)^2 = %Lg (error %Lg) should = 1 [ as doubles = %g (error=%g)]\n",ld,ld,lr,lr1,r,r1);
    }

 // sin(arctan(x))=x/sqrt(1+x^2)  - https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Inverse_trigonometric_functions     -
 lr=sinl(atanl(ld));   // in long double
 lr1=ld/sqrtl(1.0L+powl(ld,2.0L));
 r=sin(atan(d));    // in double
 r1=d/sqrt(1.0+pow(d,2.0));
 lr=fabsl(lr-lr1); // abs error for long double
 r=fabs(r-r1); // abs error for double
 if(r>0 && lr>r) // long double error should be no worse than double error, unless double error is 0 (cannot use DBL_EPSILON here as correct result is not "1" )
    {++errs;
     printf("  sin(arctan(x))=x/sqrt(1+x^2) for x=%Lg:error for long double=%Lg, error for doubles = %g\n",ld,lr,r);
    }

}
#endif

static void check_fma(void) // do some checks on fma(), fmal() amd fmaq()
{
 // start double :
 // DBL_DIG=15,DBL_DECIMAL_DIG__ = 17 , __DBL_EPSILON__ = ((double)2.22044604925031308084726333618164062e-16L)
 // use (1+e) * (1+e) = 1+2e+e^2
 // fma(a,b,c) = a*b-c with just a single rounding on the final result
 // we check fma(1+e,1+e,-(1+2e)) which should equal e^2
 char buf128[128],buf128_e[128];
 #define DBL_DECIMAL_DIG __DBL_DECIMAL_DIG__
 {
  double one_e=1.0+DBL_EPSILON;
  double one_2e=one_e+DBL_EPSILON;
  double r_direct=(one_e*one_e)-one_2e;
  double r_fma_d=fma(one_e,one_e,-one_2e);
  __float128 r_128=(__float128)(one_e)*(__float128)(one_e)-(__float128)(one_2e);// "fma" calculated exactly using f128's
  __float128 r_exact_128=(1.0f128+2.22044604925031308084726333618164062e-16f128)*(1.0f128+2.22044604925031308084726333618164062e-16f128)-(1.0f128+2.0f128*2.22044604925031308084726333618164062e-16f128);// all values f128 - so "exact"
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  quadmath_snprintf (buf128_e, sizeof buf128,"%.33Qe",r_exact_128); // 33 is full precision for full 128bit ieee value
  printf("double (DBL_EPSILON = 2.22044604925031308084726333618164062e-16 so e^2=4.9303806576313237838233035330172e-32 :\n");//  e^2 value calculated by "Windows calculator" using its maximum resolution
  printf("1+epsilon=%.*e 1+2*epsilon=%.*e\n",DBL_DECIMAL_DIG,one_e,DBL_DECIMAL_DIG,one_2e);
  printf("   direct=%.*e\n",DBL_DECIMAL_DIG,r_direct);
  printf("      fma=%.*e\n",DBL_DECIMAL_DIG,r_fma_d);
  printf("     f128=%s\n",buf128);
  printf(" all f128=%s\n",buf128_e); 
  if(r_direct==0 && r_fma_d==(double)r_128 && (double)r_128==(double)r_exact_128 && r_fma_d==4.9303806576313237838233035330172e-32) 
  	printf("  - double check passed\n");
  else
  	{printf("  - double check failed\n");  
  	 errs++;
  	}
  printf("\n");
 }
 // now long double:
 // __LDBL_DECIMAL_DIG__ 21 __LDBL_DIG__ 18
 // __LDBL_EPSILON__ 1.08420217248550443400745280086994171e-19L
 // __LDBL_MANT_DIG__ 64 __FLT128_MANT_DIG__ 113 => 2*64=128 so a float128 by itself does not have enough mantissa bits to exactly calculate the result (__DBL_MANT_DIG__ 53 so float128 is ok for that)
 // remember: you need to use my_printf() to correctly print longs...
 #define LDBL_DECIMAL_DIG __LDBL_DECIMAL_DIG__
 {
  long double one_e=1.0+LDBL_EPSILON;
  long double one_2e=one_e+LDBL_EPSILON;
  long double r_direct=(one_e*one_e)-one_2e;
  long double r_fma_d=fmal(one_e,one_e,-one_2e);
  __float128 r_128=(__float128)(one_e)*(__float128)(one_e)-(__float128)(one_2e);// "fma" calculated exactly using f128's
  // __float128 r_exact_128=(1.0f128+1.08420217248550443400745280086994171e-19f128)*(1.0f128+1.08420217248550443400745280086994171e-19f128)-(1.0f128+2.0f128*1.08420217248550443400745280086994171e-19f128);// f128 does not have enough mantissa bits to make this exact
  __float128 r_exact_128=fmaq((1.0f128+1.08420217248550443400745280086994171e-19f128),(1.0f128+1.08420217248550443400745280086994171e-19f128),-(1.0f128+2.0f128*1.08420217248550443400745280086994171e-19f128));
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  quadmath_snprintf (buf128_e, sizeof buf128,"%.33Qe",r_exact_128); // 33 is full precision for full 128bit ieee value
  my_printf("long double (LDBL_EPSILON = 1.08420217248550443400745280086994171e-19 so e^2=1.1754943508222875079687365372222e-38 :\n");//  e^2 value calculated by "Windows calculator" using its maximum resolution
  my_printf("1+epsilon=%.*Le 1+2*epsilon=%.*Le\n",LDBL_DECIMAL_DIG,one_e,LDBL_DECIMAL_DIG,one_2e);
  my_printf("   direct=%.*Le\n",LDBL_DECIMAL_DIG,r_direct);
  my_printf("     fmal=%.*Le\n",LDBL_DECIMAL_DIG,r_fma_d);
  my_printf("     f128=%s\n",buf128);// expect this to give zero
  my_printf(" all f128=%s\n",buf128_e); // using fmaq => note the library one works, so don't need to include my own version..
  if(r_direct==0 && r_fma_d==(long double)r_exact_128 && r_128==0 && r_fma_d==1.1754943508222875079687365372222e-38L) 
  	printf("  - long double check passed\n");
  else
  	{
  	 printf("  - long double check failed\n");  
  	 errs++;
  	}  	 
  printf("\n");  
  
 }
 
 // now __float128: - this is simpler as there is nothing more accurate to compare it with (other than fmaq())
 // __FLT128_DECIMAL_DIG__ 36 __FLT128_DIG__ 33
 // __FLT128_EPSILON__ 1.92592994438723585305597794258492732e-34F128

 #define FLT128_DECIMAL_DIG __FLT128_DECIMAL_DIG__
 {
  __float128 one_e=1.0+FLT128_EPSILON;
  __float128 one_2e=one_e+FLT128_EPSILON;
  __float128 r_direct=(one_e*one_e)-one_2e;
  __float128 r_fma_d=fmaq(one_e,one_e,-one_2e);

  my_printf("float128 (FLT128_EPSILON = 1.92592994438723585305597794258492732e-34 so e^2=3.7092061506874213857317352615475e-68 :\n"); //  e^2 value calculated by "Windows calculator" using its maximum resolution
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,one_e); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,one_2e);  
  my_printf("1+epsilon=%s 1+2*epsilon=%s\n",buf128,buf128_e);
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r_direct); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r_fma_d); 
  my_printf("   direct=%s\n",buf128);
  my_printf("     fmaq=%s\n",buf128_e);
  if(r_direct==0 && (long double)r_fma_d==3.7092061506874213857317352615475e-68L) 
  	printf("  - float128 check passed\n");
  else
  	{
  	 printf("  - float128 check failed\n");  
  	 errs++;
  	}  	 
  printf("\n");  
 } 
 
 // now look at double-double functions
 // start double :
 // DBL_DIG=15,DBL_DECIMAL_DIG__ = 17 , __DBL_EPSILON__ = ((double)2.22044604925031308084726333618164062e-16L)
 // use (1+e) * (1+e) = 1+2e+e^2
 // fma(a,b,c) = a*b-c with just a single rounding on the final result
 // we check fma(1+e,1+e,-(1+2e)) which should equal e^2
 {
  double one_e_h=1.0,one_e_l=DBL_EPSILON;
  double one_2e_h=1.0,one_2e_l=one_e_l+one_e_l;

  __float128 r_128=((__float128)(one_e_h)+(__float128)(one_e_l))*((__float128)(one_e_h)+(__float128)(one_e_l))-((__float128)(one_2e_h)+(__float128)(one_2e_l));// "fma" calculated exactly using f128's
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
//
  printf("double-double (DBL_EPSILON = 2.22044604925031308084726333618164062e-16 so e^2=4.9303806576313237838233035330172e-32 :\n");//  e^2 value calculated by "Windows calculator" using its maximum resolution
  //printf("1+epsilon=%.*e 1+2*epsilon=%.*e\n",DBL_DECIMAL_DIG,one_e,DBL_DECIMAL_DIG,one_2e);
  printf(" (1+e)*(1+e) - (1+2e) gives:\n");
  double r1h,r1l;
  mult_dd_dd(&r1h,&r1l,one_e_h,one_e_l,one_e_h,one_e_l); //(1+e)*(1+e)
  // printf("  mult_dd_dd =%.*e %.*e\n",DBL_DECIMAL_DIG,r1h,DBL_DECIMAL_DIG,r1l);
  sub_dd_dd(&r1h,&r1l,r1h,r1l,one_2e_h,one_2e_l);
  printf("   double-double=%.*e %.*e\n",DBL_DECIMAL_DIG,r1h,DBL_DECIMAL_DIG,r1l);
  printf("            f128=%s\n",buf128);
  // now do (1+e)*(1-e) = 1+e^2
  printf(" (1+e)*(1-e) gives (should be 1-e^2):\n");
  double r2h,r2l;
  mult_dd_dd(&r2h,&r2l,1.0,DBL_EPSILON,1.0,-DBL_EPSILON); //(1+e)*(1-e)
  printf("   double-double=%.*e %.*e\n",DBL_DECIMAL_DIG,r2h,DBL_DECIMAL_DIG,r2l);
  //r_128=((__float128)(one_e_h)+(__float128)(one_e_l))*((__float128)(one_e_h)-(__float128)(one_e_l));// calculated exactly using f128's
  r_128=(1.0f128+2.22044604925031308084726333618164062e-16f128)*(1.0f128-2.22044604925031308084726333618164062e-16f128);// calculated exactly using f128's
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  printf("            f128=%s\n",buf128);  
  // calc r_128-1
  r_128-=1.0f128;// calculated exactly using f128's
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  printf("          f128-1=%s\n",buf128);    

  if(r1h==4.9303806576313237838233035330172e-32 && r1l==0.0 && (double)r_128==(double)((__float128)r2h+(__float128)r2l - 1.0f128) && (double)r_128==-4.9303806576313237838233035330172e-32) 
  	printf("  - double-double check passed\n");
  else
  	{
  	 printf("  - double-double check failed\n");   	
  	 errs++;
  	}  	 
  printf("\n");
 } 
 // now for long double-long double 
 {
  long double one_e_h=1.0,one_e_l=LDBL_EPSILON;
  long double one_2e_h=1.0,one_2e_l=one_e_l+one_e_l;
  __float128 r_128=fmaq((1.0f128+1.08420217248550443400745280086994171e-19f128),(1.0f128+1.08420217248550443400745280086994171e-19f128),-(1.0f128+2.0f128*1.08420217248550443400745280086994171e-19f128));
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  my_printf("long double-double (LDBL_EPSILON = 1.08420217248550443400745280086994171e-19 so e^2=1.1754943508222875079687365372222e-38 :\n");//  e^2 value calculated by "Windows calculator" using its maximum resolution  

  printf(" (1+e)*(1+e) - (1+2e) gives:\n");
  long double r1h,r1l;
  mult_ldd_ldd(&r1h,&r1l,one_e_h,one_e_l,one_e_h,one_e_l); //(1+e)*(1+e)
  sub_ldd_ldd(&r1h,&r1l,r1h,r1l,one_2e_h,one_2e_l);
  my_printf("   double-double=%.*Le %.*Le\n",LDBL_DECIMAL_DIG,r1h,LDBL_DECIMAL_DIG,r1l);
  printf("            f128=%s\n",buf128);
  // now do (1+e)*(1-e) = 1+e^2
  printf(" (1+e)*(1-e) gives (should be 1-e^2):\n");
  long double r2h,r2l;
  mult_ldd_ldd(&r2h,&r2l,1.0L,LDBL_EPSILON,1.0L,-LDBL_EPSILON); //(1+e)*(1-e)
  my_printf("   double-double=%.*Le %.*Le\n",LDBL_DECIMAL_DIG,r2h,LDBL_DECIMAL_DIG,r2l);

  r_128=fmaq((1.0f128+1.08420217248550443400745280086994171e-19f128),(1.0f128-1.08420217248550443400745280086994171e-19f128),-1.0f128);// calculated exactly using f128's - need fmaq so have to calculate "-1"
  quadmath_snprintf (buf128, sizeof buf128,"%.33Qe",r_128); // 33 is full precision for full 128bit ieee value
  printf("          f128-1=%s\n",buf128);    

  if(r1h==1.1754943508222875079687365372222e-38L  && r1l==0.0 && (double)r_128==(double)((__float128)r2l) && r2h==1.0L && (double)r_128==-1.1754943508222875079687365372222e-38L) 
  	printf("  - long double-double check passed\n");
  else
	{  
  	 printf("  - long double-double check failed\n");   	
  	 errs++;
  	}  	 
  printf("\n");
 } 

 // now for f128-f128 
 // __FLT128_DECIMAL_DIG__ 36 __FLT128_DIG__ 33
 // __FLT128_EPSILON__ 1.92592994438723585305597794258492732e-34F128
 {
  __float128 one_e=1.0+FLT128_EPSILON;
  __float128 one_2e=one_e+FLT128_EPSILON;
  __float128 one_e_h=1.0f128,one_e_l=FLT128_EPSILON;
  __float128 one_2e_h=1.0f128,one_2e_l=FLT128_EPSILON+FLT128_EPSILON;
  __float128 r_direct=(one_e*one_e)-one_2e;
  __float128 r_fma_d=fmaq(one_e,one_e,-one_2e);

  my_printf("float128-double-double: (FLT128_EPSILON = 1.92592994438723585305597794258492732e-34 so e^2=3.7092061506874213857317352615475e-68 :\n"); //  e^2 value calculated by "Windows calculator" using its maximum resolution
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,one_e); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,one_2e);  
  my_printf("1+epsilon=%s 1+2*epsilon=%s\n",buf128,buf128_e);
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r_direct); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r_fma_d); 
  my_printf(" (1+e)*(1+e) - (1+2e) gives:\n"); 
  my_printf("                 direct= %s\n",buf128);
  my_printf("                   fmaq= %s\n",buf128_e);
  
  __float128 r1h,r1l;
  f128_mult_dd_dd(&r1h,&r1l,one_e_h,one_e_l,one_e_h,one_e_l); //(1+e)*(1+e)
  f128_sub_dd_dd(&r1h,&r1l,r1h,r1l,one_2e_h,one_2e_l);  
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r1h); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r1l);   
  printf("  as f128 double-double= %s , %s\n",buf128,buf128_e);
  
  // now do (1+e)*(1-e) = 1+e^2
  printf(" (1+e)*(1-e) gives (should be 1-e^2):\n");
  __float128 r2h,r2l;
  f128_mult_dd_dd(&r2h,&r2l,one_e_h,one_e_l,one_e_h,-one_e_l); //(1+e)*(1-e)
  // f128_sub_dd_dd(&r1h,&r1l,r1h,r1l,one_2e_h,one_2e_l);  
  quadmath_snprintf (buf128, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r2h); 
  quadmath_snprintf (buf128_e, sizeof buf128,"%.*Qe",FLT128_DECIMAL_DIG,r2l);   
  printf("  as f128 double-double= %s , %s\n",buf128,buf128_e);
  
  if((long double)r1h==3.7092061506874213857317352615475e-68L && r1l==0.0 &&  r2h==1.0L && (long double)r2l==-3.7092061506874213857317352615475e-68L) 
  	printf("  - double-double-float128 check passed\n");
  else
  	{
  	 printf("  - double-double-float128 check failed\n");  
  	 errs++;
  	}  	 
  printf("\n");  
 } 

}
 
 

 /* __DBL_DENORM_MIN__ ((double)4.94065645841246544176568792868221372e-324L) while__DBL_MIN__ ((double)2.22507385850720138309023271733240406e-308L) so exponents below -308 start to denormalise the mantissa 
   __DBL_MAX__ ((double)1.79769313486231570814527423731704357e+308L) so 308 is max positive exponent before overflow 
   __DBL_EPSILON__ ((double)2.22044604925031308084726333618164062e-16L)
 */
#if 1
void test_div(void)
 {// start with a simple division test use 355/113 which is a good approximation for pi, then add small constants to top & bottom to validate the resolution.
  double dh,dl;
  char buf[128],buf1[128]; 
  __float128 a128,e128;// a128 is approximation, e128 is exact result
  printf(" 355/113 (approximation to pi):\n");
  div_dd_dd(&dh,&dl,355.0,0.0,113.0,0.0);
  my_printf(" as dd =%.*e,%.*e\n",DBL_DECIMAL_DIG,dh,DBL_DECIMAL_DIG,dl);
  a128=(__float128)dh+(__float128)dl; 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,a128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf(" dd's combined to f128's=%s\n",buf);
  e128=355.0f128/113.0f128;
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("               as f128's=%s\n",buf);
  e128=a128-e128; // difference
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("              Difference=%s\n",buf);
  if(fabsq(e128)<3.0815e-33f128)
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	}  	
  	
  printf(" (355+e)/(113) (e=double epsilon):\n");
  div_dd_dd(&dh,&dl,355.0,DBL_EPSILON,113.0,0.0);
  my_printf(" as dd =%.*e,%.*e\n",DBL_DECIMAL_DIG,dh,DBL_DECIMAL_DIG,dl);
  a128=(__float128)dh+(__float128)dl; 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,a128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf(" dd's combined to f128's=%s\n",buf);
  e128=(355.0f128 + (__float128)DBL_EPSILON)/(113.0f128);
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("               as f128's=%s\n",buf);    
  e128=a128-e128; // difference
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("              Difference=%s\n",buf); 
  if(fabsq(e128)<6.163e-33f128)
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	}  	
  	  
  printf(" (355+e)/(113+e) (e=double epsilon):\n");
  div_dd_dd(&dh,&dl,355.0,DBL_EPSILON,113.0,DBL_EPSILON);
  my_printf(" as dd =%.*e,%.*e\n",DBL_DECIMAL_DIG,dh,DBL_DECIMAL_DIG,dl);
  a128=(__float128)dh+(__float128)dl; 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,a128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf(" dd's combined to f128's=%s\n",buf);
  e128=(355.0f128 + (__float128)DBL_EPSILON)/(113.0f128+(__float128)DBL_EPSILON);
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("               as f128's=%s\n",buf);  
  e128=a128-e128; // difference
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("              Difference=%s\n",buf); 
  if(fabsq(e128)<9.62965e-33f128)
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	}  	
  	  
  
  const double c=1e14; // 1e14 is the largest power of 10 that produces a change in the results for double-doubles 
  printf(" (355+e/%g)/(113+e/%g) (e=double epsilon):\n",c,c);
  div_dd_dd(&dh,&dl,355.0,DBL_EPSILON/c,113.0,DBL_EPSILON/c);
  my_printf(" as dd =%.*e,%.*e\n",DBL_DECIMAL_DIG,dh,DBL_DECIMAL_DIG,dl);
  a128=(__float128)dh+(__float128)dl; 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,a128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf(" dd's combined to f128's=%s\n",buf);
  e128=(355.0f128 + (__float128)DBL_EPSILON/(__float128)c)/(113.0f128+(__float128)DBL_EPSILON/(__float128)c);
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("               as f128's=%s\n",buf);    
  e128=a128-e128; // difference
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("              Difference=%s\n",buf);  
  if(fabsq(e128)<4.2371e-33f128)
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	}  	
  	  
  // repeat with long doubles
  printf("\n Long-double-double:\n");	  
  long double Lh,Ll;
  const long double cL=1e15L; // 1e15 is the largest power of 10 that produces a change in the results for long double-doubles 
  my_printf(" (355+e/%Lg)/(113+e/%Lg) (e=double epsilon):\n",cL,cL);
  div_ldd_ldd(&Lh,&Ll,355.0,DBL_EPSILON/cL,113.0,DBL_EPSILON/cL);
  my_printf(" as ldd =%.*Le,%.*Le\n",LDBL_DECIMAL_DIG,Lh,LDBL_DECIMAL_DIG,Ll);
  a128=(__float128)Lh+(__float128)Ll; 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,a128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf(" ldd's combined to f128's=%s\n",buf);
  e128=(355.0f128 + (__float128)DBL_EPSILON/(__float128)cL)/(113.0f128+(__float128)DBL_EPSILON/(__float128)cL);
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("                as f128's=%s\n",buf);    
  e128=a128-e128; // difference
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,e128);	 // use __FLT128_DIG__ here (33) as all these should be exact
  my_printf("               Difference=%s\n",buf);  
  if(fabsq(e128)<3.8519e-34f128)
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	} 
    	
  // basic check with double-double f128's
  __float128 fh,fl;
  f128_div_dd_dd(&fh,&fl,355.0f128,(__float128)DBL_EPSILON/(__float128)cL,113.0f128,(__float128)DBL_EPSILON/(__float128)cL); 
  quadmath_snprintf(buf, sizeof buf,"%.*Qe",__FLT128_DIG__,fh);
  quadmath_snprintf(buf1, sizeof buf1,"%.*Qe",__FLT128_DIG__,fl);  
  my_printf("  above result calculated as double-double f128\n");
  my_printf("                         =%s,%s\n",buf,buf1);
  if(fabsq(a128-fh)<=fabsq(e128))
  	printf("      Test passed\n");
  else
  	{printf("      Test failed\n");
  	 errs++;
  	}   
  printf("\n\n");
 }
#endif 
 


int main(int argc, char *argv[])
{
 printf("sizeof float=%zd double=%zd long double=%zd\n",sizeof(float),sizeof(double),sizeof(long double));

#if defined(__SIZEOF_FLOAT__)
 printf("__SIZEOF_FLOAT__ is defined and set to %d\n",__SIZEOF_FLOAT__);
#endif
#if defined(__SIZEOF_DOUBLE__)
 printf("__SIZEOF_DOUBLE__ is defined and set to %d\n",__SIZEOF_DOUBLE__);
#endif
#if defined(__SIZEOF_LONG_DOUBLE__)
 printf("__SIZEOF_LONG_DOUBLE__ is defined and is set to %d\n",__SIZEOF_LONG_DOUBLE__);
#endif
#if defined(__SIZEOF_FLOAT128__)
 printf("__SIZEOF_FLOAT128__ is defined and set to %d\n",__SIZEOF_FLOAT128__);
#endif
#if defined(__SIZEOF_INT128__)
 printf("__SIZEOF_INT128__ is defined and set to %d\n",__SIZEOF_INT128__);
#endif
 printf("number of bits in mantissa of float is FLT_MANT_DIG=%d  and max exponent is 2^FLT_MAX_EXP=%d max exponent (decimal)=FLT_MAX_10_EXP=%d\n",FLT_MANT_DIG, FLT_MAX_EXP,FLT_MAX_10_EXP );
 printf("number of bits in mantissa of double is DBL_MANT_DIG=%d  and max exponent is 2^DBL_MAX_EXP=%d max exponent (decimal)=DBL_MAX_10_EXP=%d\n",DBL_MANT_DIG, DBL_MAX_EXP,DBL_MAX_10_EXP );
 printf("number of bits in mantissa of long double is LDBL_MANT_DIG=%d  and max exponent is 2^LDBL_MAX_EXP=%d max exponent (decimal)=LDBL_MAX_10_EXP=%d\n",LDBL_MANT_DIG, LDBL_MAX_EXP,LDBL_MAX_10_EXP );
 printf(" max difference between 1 and next larger representable number for doubles is  DBL_EPSILON=%s\n",mkstr(DBL_EPSILON));
#ifdef __FLT128_MANT_DIG__
 printf("number of bits in mantissa of float128 is __FLT128_MANT_DIG__=%d  and max exponent is 2^__FLT128_MAX_EXP__=%d max exponent (decimal)=__FLT128_MAX_10_EXP__=%d\n",__FLT128_MANT_DIG__, __FLT128_MAX_EXP__,__FLT128_MAX_10_EXP__ );
#endif
#ifdef  __BORLANDC__
  printf("__BORLANDC__ is defined as\"%s\"\n",mkstr(__BORLANDC__));
#endif
#ifdef __clang__
 printf("__clang__ defined and is set to %s\n",mkstr(__clang__));
#endif
#ifdef __GNUC__
 printf("__GNUC__ defined, __GNUC__=%u __GNUC_MINOR__=%u\n",__GNUC__ ,__GNUC_MINOR__);
#endif
#ifdef __MINGW32__ 
 printf("__MINGW32__ defined and is set to %s\n",mkstr(__MINGW32__));
 #ifdef __USE_MINGW_ANSI_STDIO
  printf("__USE_MINGW_ANSI_STDIO is defined as %d\n",(int)(__USE_MINGW_ANSI_STDIO));
 #else
  printf("__USE_MINGW_ANSI_STDIO is NOT defined\n");
 #endif
 #if defined(__USE_MINGW_ANSI_STDIO) && __USE_MINGW_ANSI_STDIO==1
  printf("Using %s() for sprintf() & snprintf()\n","__mingw_...");
 #elif defined(_UCRT)
  printf("Using UCRT for sprintf() &  snprintf()\n");
 #else
  printf("Using %s() for sprintf() &  snprintf()\n","__ms_...");
 #endif
#else
 printf("__MINGW32__ is not defined\n");
#endif 
#ifdef __MINGW64__
 printf("__MINGW64__ defined and is set to %s\n",mkstr(__MINGW64__));
#endif
#ifdef _CRTBLD
 printf("_CRTBLD defined and is set to %s\n",mkstr(_CRTBLD));
#endif
#ifdef _UCRT
 printf("_UCRT defined and is set to %s\n",mkstr(_UCRT));
#endif
#ifdef __MSVCRT__
 printf("__MSVCRT__ defined and is set to %s\n",mkstr(__MSVCRT__));
#endif
#ifdef  _WIN32 /* true even on w64 */
 printf("_WIN32 defined and is set to %s\n",mkstr(_WIN32));
#endif
#ifdef  _WIN64 /* true even on w64 */
 printf("_WIN64 defined and is set to %s\n",mkstr(_WIN64));
#endif

#ifdef  __x86_64 /* true even on w64 */
 printf("__x86_64 defined and is set to %s\n",mkstr(__x86_64));
#endif


#if 1
 /* do more checks on long doubles that are normal numbers */
 // 355.0L/113.0L = PI to six decimal places - (relative error of about 8E−8)
 errs=0;
 printf("\nChecking long double functions:\n");
 long double testmld[]={-100.0L,-10.0L,-355.0L/113.0L,-1.0L,-0.5L,-0.1L,-0.01L,0.0L,0.01L,0.1L,0.5L,1.0L,355.0L/113.0L,10.0L,100.0L};
 for(int i=0;i<sizeof(testmld)/sizeof(testmld[0]);++i)
 	{
 	 check_mldbl(testmld[i]);
 	}
 printf(" long double functions: %u errors\n\n",errs);
#endif

#if 1
 printf("Checking fma functions:\n");
 check_fma(); // do some checks on fma(), fmal() amd fmaq()
 printf("\n");
#endif

 printf("Checking double double multiply and divide for accuracy\n");
 test_div();

 if(errs==0)
 	printf("All tests passed\n");
 else
    printf("%u tests failed\n",errs); 
#if defined(__BORLANDC__)
 fprintf(stderr,"Finished - press return to exit:");
 getchar();
#endif
 return 0;
}



/* now restore gcc options to those set by the user */
#if GCC_OPTIMIZE_AWARE
#pragma GCC pop_options
#endif