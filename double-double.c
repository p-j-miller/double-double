/* double-double.c
   ================
   
   double-double routines for basic doubles, long doubles (e.g. __float80) and __float128 types.
   
   The double-double representation represents numbers by a pair of native floating point types. This doubles the resolution (in mantissa bits or decimal digits),
   for example meaning that double-doubles based on the standard IEEE-754 64 bit double have 30 decimal digits (106 bits)
   compared to a single "native" double with 15 decimal digits (53bits). The exponent range is not changed.
   This approach was originally proposed by 
   D. E. Knuth in "The Art of Computer Programming, Volume II: Seminumerical Algorithms", 1st Edition, 1969, e.g. pp 237 and enhanced by amongst others 
   T. Dekker, "A floating-point technique for extending the available precision", Numerische Mathematik, 18(3) in 1971.
   See also IEEE trans. computers Vol 58 nos 7 July 2009 pp 994 "Accurate floating point product and exponentiation".
   Section 5.1 of the paper "Floating-Point Arithmetic" by Boldo, Jeannerod , Melquiond and Mukker (2023) which is at https://guillaume.melquiond.fr/doc/23-actanum.pdf ,
   gives an up to date status, and detailed error bounds for double word arithmetic of the type used here.
   
   Note TDM-GCC 9.2.2 fmaq() appears to be broken, see https://github.com/jmeubank/tdm-gcc/issues/14
    A working  fmaq() can be found at https://fossies.org/linux/gcc/libquadmath/math/fmaq.c  
	which is described as Member "gcc-9.3.0/libquadmath/math/fmaq.c" (12 Mar 2020, 9888 Bytes) of package /linux/misc/gcc-9.3.0.tar.xz:
	The corresponding header file was at https://raw.githubusercontent.com/gcc-mirror/gcc/master/libquadmath/quadmath-imp.h and 
	https://github.com/gcc-mirror/gcc/blob/master/libquadmath/quadmath-imp.h
	At least as of winlibs gcc 15.2.0 the supplied fmaq() works , but using the "local version" of fmaq() gives ~ 2* the execution speed on "typical" test programs so using it is still worthwhile!
   
   Version 1 of this file Peter Miller 25/5/2020 but parts go back a long way.
   Changes in this version (1/2026):
    - issue in f128_twosum() when both arguments are infinity fixed - I assume this effects other functions as well ?
    - issue in f128_mult_dd_dd() where it would return inf prematurely fixed
    - changed so long double function use "raw" long doubles rather than forcing them to be f80 (which is a type that's not defined in the C99 standard).
   	   - in some cases long doubles types are identical to double (i.e. 8 bytes), in others they are 12 or 16 bytes (perhaps with only 12 bytes actually used) - all of this is allowed by the C99 standard.
   	   - the code here should always work - but the resolution will obviously depend upon the underlying long double datatype.
    - u2_64toDD_f128() function added (works with u2_64.[ch] which gives a pseudo uint128 using uint64's)
    - suspected issue with fmal() - but nan_type.c test program disproved this
    - power10() functions added
    - multiply functions revised to improve their accuracy  (#define USE_ACCURATE_ALG)
    - divide algorithms revised to improve their accuracy (#define USE_ACCURATE_ALG)
    - ldd_power10() function optimised to reduce errors with large negative powers.

 main.c provides a test program.
 
The following functions are provided:
	// double-double functions based on the C double datatype gives approximately 30 significant decimal digits
	void add_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl);  // adds double double x and y to give double double "z"
	void sub_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl);  // subtracts double double y from x to give double double "z"
	void sub_dd_d(double *zh, double *zl,double xh, double xl,double y);  // subtracts double y from double double x to give double double "z"
	void mult_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl);  // multiplies double double a and b to give double double "x"
	void mult_d_dd( double *xh, double *xl,double a,double bh, double bl);  // multiplies a and double double  b to give double double "x"
	void dd_power(double *rh, double *rl,double x, unsigned int n);// raise x to nth power - return double double result , uses double double maths internally to minimise the error
	void dd_mult_power10( double *ohi, double *olo, double dh,double dl, int power );  // power can be +/-500. ASSUMES dh,dl came from a uint64 - this code will not work in general for all values dh,dl.
	void div_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl) ; // divides double double a by b to give double double "x"
	void U64toDD(uint64_t u64,double *xh,double *xl); // convert uint64 to double double 
	uint64_t ddtoU64(double xh,double xl);  // convert double double to uint64, deals with case where x may be a very small negative number

	// double-double based on two __float128 - gives approximately 66 significant decimal digits
	void f128_add_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl);  // adds double double x and y to give double double "z"
	void f128_sub_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl);  // subtracts double double y from x to give double double "z"
	void f128_mult_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl);  // multiplies double double a and b to give double double "x"
	void f128_mult_d_dd( f128_t *xh, f128_t *xl,f128_t a,f128_t bh, f128_t bl) ; // multiplies a and double double  b to give double double "x"
	void f128_dd_power(f128_t *rh, f128_t *rl,f128_t x, unsigned int n); // raise x to nth power - return double double result , uses double double maths internally to minimise the error
	void f128_power10(f128_t *rh, f128_t *rl , int n);// raise 10 to nth power (n can be positive or negative) - return f128_t double result , uses long double f128_t maths internally to give effectively perfect results - even for denormalised numbers
	void f128_mult_power10(__float128 *rh,__float128 *rl, __float128 rin_h, __float128 rin_l, int32_t rexp );  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
	void f128_div_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl);  // divides double double a by b to give double double "x"
	void u2_64toDD_f128(u2_64 u,f128_t *xh,f128_t *xl); // convert u2_64 (128 bits) to double double f128 
	void U128toDD_f128(uint128_t u128,f128_t *xh,f128_t *xl); // convert uint128 to double double f128 

	// double double based on two long doubles - resolution depends upon the underlying long double type, if this is "f80", then approximately 36 significant decimal digits are provided.
	void ld_twosum(long double *xh, long double *xl,long double a,long double b);  // adds a and b to give double double "x". 
	void add_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl);  // adds long double double x and y to give long double double "z"
	void sub_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl);  // subtracts long double double y from x to give long double double "z"
	void mult_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl);  // multiplies long double double a and b to give long double double "x"
	void mult_ld_ldd( long double *xh, long double *xl,long double a,long double bh, long double bl);  // multiplies long double a and long double double  b to give long double double "x"
	void ldd_power(long double *rh, long double *rl,long double x, unsigned int n); // raise x to nth power - return long double double result , uses long double double maths internally to minimise the error
	void ldd_power10(long double *rh, long double *rl , int n); // raise 10 to nth power (n can be positive or negative) - return long double double result , uses long double double maths internally to minimise the error
	void ldd_mult_power10(long double *rh,long double *rl,long double rin_h, long double rin_l, int32_t rexp );  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
	void div_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl);  // divides long double double a by b to give long double double "x"
	void U128toLDD(uint128_t u128,long double *xh,long double *xl); // convert uint128 to long double double 
 
  Note this code leverages (all available at https://github.com/p-j-miller )
	128 bit functionality using two uint64_t's from u2_64 
	power of 10 tables from power10
	my_printf for debugging from my_printf 
	
  It is recommended that these files (double-double.[ch]) are placed in a directory called double-double
*/   

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2020,2025,2026 Peter Miller
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/

#define USE_ACCURATE_ALG /* if defined use new, more accurate double-double multiplication and division code. Gives small accuracy improvements (but are a little slower) */ 


#define __USE_MINGW_ANSI_STDIO 1 /* So mingw uses its printf not msvcrt */

#define _FILE_OFFSET_BITS 64

#include <inttypes.h> /* defines PRI64 etc */
#include <limits.h>
#include <float.h> /* for limits for float, double , long double */
#include <stdio.h>
#include <time.h>      /* for time_t */
#include <sys/types.h> /* for off_t */
#include <string.h>
#include <stdint.h> /* for int32_t etc */
#include <stdbool.h> /* for bool, true, false */
#include <ctype.h> /* isspace() etc */
#include <math.h> /* fma() etc */


#if defined(__GNUC__) && defined(__SIZEOF_FLOAT128__)   && !defined(__BORLANDC__)      /* Builder C++ Community version 12.1 patch 1 defines __SIZEOF_FLOAT128__ but __float128's cannot be used in sensible programs due to compiler bugs */
  #include <quadmath.h> /* see https://gcc.gnu.org/onlinedocs/libquadmath/quadmath_005fsnprintf.html#quadmath_005fsnprintf - also needs quadmath library linking in - this is only available with gcc */
#endif


#include "../u2_64-128bits-with-two-u64/u2_64.h" /* this is a 128bit unsigned type that only used uint64_t types so needs no special compiler support */

#include "double-double.h"

#include "../power10/table_bin_10.h" /* "compressed table" created by power10.c - used by ldd_power10 & f128power10() */

#include "../my_printf/my_printf.h" /* printf that works around issues with standard print if necessary */

#define nos_elements_in(x) (sizeof(x)/(sizeof(x[0]))) /* number of elements in x , max index is 1 less than this as we index 0... */


/* double and double double functions 

*/


/* code below cannot be compiled with -Ofast as this makes the compiler break some C rules that we need, so make sure of this here */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) || defined(__clang__)
 #pragma GCC push_options
 #pragma GCC optimize ("-O3") /* cannot use Ofast, normally -O3 is OK. Note macro expansion does not work here ! */
 #if defined(_WIN32) && !defined(_WIN64)
  #pragma GCC target("sse2")
 #endif
#endif

static inline void fasttwosum(double *xh, double *xl,double a,double b)  // adds a and b to give double double "x". Requires |a| >= |b|
{// need to check for infinity as otherwise this function will return nan when given inf and inf as arguments 
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 double s=a+b;
 *xh=s;
 *xl=isfinite(s)? b-(s-a) : 0.0;
}
static inline void twosum(double *xh, double *xl,double a,double b)  // adds a and b to give double double "x". 
{// need to check for infinity as otherwise this function will return nan when given inf and inf as arguments 
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 double s=a+b;
 *xh=s;
 if(isfinite(s))
 	{double z=s-a;
 	 *xl=(a-(s-z))+(b-z);
 	}
 else
 	*xl=0;
}

 
static inline void twomult(double *xh, double *xl,double a,double b)  // multiplies  a and b to give double double "x"
{double x,y;
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 x=a*b;
 y=fma(a,b,-x);
 //if(isnan(y) || isnan(x)) fprintf(stderr,"twomult: %g*%g x=%g y=%g\n",a,b,x,y);
 *xh=x;
 *xl=y;
}

void add_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl)  // adds double double x and y to give double double "z"
{double sh,sl,th,tl,vh,vl,c,w;
twosum(&sh,&sl,xh,yh);
twosum(&th,&tl,xl,yl);
c=sl+th;
fasttwosum(&vh,&vl,sh,c);
w=tl+vl;
fasttwosum(zh,zl,vh,w);	
}

void sub_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl)  // subtracts double double y from x to give double double "z"
{
 add_dd_dd(zh,zl,xh,xl,-yh,-yl);
}

void sub_dd_d(double *zh, double *zl,double xh, double xl,double y)  // subtracts double y from double double x to give double double "z"
{double sh,sl,c;
 twosum(&sh,&sl,xh,-y);
 c=sl+xl;
 fasttwosum(zh,zl,sh,c);	
}


#ifdef USE_ACCURATE_ALG
// my new versions -  more accurate, but slower
void mult_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl)  // multiplies double double a and b to give double double "x"
{// ah*bh + ( ah*bl + al*bh ) + al*bl
 //  m1        m2      m3        m4
 // note the middle pair are similar order of magnitude and in the worse case can completely cancel 
 // e.g. (1+e)*(1-e) = 1-e^2  [ ie ( ah*bl + al*bh ) = e-e = 0 ]
 double t1,t2,t3;
 twomult(&t1,&t2,ah,bh); // m1
 t3=((ah*bl)+(al*bh));// (m2 + m3) 
 add_dd_dd(xh,xl,t1,t2,t3,al*bl);// add m4
 if(!isfinite(*xh))
	{
	 if( isfinite(fma(ah,bh,( ((ah*bl)+(al*bh)) +al*bl)) ))
	 	{ // trap incorrect overflow and return highest non infinite value [ can happen as ((ah*bl)+(al*bh))+al*bl can be negative ]
	 	  *xh=DBL_MAX;
	 	  *xl=0;
	 	  return;
	 	}
	 else if(isinf(t1))
	 	{*xh=t1;
	 	 *xl=0.0;
	 	}
	}
}
  
void mult_d_dd( double *xh, double *xl,double a,double bh, double bl)  // multiplies a and double double  b to give double double "x"
{// a*(bh+bl) = a*bh + a*bl
 //              m1     m2
 double t1,t2,t3;
 twomult(&t1,&t2,a,bh);// m1
 t3=(a*bl)+t2; // m2+ lower double from m1
 twosum(xh,xl,t1,t3);
}

#else /* original versions */
void mult_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl)  // multiplies double double a and b to give double double "x"
{double t1,t2,t3;
 twomult(&t1,&t2,ah,bh);
 //if(isnan(t1) || isnan(t2)) fprintf(stderr,"mult_dd_dd 1: (%g,%g)*(%g,%g) t1=%g t2=%g\n",ah,al,bh,bl,t1,t2); 
 t3=((ah*bl)+(al*bh))+t2;
 //if(isnan(t3) ) fprintf(stderr,"mult_dd_dd 2: (%g,%g)*(%g,%g) t1=%g t2=%g t3=%g\n",ah,al,bh,bl,t1,t2,t3); 
 twosum(xh,xl,t1,t3);
 //if(isnan(*xh) || isnan(*xl)) fprintf(stderr,"mult_dd_dd 3: (%g,%g)*(%g,%g) t1=%g t2=%g t3=%g xh=%g xl=%g\n",ah,al,bh,bl,t1,t2,t3,*xh,*xl); 
}

void mult_d_dd( double *xh, double *xl,double a,double bh, double bl)  // multiplies a and double double  b to give double double "x"
{double t1,t2,t3;
 twomult(&t1,&t2,a,bh);
 t3=(a*bl)+t2;
 twosum(xh,xl,t1,t3);
}
#endif

void dd_power(double *rh, double *rl,double x, unsigned int n)// raise x to nth power - return double double result , uses double double maths internally to minimise the error
/* even though this requires an initial (integer) loop to find lt, it is faster than starting from lsb as that needs two calls to mult_dd_dd() rather than one mult_dd_dd() and one mult_d_dd() */
{int lt=n,t;
 double hi=1.0,lo=0.0;
 while( (t=lt&(lt-1)) != 0) lt=t; // find msb of n
 // printf("dd_power(%g,%u [0x%x]): lt=%u [0x%x]\n",x,n,n,lt,lt);
 for(t=lt;t>0;t>>=1)
 	{// 1 bit at a time from msb
 	 mult_dd_dd(&hi,&lo,hi,lo,hi,lo);// r=r*r;
 	 if(t&n) mult_d_dd(&hi,&lo,x,hi,lo);//r=r*x;
 	}
 *rh=hi;
 if(!isfinite(hi)) *rl=0;
 else *rl=lo;
}

 /* use double double  for power 10 - note we don't use a compressed table for double-doubles as the full table is quite small.
    This function also avoids a divide for large negative powers (when we could have issues with denormalised mantissa's) by multiplying by a power of 2, then multiplying by 10-x then dividing by a power of 2 (which is fast and exact).
 */

void dd_mult_power10( double *ohi, double *olo, double dh,double dl, int power )  // power can be +/-500. ASSUMES dh,dl came from a uint64 - this code will not work in general for all values dh,dl.
{ 
  double th,tl; // power10
  double xh,xl;
  if(power==0)
  	{*ohi=dh; /* 10^0=1 so can just return d - we expect this special case to be relatively common so worth doing */
  	 *olo=dl;
    }
  else if(power<-200) // need to do in 2 multiplies rather than 1 as the max -ve exponent for a normalised double is -308. Picked -200 so we keep full precision of double double
    {if((-power)-200 >= sizeof(NegPowerOf10_hi)/sizeof(NegPowerOf10_hi[0]))
    	{*ohi=0.0; // underflow
    	 *olo=0.0;
    	 return;
    	}
	 th=NegPowerOf10_hi[200];
	 tl=NegPowerOf10_lo[200];   
     mult_dd_dd(&xh,&xl,th,tl,dh,dl);// xh/l=th/l*d
     // there is a risk that the next multiply will end up with denormalised numbers and so can be inexact, so multiply by 2^100 first, do multiply then divide back again
     // we have just divided by 10^200 so the multiply by 2^100 (1.3e30) cannot overflow
     xh=ldexp(xh,100);// multiply by 2^100
     xl=ldexp(xl,100);
	 th=NegPowerOf10_hi[(-power)-200];
	 tl=NegPowerOf10_lo[(-power)-200];  
	 mult_dd_dd(&xh,&xl,th,tl,xh,xl);    	 
     xh=ldexp(xh,-100);// multiply by 2^-100
     xl=ldexp(xl,-100);
	 *ohi=xh;
	 *olo=xl;	 
  	}
  else if(power<0) 	
  	{// can do with one double-double multiply
	 th=NegPowerOf10_hi[-power];
	 tl=NegPowerOf10_lo[-power];
	 mult_dd_dd(ohi,olo,th,tl,dh,dl);
  	}
	
#define _ATOF_MAXPOWER1 288 /* when getting near the max (e308) split into 2 to avoid an overflow. Test suite fails at 290, is OK for 287 to 289 and fails again at 286 - so set to 288 */
  else if (power>_ATOF_MAXPOWER1) // need to do in 2 multiplies rather than 1 as the max exponent for a double is 308 and we may need to multiply by 10^350 here. Picked 290 as thats 18 sf from 308 so keeping full precision of double double
  	{if(power-_ATOF_MAXPOWER1 >= sizeof(PosPowerOf10_hi)/sizeof(PosPowerOf10_hi[0]))
    	{*ohi=INFINITY; // overflow
    	 *olo=0.0;
    	 return;
    	}
	 th=PosPowerOf10_hi[_ATOF_MAXPOWER1];
	 tl=PosPowerOf10_lo[_ATOF_MAXPOWER1];   
     mult_dd_dd(&xh,&xl,th,tl,dh,dl);// xh/l=th/l*d
	 th=PosPowerOf10_hi[power-_ATOF_MAXPOWER1];
	 tl=PosPowerOf10_lo[power-_ATOF_MAXPOWER1];  
	 mult_dd_dd(ohi,olo,th,tl,xh,xl);   
	 if(!isfinite(*ohi) || !isfinite(*olo))
	 	{// result overflows - return infinity (note they could be set to NAN at this point)
	 	 *ohi=INFINITY;
	 	 *olo=0.0;
	 	}	 	 
  	}
  else
    {// can do with one double-double multiply
	 th=PosPowerOf10_hi[power];
	 tl=PosPowerOf10_lo[power];	
	 mult_dd_dd(ohi,olo,th,tl,dh,dl);
	}	
}

#ifdef USE_ACCURATE_ALG
#if defined(__SIZEOF_FLOAT128__)   && !defined(__BORLANDC__)      /* Builder C++ Community version 12.1 patch 1 defines __SIZEOF_FLOAT128__ but __float128's cannot be used in sensible programs due to compiler bugs */
void div_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl)  // divides double double a by b to give double double "x" 
{ // if we have __float128 then there is a simple implementation thats probably faster than the code below as the mantissa of a __float has enough bits to accurately hold a double-double, and the float128 exponent range is much larger
  // this gives identical results to the implementation below with the test program.
 __float128 r;
 double rh,rl;
 r=((__float128)ah+(__float128)al)/((__float128)bh+(__float128)bl);// do the actual divide as float128
 rh=r;
 if(isnormal(rh))
 	rl=r-rh;
 else rl=0.0;
 *xh=rh;
 *xl=rl;
}
#else
/* 2 iterations of 2nd order iteration */
/* results :
	 355/113 (approximation to pi):
	 as dd =3.14159292035398252e+00,-2.20079608421269966e-16
	 dd's combined to f128's=3.141592920353982300884955752212392e+00
	               as f128's=3.141592920353982300884955752212389e+00
	              Difference=3.081487911019577364889564708135884e-33
	 (355+e)/(113) (e=double epsilon):
	 as dd =3.14159292035398252e+00,-2.18114611917508624e-16
	 dd's combined to f128's=3.141592920353982302849952255973734e+00
	               as f128's=3.141592920353982302849952255973728e+00
	              Difference=6.162975822039154729779129416271767e-33
	 (355+e)/(113+e) (e=double epsilon):
	 as dd =3.14159292035398208e+00,2.19801378827817027e-16
	 dd's combined to f128's=3.141592920353982296676733151236769e+00
	               as f128's=3.141592920353982296676733151236778e+00
	              Difference=-9.629649721936179265279889712924637e-33
	 (355+e/1e+14)/(113+e/1e+14) (e=double epsilon):
	 as dd =3.14159292035398252e+00,-2.20079608421270015e-16
	 dd's combined to f128's=3.141592920353982300884955752212343e+00
	               as f128's=3.141592920353982300884955752212347e+00
	              Difference=-4.237045877651918876723151473686840e-33
*/
               
void div_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl)  // divides double double a by b to give double double "x" - this is  slightly more accurate than the basic version below
{ // two *, 3 /
 double q1, q2, q3;
 double rh,rl,r1h,r1l;  
 if(!isfinite(ah))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!isfinite(bh))
 	{if(isinf(bh)) 
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	} 
 q1=ah/bh; // initial estimate
 mult_d_dd(&rh,&rl,q1,bh,bl); // r=q1*b
 sub_dd_dd(&rh,&rl,ah,al,rh,rl);// r=a-r ie r=a-q1*b
 q2=rh/bh;
 mult_d_dd(&r1h,&r1l,q2,bh,bl); // r1=q2*b
 sub_dd_dd(&rh,&rl,rh,rl,r1h,r1l);// r=r-(q2*b)
 q3=rh/bh;
 fasttwosum(&q1,&q2,q1,q2);
 add_dd_dd(xh,xl,q1,q2,q3,0.0);
} 
#endif

#else /* original version - less accurate */
/* results:
	 355/113 (approximation to pi):
	 as dd =3.14159292035398252e+00,-2.20079608421269966e-16
	 dd's combined to f128's=3.141592920353982300884955752212392e+00
	               as f128's=3.141592920353982300884955752212389e+00
	              Difference=3.081487911019577364889564708135884e-33 => same as above
	 (355+e)/(113) (e=double epsilon):
	 as dd =3.14159292035398252e+00,-2.18114611917508624e-16
	 dd's combined to f128's=3.141592920353982302849952255973734e+00
	               as f128's=3.141592920353982302849952255973728e+00
	              Difference=6.162975822039154729779129416271767e-33 => same as above
	 (355+e)/(113+e) (e=double epsilon):
	 as dd =3.14159292035398208e+00,2.19801378827817027e-16
	 dd's combined to f128's=3.141592920353982296676733151236769e+00
	               as f128's=3.141592920353982296676733151236778e+00
	              Difference=-9.629649721936179265279889712924637e-33 => same as above
	 (355+e/1e+14)/(113+e/1e+14) (e=double epsilon):
	 as dd =3.14159292035398252e+00,-2.20079608421269991e-16
	 dd's combined to f128's=3.141592920353982300884955752212368e+00
	               as f128's=3.141592920353982300884955752212347e+00
	              Difference=2.041485741050470004239336619140023e-32 => this is much higher than those above (5*) which were -4.237045877651918876723151473686840e-33

*/
void div_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl)  // divides double double a by b to give double double "x"
{ // 3 *, 2 /
 double hi;
 double ch,cl,uh,ul;
  if(!isfinite(ah))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!isfinite(bh))
 	{if(isinf(bh)) 
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	} 
 ch=ah/bh; // initial estimate
 uh=ch*bh;
 ul=fma(ch,bh,-uh);
 cl=(((ah-uh)-ul)+al-ch*bl)/bh;
 hi=ch+cl;
 fasttwosum(xh,xl,hi,cl+(ch-hi)); /* ensure result is "normalised" */
}
#endif

void U64toDD(uint64_t u64,double *xh,double *xl) /* convert uint64 to double double */
{uint64_t u64_l=u64 & UINT64_C(0xffffffff );// bottom 32 bits [ double has a 52 bit mantissa so 32 bits can fit exactlty  ]
 uint64_t u64_h=u64^u64_l;// upper 32 bits
 if(u64_h==0)
 	{*xh=u64; // <= 32 bits so will fit into one double
 	 *xl=0;
 	 return;
 	}
 twosum(xh,xl,u64_h,u64_l);// combine two halves of u64 to make the double double could probably be fasttwosum() but lets be sure	
}


uint64_t ddtoU64(double xh,double xl)  // convert double double to uint64, deal with case where x may be a very small negative number
{ 
  double t;
  uint64_t ob;
  if(xh<0.0)
  	{t=rint(xh+xl);
  	 if(t<0) return 0;
  	  return t;
  	}
  ob = (uint64_t)xh;	
  t = ( xh - (double)ob );
  t += xl;
  ob += (uint64_t)rint(t);
  return ob;
}

#if defined(__SIZEOF_FLOAT128__)   && !defined(__BORLANDC__)      /* Builder C++ Community version 12.1 patch 1 defines __SIZEOF_FLOAT128__ but __float128's cannot be used in sensible programs due to compiler bugs */

/* double double functions using flt128 so in theory give ~ 66 significant digits.

   
*/

/* following are predefined by gcc
__FLT128_MAX_10_EXP__ 4932
__FLT128_DENORM_MIN__ 6.47517511943802511092443895822764655e-4966F128
__FLT128_MIN_EXP__ (-16381)
__FLT128_MIN_10_EXP__ (-4931)
__FLT128_MANT_DIG__ 113
__FLT128_HAS_INFINITY__ 1
__FLT128_MAX_EXP__ 16384
__FLT128_HAS_DENORM__ 1
__FLT128_DIG__ 33
__FLT128_MAX__ 1.18973149535723176508575932662800702e+4932F128
__SIZEOF_FLOAT128__ 16
__FLT128_MIN__ 3.36210314311209350626267781732175260e-4932F128
__FLT128_HAS_QUIET_NAN__ 1
__FLT128_EPSILON__ 1.92592994438723585305597794258492732e-34F128
__FLT128_DECIMAL_DIG__ 36


__SIZEOF_INT128__ 16

These should be referenced as :

	FLT128_MAX: largest finite number
	FLT128_MIN: smallest positive number with full precision
	FLT128_EPSILON: difference between 1 and the next larger representable number 
	FLT128_DENORM_MIN: smallest positive denormalized number
	FLT128_MANT_DIG: number of digits in the mantissa (bit precision)
	FLT128_MIN_EXP: maximal negative exponent
	FLT128_MAX_EXP: maximal positive exponent
	FLT128_DIG: number of decimal digits in the mantissa
	FLT128_MIN_10_EXP: maximal negative decimal exponent
	FLT128_MAX_10_EXP: maximal positive decimal exponent
	
*/


static inline void f128_fasttwosum(f128_t *xh, f128_t *xl,f128_t a,f128_t b)  // adds a and b to give double double "x". Requires |a| >= |b|
{if(isinfq(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinfq(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 f128_t s=a+b;
 *xh=s;
 *xl=finiteq(s)? b-(s-a) : 0.0;
}

static inline void f128_twosum(f128_t *xh, f128_t *xl,f128_t a,f128_t b)  // adds a and b to give double double "x". 
{if(isinfq(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinfq(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 f128_t s=a+b;
 *xh=s;
 if(finiteq(s))
 	{f128_t z=s-a;
 	 *xl=(a-(s-z))+(b-z);
 	}
 else
 	*xl=0;
}


static inline void f128_twomult(f128_t *xh, f128_t *xl,f128_t a,f128_t b)  // multiplies  a and b to give double double "x"
{f128_t x,y;
 if(isinfq(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinfq(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 x=a*b;
 y=fmaq(a,b,-x);
 *xh=x;
 *xl=y;
}


 void f128_add_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl)  // adds double double x and y to give double double "z"
{f128_t sh,sl,th,tl,vh,vl,c,w;
f128_twosum(&sh,&sl,xh,yh);
f128_twosum(&th,&tl,xl,yl);
c=sl+th;
f128_fasttwosum(&vh,&vl,sh,c);
w=tl+vl;
f128_fasttwosum(zh,zl,vh,w);	
}

 void f128_sub_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl)  // subtracts double double y from x to give double double "z"
{
 f128_add_dd_dd(zh,zl,xh,xl,-yh,-yl);
}

#ifdef USE_ACCURATE_ALG
// my new versions - are more accurate
void f128_mult_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl)  // multiplies double double a and b to give double double "x"
{// ah*bh + ( ah*bl + al*bh ) + al*bl
 //  m1        m2      m3        m4
 // note the middle pair are similar order of magnitude and in the worse case can completely cancel 
 // e.g. (1+e)*(1-e) = 1-e^2  [ ie ( ah*bl + al*bh ) = e-e = 0 ]
 f128_t t1,t2,t3;
 f128_twomult(&t1,&t2,ah,bh); // m1
 t3=((ah*bl)+(al*bh));// (m2 + m3) 
 f128_add_dd_dd(xh,xl,t1,t2,t3,al*bl);// add m4

 if(!finiteq(*xh))
	{
	 if( finiteq(fmaq(ah,bh,( ((ah*bl)+(al*bh)) +al*bl)) ))
	 	{ // trap incorrect overflow and return highest non infinite value [ can happen as ((ah*bl)+(al*bh)) +al*bl can be negative ]
	 	  // printf("\n\n f128_mult_dd_dd t1=inf, t3=inf but ah*bh != inf\n");
	 	  *xh=FLT128_MAX;
	 	  *xl=0;
	 	  return;
	 	}
	 else if(isinfq(t1))
	 	{*xh=t1;
	 	 *xl=0.0;
	 	}
	}
}

void f128_mult_d_dd( f128_t *xh, f128_t *xl,f128_t a,f128_t bh, f128_t bl)  // multiplies a and double double  b to give double double "x"
{// a*(bh+bl) = a*bh + a*bl
 //              m1     m2
 f128_t t1,t2,t3;
 f128_twomult(&t1,&t2,a,bh);// m1
 t3=(a*bl)+t2;// m2+ lower double from m1
 f128_twosum(xh,xl,t1,t3);
}

#else
void f128_mult_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl)  // multiplies double double a and b to give double double "x"
{f128_t t1,t2,t3;
 f128_twomult(&t1,&t2,ah,bh);
 t3=((ah*bl)+(al*bh))+t2;
 if(isinfq(t1) && isinfq(t3) && finiteq(fmaq(ah,bh,( ((ah*bl)+(al*bh)) +al*bl)) ))
 	{ // trap incorrect overflow and return highest non infinite value [ can happen as ((ah*bl)+(al*bh)) can be negative ]
 	  // printf("\n\n f128_mult_dd_dd t1=inf, t3=inf but ah*bh != inf\n");
 	  *xh=FLT128_MAX;
 	  *xl=0;
 	  return;
 	}
 f128_twosum(xh,xl,t1,t3);
}

 void f128_mult_d_dd( f128_t *xh, f128_t *xl,f128_t a,f128_t bh, f128_t bl)  // multiplies a and double double  b to give double double "x"
{f128_t t1,t2,t3;
 f128_twomult(&t1,&t2,a,bh);
 t3=(a*bl)+t2;
 f128_twosum(xh,xl,t1,t3);
}
#endif

 void f128_dd_power(f128_t *rh, f128_t *rl,f128_t x, unsigned int n)// raise x to nth power - return double double result , uses double double maths internally to minimise the error
{int lt=n,t;
 f128_t hi=1.0Q,lo=0.0Q;
 while( (t=lt&(lt-1)) != 0) lt=t; // find msb of n
 // printf("dd_power(%g,%u [0x%x]): lt=%u [0x%x]\n",x,n,n,lt,lt);
 for(t=lt;t>0;t>>=1)
 	{// 1 bit at a time from msb
 	 f128_mult_dd_dd(&hi,&lo,hi,lo,hi,lo);// r=r*r;
 	 if(t&n) f128_mult_d_dd(&hi,&lo,x,hi,lo);//r=r*x;
 	}
 *rh=hi;
 if(!finiteq(hi)) *rl=0;
 else *rl=lo;
}


void f128_power10(f128_t *rh, f128_t *rl , int n)// raise 10 to nth power (n can be positive or negative) - return f128_t double result , uses long double f128_t maths internally to give effectively perfect results - even for denormalised numbers
{int offset=0,index;
bool first=true;
 if(n==0) {*rh=1.0Q; *rl=0.0Q;}
 else if(n>0)
 	{
	 if(n>FLT128_MAX_10_EXP)
	 	{*rh=INFINITY;// avoids having a special check below for going of the end of the table max 10 exponent for long doubles is 4932 so this function will return infinity well before 8191 (the limit from the table structure)
	 	 *rl=0.0;// this matches what we have used elsewhere (in double-double.c) for overflow of "double-double"
		 return;
		}
	 while(n)
	 	{
	 	 index=offset+(((n)&0xf)<<1);
	 	 if(n&0xf) 
		  	{if(first) // just use assign for 1st value, then need to use multiply
				{first=false;
				 *rh=F128_DD_pos_pow10_compressed[index]; 
	 			 *rl=F128_DD_pos_pow10_compressed[index+1]; 
	 			}
			 else 
				{f128_mult_dd_dd(rh,rl,*rh,*rl,F128_DD_pos_pow10_compressed[index],F128_DD_pos_pow10_compressed[index+1]); // only multiply if non zero */
				}
			}
	 	 n>>=4;
	 	 offset+=30; /* 2 values (hi & lo) per "entry" */
	 	}
	}
 else
 	{// n< 0
	 if(n<-8191)
	 	{*rh=*rl=0;// avoids having a special check below for going of the end of the table max 10 exponent for long doubles is 4932 so this function will return 0 well before -8191
		 return;
		}
	 n= -n; // make n positive, the remainder of this function is identical to the +ve case - but using the negative table rather than the positive table.
	 while(n)
	 	{
	 	 index=offset+(((n)&0xf)<<1);
	 	 if(n&0xf) 
		  	{if(first) // just use assign for 1st value, then need to use multiply
				{first=false;
				 *rh=F128_DD_neg_pow10_compressed[index]; 
	 			 *rl=F128_DD_neg_pow10_compressed[index+1]; 
	 			}
			 else 
#if 0
				{f128_div_dd_dd(rh,rl,*rh,*rl,F128_DD_pos_pow10_compressed[index],F128_DD_pos_pow10_compressed[index+1]); // only divide if non zero, divide is probably more accurate - but slower */
				}
#else		 
				{f128_mult_dd_dd(rh,rl,*rh,*rl,F128_DD_neg_pow10_compressed[index],F128_DD_neg_pow10_compressed[index+1]); // only multiply if non zero */
				}
#endif				
			}
	 	 n>>=4;
	 	 offset+=30; /* 2 values (hi & lo) per "entry" */
	 	}
	} 
 return;
}

 void f128_mult_power10(__float128 *rh,__float128 *rl, __float128 rin_h, __float128 rin_l, int32_t rexp )  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
 {
  *rh=rin_h; *rl=rin_l; 	
  if(rexp>0)
  	{
 	 // here we allow multiplication by up to 10^2*maxExponent which is by far enough
	 int exp=rexp;	 
	 if(rexp>4096) // 4096 is both "exact" and minimises the mumber of multiplies required
	 	{
 		 rexp=4096;
 		 exp-=4096; // any excess which we will also need to divide by (if its > 0)
 		}
 	  else exp=0;  
	  f128_t tenh,tenl; 	
 	  f128_power10(&tenh,&tenl,rexp); // 10^rexp	     	  
 	  f128_mult_dd_dd(rh,rl,*rh,*rl,tenh,tenl); 
 	  if(exp>0)
 	  	{  		   	  		
		 f128_power10(&tenh,&tenl,exp); // 10^exp					  
 	  	 f128_mult_dd_dd(rh,rl,*rh,*rl,tenh,tenl); 
 	    } 	    		 
	}
 else if(rexp<0)
 	{// need to take care here as mantissa is > 1 so even dividing by 10^maxExponent may not be enough, here we allow division by up to 10^2*maxExponent which is by far enough
	 rexp= -rexp;
	 int exp=rexp;	 
	 if(rexp>4096) // 4096 is both "exact" and minimises the mumber of multiplies required
	 	{
 		 rexp=4096;
 		 exp-=4096; // any excess which we will also need to divide by (if its > 0)
 		}
 	  else exp=0;  
	  f128_t tenh,tenl; 	
 	  f128_power10(&tenh,&tenl,-rexp); // 10^-rexp	     	  
 	  f128_mult_dd_dd(rh,rl,*rh,*rl,tenh,tenl); // negative exponent means we divide by powers of 10
 	  if(exp>0)
 	  	{  		   	  		
		 f128_power10(&tenh,&tenl,exp); // 10^exp					  
 	  	 f128_div_dd_dd(rh,rl,*rh,*rl,tenh,tenl); // negative exponent means we divide by powers of 10 - so divide by the rest of the exponent - this must be a divide for accuracy.
 	    } 	    
	}	
 // else special case, rexp==0 - nothing to do
}

#ifdef USE_ACCURATE_ALG
/* my implementation based on 2nd order formula for reciprocal in "Multiplicative Iteration for Reciprocals", Prof. W. Kahan, https://people.eecs.berkeley.edu/~wkahan/CS279/recip.pdf 
   for r=1/x a more accurate approximation R is given by:
   d=(1-x*r)
   R=d*r+r
   2 iterations of 2nd order - arranged so we approximate a/b rather than calculate 1/b then multiply by a (so d=a-x*r)
*/
void f128_div_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl)  // divides double double a by b to give double double "x"
{
 f128_t q1, q2, q3;
 f128_t rh,rl,r1h,r1l; 
 if(!finiteq(ah) || (ah==0.0Q && al==0.0Q))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!finiteq(bh))
 	{if(isinfq(bh))  
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	}
 q1=ah/bh; // initial estimate
 f128_mult_d_dd(&rh,&rl,q1,bh,bl); // r=q1*b
 f128_sub_dd_dd(&rh,&rl,ah,al,rh,rl);// r=a-r ie r=a-q1*b
 q2=rh/bh;
 f128_mult_d_dd(&r1h,&r1l,q2,bh,bl); // r1=q2*b
 f128_sub_dd_dd(&rh,&rl,rh,rl,r1h,r1l);// r=r-(q2*b)
 q3=rh/bh;
 f128_fasttwosum(&q1,&q2,q1,q2);
 f128_add_dd_dd(xh,xl,q1,q2,q3,0.0);
}	

#else
void f128_div_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl)  // divides double double a by b to give double double "x"
{f128_t hi;
 f128_t ch,cl,uh,ul;
 if(!finiteq(ah))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!finiteq(bh))
 	{if(isinfq(bh))  
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	}
 ch=ah/bh; // initial estimate
 uh=ch*bh;
 ul=fmaq(ch,bh,-uh);
 cl=(((ah-uh)-ul)+al-ch*bl)/bh;
 hi=ch+cl;
 f128_fasttwosum(xh,xl,hi,cl+(ch-hi)); /* ensure result is "normalised" */ 
}
#endif

#if defined( __SIZEOF_INT128__)  /* only allow if compiler also supports  __int128 */
  void U128toDD_f128(uint128_t u128,f128_t *xh,f128_t *xl) /* convert uint128 to double double f128 */
{uint128_t u128_l=u128 & UINT64_C(0xffffffffffffffff );// bottom 64 bits [ f128 has a 112 bit mantissa so 64 bits can fit exactlty  ]
 uint128_t u128_h=u128^u128_l;// upper 64 bits
 if(u128_h==0)
 	{*xh=u128; // <= 64 bits so will fit into one f128
 	 *xl=0;
 	 return;
 	}
 f128_twosum(xh,xl,u128_h,u128_l);// combine two halves of u128 to make the double double could probably be f128_fasttwosum() but lets be sure	
}
#endif

void u2_64toDD_f128(u2_64 u,f128_t *xh,f128_t *xl) /* convert u2_64 (128 bits) to double double f128 */
{
 if(u.hi==0)
 	{*xh=u.lo; // <= 64 bits so will fit into one f128
 	 *xl=0;
 	 return;
 	}
 f128_t fh=u2_64_to_flt128(u64_to_u2_64(u.hi,0));
 f128_t fl=u2_64_to_flt128(u64_to_u2_64(0,u.lo));
 f128_twosum(xh,xl,fh,fl);// combine two halves of u2_64 to make the double double could probably be f128_fasttwosum() but lets be sure	
}

#endif


/* double double functions using long double so in theory give ~ 36 significant digits.  
*/
static inline void ld_fasttwosum(long double *xh, long double *xl,long double a,long double b)  // adds a and b to give double double "x". Requires |a| >= |b|
{ // need to check for infinity as otherwise this function will return nan when given inf and inf as arguments 
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 long double s=a+b;
 *xh=s;
 *xl=isfinite(s)? b-(s-a) : 0.0;
}

void ld_twosum(long double *xh, long double *xl,long double a,long double b)  // adds a and b to give double double "x". 
{ // need to check for infinity as otherwise this function will return nan when given inf and inf as arguments 
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 long double s=a+b;
 *xh=s;
 if(isfinite(s))
 	{long double z=s-a;
 	 *xl=(a-(s-z))+(b-z);
 	}
 else
 	*xl=0;
}

static inline void ld_twomult(long double *xh, long double *xl,long double a,long double b)  // multiplies  a and b to give double double "x"
{long double x,y;
 if(isinf(a))
 	{*xh=a;
 	 *xl=0;
 	 return;
 	}
 else if(isinf(b))
 	{*xh=b;
	 *xl=0;
	 return;
	}
 x=a*b;
 y=fmal(a,b,-x);
 *xh=x;
 *xl=y;
}


void add_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl)  // adds long double double x and y to give long double double "z"
{long double sh,sl,th,tl,vh,vl,c,w;
ld_twosum(&sh,&sl,xh,yh);
ld_twosum(&th,&tl,xl,yl);
c=sl+th;
ld_fasttwosum(&vh,&vl,sh,c);
w=tl+vl;
ld_fasttwosum(zh,zl,vh,w);	
}

void sub_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl)  // subtracts long double double y from x to give long double double "z"
{
 add_ldd_ldd(zh,zl,xh,xl,-yh,-yl);
}

#ifdef USE_ACCURATE_ALG
// my new versions - are more accurate, but slower
void mult_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl)  // multiplies long double double a and b to give long double double "x"
{// ah*bh + ( ah*bl + al*bh ) + al*bl
 //  m1        m2      m3        m4
 // note the middle pair are similar order of magnitude and in the worse case can completely cancel 
 // e.g. (1+e)*(1-e) = 1-e^2  [ ie ( ah*bl + al*bh ) = e-e = 0 ]
 long double t1,t2,t3;
 ld_twomult(&t1,&t2,ah,bh); // m1
 t3=((ah*bl)+(al*bh));// (m2 + m3) 
 add_ldd_ldd(xh,xl,t1,t2,t3,al*bl);// add m4 

 if(!isfinite(*xh))
	{
	 if( isfinite(fmal(ah,bh, ( ((ah*bl)+(al*bh)) +al*bl)) ))
	 	{ // trap incorrect overflow and return highest non infinite value [ can happen as ((ah*bl)+(al*bh))+al*bl can be negative ]
	 	  *xh=LDBL_MAX;
	 	  *xl=0;
	 	  return;
	 	}
	 else if(isinf(t1))
	 	{*xh=t1;
	 	 //my_printf(" t1 infinite ( ((ah*bl)+(al*bh)) +al*bl)) )=%Lg\n", ((ah*bl)+(al*bh)) +al*bl); // t3 is nan if t1 is inf
	 	 *xl=0.0;
	 	}
	}
}

void mult_ld_ldd( long double *xh, long double *xl,long double a,long double bh, long double bl)  // multiplies long double a and long double double  b to give long double double "x"
{// a*(bh+bl) = a*bh + a*bl
 //              m1     m2
 long double t1,t2,t3;
 ld_twomult(&t1,&t2,a,bh);// m1
 t3=(a*bl)+t2; // m2+ lower double from m1
 ld_twosum(xh,xl,t1,t3);
}

#else
void mult_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl)  // multiplies long double double a and b to give long double double "x"
{long double t1,t2,t3;
 ld_twomult(&t1,&t2,ah,bh);
 t3=((ah*bl)+(al*bh))+t2;
 ld_twosum(xh,xl,t1,t3);
 //if(isinf(t1) && isinf(t3) && !isinf(fmal(ah,bh,((ah*bl)+(al*bh)))))
 if(!isfinite(*xh))
	{
	 if( isfinite(fmal(ah,bh, (ah*bl)+(al*bh)) ))
	 	{ // trap incorrect overflow and return highest non infinite value [ can happen as ((ah*bl)+(al*bh))+al*bl can be negative ]
	 	  *xh=LDBL_MAX;
	 	  *xl=0;
	 	  return;
	 	}
	 else if(isinf(t1))
	 	{*xh=t1;
	 	 //my_printf(" t1 infinite ( ((ah*bl)+(al*bh)) +al*bl)) )=%Lg\n", ((ah*bl)+(al*bh)) +al*bl); // t3 is nan if t1 is inf
	 	 *xl=0.0;
	 	}
	}
}

void mult_ld_ldd( long double *xh, long double *xl,long double a,long double bh, long double bl)  // multiplies long double a and long double double  b to give long double double "x"
{long double t1,t2,t3;
 ld_twomult(&t1,&t2,a,bh);
 t3=(a*bl)+t2;
 ld_twosum(xh,xl,t1,t3);
}
#endif

void ldd_power(long double *rh, long double *rl,long double x, unsigned int n)// raise x to nth power - return long double double result , uses long double double maths internally to minimise the error
{int lt=n,t;
 long double hi=1.0L,lo=0.0L;
 while( (t=lt&(lt-1)) != 0) lt=t; // find msb of n
 // printf("dd_power(%g,%u [0x%x]): lt=%u [0x%x]\n",x,n,n,lt,lt);
 for(t=lt;t>0;t>>=1)
 	{// 1 bit at a time from msb
 	 mult_ldd_ldd(&hi,&lo,hi,lo,hi,lo);// r=r*r;
 	 if(t&n) mult_ld_ldd(&hi,&lo,x,hi,lo);//r=r*x;
 	}
 *rh=hi;
 if(!isfinite(hi)) *rl=0;
 else *rl=lo;
}


#if LDBL_MAX_10_EXP==4932 /* these two power10 functions are only available for "true" long doubles ) if long double=double there are replacement functions below (#else) */
void ldd_power10(long double *rh, long double *rl , int n)// raise 10 to nth power (n can be positive or negative) - return long double double result , uses long double double maths internally to give effectively perfect results - even for denormalised numbers
{int offset=0,index;
bool first=true;
 if(n==0) {*rh=1.0L; *rl=0.0L;}
 else if(n>0)
 	{
	 if(n>LDBL_MAX_10_EXP)
	 	{*rh=INFINITY;// avoids having a special check below for going of the end of the table max 10 exponent for long doubles is 4932 so this function will return infinity well before 8191
	 	 *rl=0.0;// this matches what we have used elsewhere (in double-double.c) for overflow of "double-double"
		 return;
		}
	 while(n)
	 	{
	 	 index=offset+(((n)&0xf)<<1);
	 	 if(n&0xf) 
		  	{if(first) // just use assign for 1st value, then need to use multiply
				{first=false;
				 *rh=LD_DD_pos_pow10_compressed[index]; 
	 			 *rl=LD_DD_pos_pow10_compressed[index+1]; 
	 			}
			 else 
				{mult_ldd_ldd(rh,rl,*rh,*rl,LD_DD_pos_pow10_compressed[index],LD_DD_pos_pow10_compressed[index+1]); // only multiply if non zero */
				}
			}
	 	 n>>=4;
	 	 offset+=30; /* 2 values (hi & lo) per "entry" */
	 	}
	}
 else
 	{// n< 0
	 if(n<-8191)
	 	{*rh=*rl=0;// avoids having a special check below for going of the end of the table max 10 exponent for long doubles is 4932 so this function will return 0 well before -8191
		 return;
		}
	 n= -n; // make n positive, the remainder of this function is identical to the +ve case - but using the negative table rather than the positive table.
	 while(n)
	 	{
	 	 index=offset+(((n)&0xf)<<1);
	 	 if(n&0xf) 
		  	{if(first) // just use assign for 1st value, then need to use multiply
				{first=false;
				 *rh=LD_DD_neg_pow10_compressed[index]; 
	 			 *rl=LD_DD_neg_pow10_compressed[index+1]; 
	 			}
			 else 
				{div_ldd_ldd(rh,rl,*rh,*rl,LD_DD_pos_pow10_compressed[index],LD_DD_pos_pow10_compressed[index+1]); // only divide (by +ve power10) if non zero - divide is more accurate than multiply by negative power10 
				}				
			}
	 	 n>>=4;
	 	 offset+=30; /* 2 values (hi & lo) per "entry" */
	 	}
	} 
 return;
}

 /*  calculate dr=double-double (long double)r*powl(10,rexp), using double double long double maths for accuracy  */
void ldd_mult_power10(long double *rh,long double *rl,long double rin_h, long double rin_l, int32_t rexp )  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
 {
  *rh=rin_h; *rl=rin_l; 
  if(rexp>0)
  	{
 	 // here we allow multiplication by up to 10^2*maxExponent which is by far enough
	 int exp=rexp;	 
	 if(rexp>4096) /* 10^4096 using ldd_power10() is "exact" (and fast) so this is a good one to duplicate */
	 	{
 		 rexp=4096;
 		 exp-=4096; // any excess which we will also need to multiply by (if its > 0)
 		}
 	  else exp=0;  
	  long double tenh,tenl; 	
 	  ldd_power10(&tenh,&tenl,rexp); // 10^rexp	     	  
 	  mult_ldd_ldd(rh,rl,*rh,*rl,tenh,tenl); 
 	  if(exp>0)
 	  	{  		   	  		
		 ldd_power10(&tenh,&tenl,exp); // 10^exp					  
 	  	 mult_ldd_ldd(rh,rl,*rh,*rl,tenh,tenl); 
 	    } 	    		 
	}
 else if(rexp<0)
 	{// need to take care here as mantissa is > 1 so even dividing by 10^maxExponent may not be enough, here we allow division by up to 10^2*maxExponent which is by far enough
 	 // using ldd_power10() with negative exponents would save 1 or 2 divides here - but gives a few round the loop errors on the test program, so we have to leave the divides...
	 rexp= -rexp;
	 int exp=rexp;	 
	 if(rexp>4096) /* 10^4096 using ldd_power10() is "exact" (and fast) so this is a good one to duplicate */
	 	{
 		 rexp=4096;
 		 exp-=4096; // any excess which we will also need to divide by (if its > 0)
 		}
 	  else exp=0;  
	  long double tenh,tenl; 
 	  /* avoid divides unless exp very big (doing both with mults gives 3 errors round the loop in the test program ) */
 	  ldd_power10(&tenh,&tenl,-rexp); // 10^-rexp	     	  
 	  mult_ldd_ldd(rh,rl,*rh,*rl,tenh,tenl); // negative exponent means we divide by powers of 10 	  
 	  if(exp>0)
 	  	{  			   	   	  		
		 ldd_power10(&tenh,&tenl,exp); // 10^exp					  
 	  	 div_ldd_ldd(rh,rl,*rh,*rl,tenh,tenl); // negative exponent means we divide by powers of 10 - so divide by the rest of the exponent	  	 
 	    } 	     
	}	
 // else special case, rexp==0 - nothing to do
}
#else
/* we assume long double==double, and so leaverage the double-double function dd_mult_power10() */
// void dd_mult_power10( double *ohi, double *olo, double dh,double dl, int power )  // power can be +/-500. ASSUMES dh,dl came from a uint64 - this code will not work in general for all values dh,dl.
void ldd_power10(long double *rh, long double *rl , int n)// raise 10 to nth power (n can be positive or negative) - return long double double result , uses long double double maths internally to give effectively perfect results - even for denormalised numbers
{dd_mult_power10((double *)rh,(double *)rl,1.0,0.0,(int)n);
}

void ldd_mult_power10(long double *rh,long double *rl,long double rin_h, long double rin_l, int32_t rexp )  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 500
{
 dd_mult_power10((double *)rh,(double *)rl,(double)rin_h,(double)rin_l,(int)rexp);
}
#endif

#ifdef USE_ACCURATE_ALG
/* my implementation based on 2nd order formula for reciprocal in "Multiplicative Iteration for Reciprocals", Prof. W. Kahan, https://people.eecs.berkeley.edu/~wkahan/CS279/recip.pdf 
   for r=1/x a more accurate approximation R is given by:
   d=(1-x*r)
   R=d*r+r
   2 iterations done
*/
void div_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl)  // divides long double double a by b to give lomg double double "x"
{long double q1, q2, q3;
 long double rh,rl,r1h,r1l;  
 if(!isfinite(ah) || (ah==0.0L && al==0.0L))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!isfinite(bh))
 	{if(isinf(bh)) 
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	} 
 q1=ah/bh; // initial estimate
 mult_ld_ldd(&rh,&rl,q1,bh,bl); // r=q1*b
 sub_ldd_ldd(&rh,&rl,ah,al,rh,rl);//  r=a-q1*b
 q2=rh/bh;
 mult_ld_ldd(&r1h,&r1l,q2,bh,bl); // r1=q2*b
 sub_ldd_ldd(&rh,&rl,rh,rl,r1h,r1l);// r=r-(q2*b)
 q3=rh/bh;
 ld_fasttwosum(&q1,&q2,q1,q2);
 add_ldd_ldd(xh,xl,q1,q2,q3,0.0L);
}	
#else
void div_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl)  // divides long double double a by b to give lomg double double "x"
{long double hi;
 long double ch,cl,uh,ul;
 if(!isfinite(ah))
 	{*xh=ah/bh; // assume rules for this make sense
 	 *xl=0;
 	 return;
 	}
 else if(!isfinite(bh))
 	{if(isinf(bh)) 
 		{// special case: x/infinity=0
 		 *xh=0;
 		 *xl=0;
 		 return;
 		}
 	 else
	 	{*xh=ah/bh;// assume rules for this make sense
	 	 *xl=0;
	 	 return;
	 	}
	} 
 ch=ah/bh; // initial estimate
 uh=ch*bh;
 ul=fmal(ch,bh,-uh);
 cl=(((ah-uh)-ul)+al-ch*bl)/bh;
 hi=ch+cl;
 ld_fasttwosum(xh,xl,hi,cl+(ch-hi)); /* ensure result is "normalised" */
}
#endif

#if defined(__SIZEOF_INT128__) && !defined(__BORLANDC__) /* only allow if compiler supports __float128 & __int128 */
#if 0/* both options here give identical results with the current test program */
void U128toLDD(uint128_t u128,long double *xh,long double *xl) /* convert uint128 to long double double */
{ int128_t r; // needs to be signed, cannot overflow as x takes most of the bits in u128 before assignment to r
 long double x=u128;
 r=u128-(uint128_t)x;
 ld_twosum(xh,xl,x,r);// combine two halves of u128 to make the double double
}
#else
void U128toLDD(uint128_t u128,long double *xh,long double *xl) /* convert uint128 to long double double */
{uint128_t u128_l=u128 & UINT64_C(0xffffffffffffffff );// bottom 64 bits [ f80 has a 64 bit mantissa so 64 bits can fit exactly  ]
 uint128_t u128_h=u128^u128_l;// upper 64 bits
 if(u128_h==0)
 	{*xh=u128; // <= 64 bits so will fit into one f80
 	 *xl=0;
 	 return;
 	}
 ld_twosum(xh,xl,u128_h,u128_l);// combine two halves of u128 to make the double double
}
#endif
#endif

/* now restore gcc options to those set by the user */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) || defined(__clang__)
#pragma GCC pop_options
#endif


