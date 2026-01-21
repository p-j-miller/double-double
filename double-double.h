/* double-double.h
   ================
   
   double-double routines for basic doubles, long doubles (e.g. __float80) and __float128 types.
   
   
   This file Peter Miller 25/5/2020 but parts go back a long way beyond that.
   
   1v0 - 1st version of this file.
   Changes to 2v0 1/2026
   - bracketed by __DOUBLE_DOUBLE_H. Added ddtoU64()
   - issue in f128_mult_dd_dd() where it would return inf prematurely fixed
   - changed so long double function use "raw" long doubles rather than forcing them to be f80 (which is a type that's not defined in the C99 standard).
   - in some cases long doubles types are identical to double (i.e. 8 bytes), in others they are 10, 12 or 16 bytes (perhaps with only 12 bytes actually used) - all of this is allowed by the C99 standard.
   - the code here should always work - but the resolution will obviously depend upon the underlying long double datatype.
	- This is not actually as change in functionality as previously there was a typedef:
			typedef long double f80_t;
		but this was misleading as it suggested the routines always provided an 80 bit datatype	
   - u2_64toDD_f128() function added (works with u2_64.[ch] which gives a pseudo uint128 using uint64's)
   - power10() functions added	

 Note this code leverages (all available at https://github.com/p-j-miller )
	128 bit functionality using two uint64_t's from u2_64 
	power of 10 tables from power10
	my_printf for debugging from my_printf   
	  	  	  
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
 
#ifndef _double_double_h
#define _double_double_h
 #include "../u2_64-128bits-with-two-u64/u2_64.h" // type for 128 bit unsigned without needing compiler support
 #ifdef __cplusplus
 extern "C" {
 #endif 
#if defined(__SIZEOF_FLOAT128__) /* only allow if compiler supports __float128 & __int128 */
 typedef __float128 f128_t;
#endif
#if defined( __SIZEOF_INT128__)
  typedef __uint128_t uint128_t; // same format as stdint.h
  typedef __int128_t int128_t;
#endif

	// double double functions.
	void add_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl);  // adds double double x and y to give double double "z"
	void sub_dd_dd(double *zh, double *zl,double xh, double xl,double yh, double yl);  // subtracts double double y from x to give double double "z"
	void sub_dd_d(double *zh, double *zl,double xh, double xl,double y);  // subtracts double y from double double x to give double double "z"
	void mult_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl);  // multiplies double double a and b to give double double "x"
	void mult_d_dd( double *xh, double *xl,double a,double bh, double bl);  // multiplies a and double double  b to give double double "x"
	void dd_power(double *rh, double *rl,double x, unsigned int n);// raise x to nth power - return double double result , uses double double maths internally to minimise the error
	void dd_mult_power10( double *ohi, double *olo, double dh,double dl, int power );  // power can be +/-500. ASSUMES dh,dl came from a uint64 - this code will not work in general for all values dh,dl.
	void div_dd_dd( double *xh, double *xl,double ah, double al,double bh, double bl) ; // divides double double a by b to give double double "x"
	void U64toDD(uint64_t u64,double *xh,double *xl); /* convert uint64 to double double */
	uint64_t ddtoU64(double xh,double xl);  // convert double double to uint64, deal with case where x may be a very small negative number

	#if defined(__SIZEOF_FLOAT128__) /* only allow if compiler supports __float128  */
	// double-double based on two __float128 
	void f128_add_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl);  // adds double double x and y to give double double "z"
	void f128_sub_dd_dd(f128_t *zh, f128_t *zl,f128_t xh, f128_t xl,f128_t yh, f128_t yl);  // subtracts double double y from x to give double double "z"
	void f128_mult_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl);  // multiplies double double a and b to give double double "x"
	void f128_mult_d_dd( f128_t *xh, f128_t *xl,f128_t a,f128_t bh, f128_t bl) ; // multiplies a and double double  b to give double double "x"
	void f128_dd_power(f128_t *rh, f128_t *rl,f128_t x, unsigned int n); // raise x to nth power - return double double result , uses double double maths internally to minimise the error
	void f128_power10(f128_t *rh, f128_t *rl , int n);// raise 10 to nth power (n can be positive or negative) - return f128_t double result , uses long double f128_t maths internally to give effectively perfect results - even for denormalised numbers
	void f128_mult_power10(__float128 *rh,__float128 *rl, __float128 rin_h, __float128 rin_l, int32_t rexp );  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
	void f128_div_dd_dd( f128_t *xh, f128_t *xl,f128_t ah, f128_t al,f128_t bh, f128_t bl);  // divides double double a by b to give double double "x"
	void u2_64toDD_f128(u2_64 u,f128_t *xh,f128_t *xl); /* convert u2_64 (128 bits) to double double f128 */
	 #if defined( __SIZEOF_INT128__)  /* only allow if compiler also supports  __int128 */
	  void U128toDD_f128(uint128_t u128,f128_t *xh,f128_t *xl); /* convert uint128 to double double f128 */
	 #endif
	#endif
	// double double based on two long doubles 
	void ld_twosum(long double *xh, long double *xl,long double a,long double b);  // adds a and b to give double double "x". 
	void add_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl);  // adds long double double x and y to give long double double "z"
	void sub_ldd_ldd(long double *zh, long double *zl,long double xh, long double xl,long double yh, long double yl);  // subtracts long double double y from x to give long double double "z"
	void mult_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl);  // multiplies long double double a and b to give long double double "x"
	void mult_ld_ldd( long double *xh, long double *xl,long double a,long double bh, long double bl);  // multiplies long double a and long double double  b to give long double double "x"
	void ldd_power(long double *rh, long double *rl,long double x, unsigned int n); // raise x to nth power - return long double double result , uses long double double maths internally to minimise the error
	void ldd_power10(long double *rh, long double *rl , int n); // raise 10 to nth power (n can be positive or negative) - return long double double result , uses long double double maths internally to minimise the error
	void ldd_mult_power10(long double *rh,long double *rl,long double rin_h, long double rin_l, int32_t rexp );  // rh/rl=rin_h/l*10^rexp  rexp can be +/- 
	void div_ldd_ldd( long double *xh, long double *xl,long double ah, long double al,long double bh, long double bl);  // divides long double double a by b to give long double double "x"
	#ifdef __SIZEOF_INT128__ /* only allow if compiler supports __float128 & __int128 */
	 void U128toLDD(uint128_t u128,long double *xh,long double *xl); /* convert uint128 to long double double */
	#endif

 #ifdef __cplusplus
    }
 #endif

#endif // ifndef _double_double_h
