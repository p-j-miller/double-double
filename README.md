# double-double
double-double routines for basic doubles, long doubles (e.g. __float80) and __float128 types

 The double-double representation represents numbers by a pair of native floating point types. This doubles the resolution (in mantissa bits or decimal digits),
   for example meaning that double-doubles based on the standard IEEE-754 64 bit double have 30 decimal digits (106 bits)
   compared to a single "native" double with 15 decimal digits (53bits). The exponent range is not changed.
   This approach was originally proposed by 
   D. E. Knuth in "The Art of Computer Programming, Volume II: Seminumerical Algorithms", 1st Edition, 1969, e.g. pp 237 and enhanced by amongst others 
   T. Dekker, "A floating-point technique for extending the available precision", Numerische Mathematik, 18(3) in 1971.
   See also IEEE trans. computers Vol 58 nos 7 July 2009 pp 994 "Accurate floating point product and exponentiation".
   Section 5.1 of the paper "Floating-Point Arithmetic" by Boldo, Jeannerod , Melquiond and Mukker (2023) which is at https://guillaume.melquiond.fr/doc/23-actanum.pdf ,
   gives an up to date status, and detailed error bounds for double word arithmetic of the type used here.

# Functions Provided
This library provides +,-,*,/ and power functions as well as some conversion functions.
~~~
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
~~~
# Installation
 It is recommended that these files (double-double.[ch] and main.c) are placed in a directory called double-double

Note this code leverages (all available at https://github.com/p-j-miller )
~~~
	128 bit functionality using two uint64_t's from u2_64 
	power of 10 tables from power10
	my_printf for debugging from my_printf
~~~
main.c is a test program which can be compiled under Windows using winlibs gcc 15.2.0 (the changes for Linux should be obvious) with the command:
~~~
C:\winlibs\winlibs-x86_64-posix-seh-gcc-15.2.0-mingw-w64ucrt-13.0.0-r2\mingw64\bin\gcc -Wall -m64 -fexcess-precision=standard -Ofast  -std=gnu99 -I. main.c ../double-double/double-double.c ../u2_64-128bits-with-two-u64/u2_64.c ../my_printf/my_printf.c  ../fma/fmaq.c -lquadmath -static -o test.exe
~~~
The full expected output from the test program is shown in main.c, but it's output should finish with the line "All tests passed".

Using fmaq.c is optional (the tests should be passed without it), but using it seems to improve the execution speed of the __float128 functions.

