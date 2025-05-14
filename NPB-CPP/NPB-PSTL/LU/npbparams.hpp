/* CLASS = A */
					/*
					  c  This file is generated automatically by the setparams utility.
					  c  It sets the number of processors and the class_npb of the NPB
					  c  in this directory. Do not modify it by hand.
					 */

/* full problem size */
#define	ISIZ1	64
#define	ISIZ2	64
#define	ISIZ3	64
/* number of iterations and how often to print the norm */
#define	ITMAX_DEFAULT	250
#define	INORM_DEFAULT	250
#define	DT_DEFAULT	2.0
#define	CONVERTDOUBLE	FALSE
#define COMPILETIME "04 May 2025"
#define NPBVERSION "4.1"
#define COMPILERVERSION "13.3.0"
#define CS1 "g++ -std=c++20"
#define CS2 "$(CC)"
#define CS3 "-lm -L$(OTBB_ROOT)/lib/ -ltbb "
#define CS4 "-I../common -I$(OTBB_ROOT)/include/ "
#define CS5 "-O3"
#define CS6 "-O3"
#define CS7 "randdp"
