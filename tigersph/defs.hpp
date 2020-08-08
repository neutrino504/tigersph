#pragma once



#define NDIM    3
#define NCHILD  2

#define ERROR() printf( "ERROR %s %i\n", __FILE__, __LINE__)


//#define USE_AVX512
#define USE_AVX2
//#define USE_AVX




#ifdef INDIRECT_DOUBLE
using ireal = double;
#define simd_real simd_double
#define SIMD_REAL_LEN SIMD_DOUBLE_LEN
#else
using ireal = float;
#define simd_real  simd_float
#define SIMD_REAL_LEN SIMD_FLOAT_LEN
#endif
