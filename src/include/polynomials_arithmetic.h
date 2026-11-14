#ifndef POLYNOMIALS_ARITHMETIC_H
#define POLYNOMIALS_ARITHMETIC_H

#include "tfhe_core.h"
#include "polynomials.h"

///@file
///@brief This file contains the operations on polynomials


/**  TorusPolynomial = 0 */
EXPORT void torusPolynomialClear(TorusPolynomial* result);

/**  TorusPolynomial = random */
EXPORT void torusPolynomialUniform(TorusPolynomial* result);

/**  TorusPolynomial = TorusPolynomial */
EXPORT void torusPolynomialCopy(TorusPolynomial* result, const TorusPolynomial*  sample);

/**  TorusPolynomial + TorusPolynomial */
EXPORT void torusPolynomialAdd(TorusPolynomial* result, const TorusPolynomial* poly1, const TorusPolynomial* poly2);

/**  TorusPolynomial += TorusPolynomial */
EXPORT void torusPolynomialAddTo(TorusPolynomial* result, const TorusPolynomial* poly2);


/**  TorusPolynomial - TorusPolynomial */
EXPORT void torusPolynomialSub(TorusPolynomial* result, const TorusPolynomial* poly1, const TorusPolynomial* poly2);

/**  TorusPolynomial -= TorusPolynomial */
EXPORT void torusPolynomialSubTo(TorusPolynomial* result, const TorusPolynomial* poly2);

/**  TorusPolynomial + p*TorusPolynomial */
EXPORT void torusPolynomialAddMulZ(TorusPolynomial* result, const TorusPolynomial* poly1, int32_t p, const TorusPolynomial* poly2);

/**  TorusPolynomial += p*TorusPolynomial */
EXPORT void torusPolynomialAddMulZTo(TorusPolynomial* result, const int32_t p, const TorusPolynomial* poly2);

/**  TorusPolynomial - p*TorusPolynomial */
EXPORT void torusPolynomialSubMulZ(TorusPolynomial* result, const TorusPolynomial* poly1, const int32_t p, const TorusPolynomial* poly2);

/**  TorusPolynomial -= p*TorusPolynomial */
EXPORT void torusPolynomialSubMulZTo(TorusPolynomial* result, int32_t p, const TorusPolynomial* poly2);

/**  TorusPolynomial = (X^a-1)* TorusPolynomial */
EXPORT void torusPolynomialMulByXaiMinusOne(TorusPolynomial* result, int32_t a, const TorusPolynomial* source);

/** result= X^{a}*source */
EXPORT void torusPolynomialMulByXai(TorusPolynomial* result, int32_t a, const TorusPolynomial* source);

/**  Norme Euclidienne d'un IntPolynomial */
EXPORT double intPolynomialNormSq2(const IntPolynomial* poly);

/**  IntPolynomial = 0 */
EXPORT void intPolynomialClear(IntPolynomial* result);

/**  result = source */
EXPORT void intPolynomialCopy(IntPolynomial* result, const IntPolynomial* source);

/**  result += poly1 */
EXPORT void intPolynomialAddTo(IntPolynomial* result, const IntPolynomial* poly1);

/**  result = (X^ai-1) * source */
EXPORT void intPolynomialMulByXaiMinusOne(IntPolynomial* result, int32_t ai, const IntPolynomial* source);

/**  ILA: Norme infini de la distance entre deux TorusPolynomial */
EXPORT double torusPolynomialNormInftyDist(const TorusPolynomial* poly1, const TorusPolynomial* poly2);

// Norme 2 d'un IntPolynomial
EXPORT double intPolynomialNorm2sq(const IntPolynomial* poly);

// Norme infini de la distance entre deux IntPolynomial
EXPORT double intPolynomialNormInftyDist(const IntPolynomial* poly1, const IntPolynomial* poly2);




/**
 * This is the naive external multiplication of an integer polynomial
 * with a torus polynomial. (this function should yield exactly the same
 * result as the karatsuba or fft version, but should be slower) 
 */
EXPORT void torusPolynomialMultNaive(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);


/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly)
 * by using Karatsuba
 * WARNING: N must be a power of 2 to use this function. Else, the
 * behaviour is unpredictable
 */
EXPORT void torusPolynomialMultKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * result += poly1 * poly2 (via karatsuba)
 * WARNING: N must be a power of 2 to use this function. Else, the
 * behaviour is unpredictable
 */
EXPORT void torusPolynomialAddMulRKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * result -= poly1 * poly2 (via karatsuba)
 * WARNING: N must be a power of 2 to use this function. Else, the
 * behaviour is unpredictable
 */
EXPORT void torusPolynomialSubMulRKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly)
 * by using Toom-Cook 3
 */
EXPORT void torusPolynomialMultToomCook3(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly)
 * by using Toom-Cook 4
 */
EXPORT void torusPolynomialMultToomCook4(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly)
 * by using Toom-Cook 5
 */
EXPORT void torusPolynomialMultToomCook5(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly) mod X^N+1
 * by using striding Toom-Cook 4 (interleaved split, negacyclic sub-products) with
 * Bodrato's interpolation sequence
 */
EXPORT void torusPolynomialMultToomCook4Strided(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * result = sum_{j<tau} poly1[j] * poly2[j] (mod X^N+1), by using Toom-Cook 3 (resp. 4).
 * If lazy is true, the interpolation is done only once, on the sum of the products in the
 * evaluation space (instead of once per product). Both modes give the same result.
 */
EXPORT void torusPolynomialMultSumToomCook3(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy);
EXPORT void torusPolynomialMultSumToomCook4(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy);

/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly) mod X^N+1
 * by using the Number Theoretic Transform (3 primes + CRT). Unlike the FFT, it is exact.
 * N must be a power of 2 (at most 2^22).
 */
EXPORT void torusPolynomialMultNTT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/**
 * result = sum_{j<tau} poly1[j] * poly2[j] (mod X^N+1), by using the NTT.
 * If lazy is true, the interpolation (inverse NTT + CRT) is done only once, on the sum of the
 * products in the evaluation space. Both modes give the same exact result, provided that
 * tau * N * max|poly1| < 2^58.
 */
EXPORT void torusPolynomialMultSumNTT(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy);

/**
 * Toom-Cook interpolation variants (toom-variants.cpp), used to compare interpolation strategies
 * under identical evaluation, workspace and recomposition code. "base" is the recursion cutoff
 * (naive product for sizes <= base). Returns 0 on success, -1 if the recursion depth would exceed
 * the 32-bit precision budget (nothing is computed), -2 for an unknown variant.
 */
enum ToomVariant {
    TV_T3_OWN = 0,      // classical Toom-3, own interpolation sequence
    TV_T3_BODRATO = 1,  // classical Toom-3, Bodrato's interpolation sequence
    TV_T3_GENERIC = 2,  // classical Toom-3, generic (inverse Vandermonde) interpolation
    TV_T4_BODRATO = 3,  // classical Toom-4, Bodrato's interpolation sequence
    TV_T4_GENERIC = 4,  // classical Toom-4, generic interpolation
    TV_T4_STRIDED = 5,  // striding Toom-4 (negacyclic sub-products), Bodrato's sequence
    TV_T5_GENERIC = 6   // classical Toom-5, generic interpolation
};
EXPORT int32_t torusPolynomialMultToomVariant(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2, const int32_t variant, const int32_t base);
/** number of recursion levels used by a variant for size N and cutoff base */
EXPORT int32_t toomVariantDepth(const int32_t variant, const int32_t N, const int32_t base);

//#define torusPolynomialMulR torusPolynomialMultKaratsuba
//#define torusPolynomialAddMulR torusPolynomialAddMulRKaratsuba
//#define torusPolynomialSubMulR torusPolynomialSubMulRKaratsuba

#define torusPolynomialMulR torusPolynomialMultFFT
#define torusPolynomialAddMulR torusPolynomialAddMulRFFT
#define torusPolynomialSubMulR torusPolynomialSubMulRFFT


#endif //POLYNOMIALS_ARITHMETIC_H
