#ifndef NUSSBAUMER_H
#define NUSSBAUMER_H

///@file
///@brief Two-level Nussbaumer(16, 8) transform + Winograd-Schoolbook
///       multiplication, as a replacement for the FFT in the external product.
///
/// A polynomial mod X^1024+1 is mapped to 512 polynomials of size 8:
///
///   1 product of size 1024                              (mod X^1024+1)
///      | level 1: 16 pieces of 64 coefs, Y = X^16, length-32 transform
///      |          with root Y^4 in Z[Y]/(Y^64+1)
///   32 products of size 64                              (mod Y^64+1)
///      | level 2: 8 pieces of 8 coefs, Z = Y^8, length-16 transform
///      |          with root Z in Z[Z]/(Z^8+1)
///   512 products of size 8                              (mod Z^8+1)
///
/// The 512 small products are computed with the Winograd-Schoolbook F(2,3)
/// algorithm. Products can be accumulated in the transform domain, and the
/// inverse transform applied once at the end (lazy interpolation).
///
/// The inverse transform divides by 32*16 = 2^9, which is not invertible
/// mod 2^32. Transform-domain values are therefore 64-bit words (arithmetic
/// mod 2^64), and the division is one exact shift at the very end of
/// nussbaumerInverse. All results are bit-exact mod 2^32.

#include "tfhe_core.h"
#include "polynomials.h"
#include "tlwe.h"
#include "tgsw.h"
#include "lwekeyswitch.h"
#include "lwebootstrappingkey.h"

/** one coefficient in the Nussbaumer domain */
typedef uint64_t NussWord;

#define NUSS_N 1024         ///< size of the input polynomials (mod X^N+1)
#define NUSS_L1_PIECES 16   ///< level 1: number of pieces
#define NUSS_L1_POINTS 32   ///< level 1: transform length (= 2 * pieces)
#define NUSS_L1_SIZE 64     ///< level 1: size of each piece / product
#define NUSS_L2_PIECES 8    ///< level 2: number of pieces
#define NUSS_L2_POINTS 16   ///< level 2: transform length (= 2 * pieces)
#define NUSS_L2_SIZE 8      ///< level 2: size of each piece / product
#define NUSS_NB_PRODUCTS (NUSS_L1_POINTS * NUSS_L2_POINTS)  ///< 512 small products
#define NUSS_SIZE (NUSS_NB_PRODUCTS * NUSS_L2_SIZE)         ///< 4096 words per polynomial
#define NUSS_LOG_SCALE 9    ///< the inverse transform divides by 2^9 = 32 * 16

// ---------------------------------------------------------------------------
// nussbaumer.cpp: the transform
// ---------------------------------------------------------------------------

/**
 * Forward Nussbaumer(16, 8) transform.
 * @param result NUSS_SIZE words (the 512 size-8 polynomials, one after the other)
 * @param poly NUSS_N coefficients (the coefs of an IntPolynomial or of a TorusPolynomial)
 */
EXPORT void nussbaumerForward(NussWord *result, const int32_t *poly);

/**
 * Inverse Nussbaumer(16, 8) transform, including the reduction mod X^N+1
 * and the final division by 2^NUSS_LOG_SCALE.
 * @param result NUSS_N torus coefficients
 * @param source NUSS_SIZE words (a product, or a sum of products)
 */
EXPORT void nussbaumerInverse(Torus32 *result, const NussWord *source);

// ---------------------------------------------------------------------------
// multiplication.cpp: multiplication in the Nussbaumer domain
// ---------------------------------------------------------------------------

/**
 * result = a * b in the Nussbaumer domain: 512 products of size 8 mod Z^8+1,
 * each computed with Winograd-Schoolbook F(2,3).
 * a is the "tile" operand and b the "filter" operand of the paper.
 * result must not overlap a or b. All three have NUSS_SIZE words.
 */
EXPORT void nussbaumerPointwiseMulWS(NussWord *result, const NussWord *a, const NussWord *b);

/**
 * accum += poly1 * poly2, staying in the Nussbaumer domain (no inverse transform).
 * poly1 is transformed here; poly2 is already in the Nussbaumer domain.
 * @param accum ncols accumulators of NUSS_SIZE words, one after the other
 * @param poly1 polynomial in coefficient form (N = NUSS_N)
 * @param poly2Nuss ncols transformed polynomials of NUSS_SIZE words, one after the other
 * @param ncols number of polynomials poly1 is multiplied with (poly1 is
 *              transformed only once); use 1 for a single product
 */
EXPORT void torusPolynomialAddMulRNussbaumer(NussWord *accum, const IntPolynomial *poly1,
                                             const NussWord *poly2Nuss, int32_t ncols);

/**
 * Standalone multiplication: result = poly1 * poly2 mod X^N+1
 * (forward transform of both, pointwise product, inverse transform).
 */
EXPORT void torusPolynomialMulNussWS(TorusPolynomial *result, const IntPolynomial *poly1,
                                     const TorusPolynomial *poly2);

// ---------------------------------------------------------------------------
// tgsw-fft-operations.cpp: TGSW samples in the Nussbaumer domain
// (counterpart of TGswSampleFFT)
// ---------------------------------------------------------------------------

struct TGswSampleNussbaumer {
    /// (k+1)*l TLWE rows of k+1 polynomials, each of NUSS_SIZE words:
    /// polynomial j of row p starts at data + (p*(k+1) + j) * NUSS_SIZE
    NussWord *data;
    int32_t k;
    int32_t l;
};
typedef struct TGswSampleNussbaumer TGswSampleNussbaumer;

EXPORT TGswSampleNussbaumer *new_TGswSampleNussbaumer(const TGswParams *params);
EXPORT TGswSampleNussbaumer *new_TGswSampleNussbaumer_array(int32_t nbelts, const TGswParams *params);
EXPORT void delete_TGswSampleNussbaumer(TGswSampleNussbaumer *obj);
EXPORT void delete_TGswSampleNussbaumer_array(int32_t nbelts, TGswSampleNussbaumer *obj);

/** transforms every polynomial of a TGSW sample (counterpart of tGswToFFTConvert) */
EXPORT void tGswToNussbaumerConvert(TGswSampleNussbaumer *result, const TGswSample *source,
                                    const TGswParams *params);

/**
 * External product: accum = gsw (*) accum (counterpart of tGswFFTExternMulToTLwe).
 * All the kpl products are accumulated in the Nussbaumer domain; the inverse
 * transform is applied only once per output polynomial.
 */
EXPORT void tGswNussbaumerExternMulToTLwe(TLweSample *accum, const TGswSampleNussbaumer *gsw,
                                          const TGswParams *params);

// ---------------------------------------------------------------------------
// lwe-bootstrapping-functions-nussbaumer.cpp: bootstrapping
// (counterpart of LweBootstrappingKeyFFT and the tfhe_*_FFT functions)
// ---------------------------------------------------------------------------

struct LweBootstrappingKeyNussbaumer {
    const LweParams *in_out_params;
    const TGswParams *bk_params;
    const TLweParams *accum_params;
    const LweParams *extract_params;
    const TGswSampleNussbaumer *bkNussbaumer; ///< the bootstrapping key, in the Nussbaumer domain
    const LweKeySwitchKey *ks; ///< NOT a copy: points to the keyswitch key of the source LweBootstrappingKey
};
// (the LweBootstrappingKeyNussbaumer typedef is in tfhe_core.h)

/** bk must outlive the returned key (the keyswitch key is shared, not copied) */
EXPORT LweBootstrappingKeyNussbaumer *new_LweBootstrappingKeyNussbaumer(const LweBootstrappingKey *bk);
EXPORT void delete_LweBootstrappingKeyNussbaumer(LweBootstrappingKeyNussbaumer *obj);

EXPORT void tfhe_blindRotate_Nussbaumer(TLweSample *accum, const TGswSampleNussbaumer *bk,
                                        const int32_t *bara, const int32_t n, const TGswParams *bk_params);
EXPORT void tfhe_blindRotateAndExtract_Nussbaumer(LweSample *result, const TorusPolynomial *v,
                                                  const TGswSampleNussbaumer *bk, const int32_t barb,
                                                  const int32_t *bara, const int32_t n,
                                                  const TGswParams *bk_params);
EXPORT void tfhe_bootstrap_woKS_Nussbaumer(LweSample *result, const LweBootstrappingKeyNussbaumer *bk,
                                           Torus32 mu, const LweSample *x);
EXPORT void tfhe_bootstrap_Nussbaumer(LweSample *result, const LweBootstrappingKeyNussbaumer *bk,
                                      Torus32 mu, const LweSample *x);

#endif //NUSSBAUMER_H
