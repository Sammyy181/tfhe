#ifndef TOOM4_LAGRANGEHALFC_IMPL_H
#define TOOM4_LAGRANGEHALFC_IMPL_H

/**
 * Toom-4 "striding" transform processor (Bernstein; Bermudo Mera et al.), with Bodrato's
 * interpolation sequence, used in place of the FFT.
 *
 * The FFT processors represent a polynomial mod X^N+1 by its values at the roots of X^N+1.
 * Here the "transform domain" is the full Toom-4 evaluation tree instead: the polynomial is
 * split as a(x) = sum_{i<4} x^i A_i(x^4), the 4 pieces are evaluated at the 7 points
 * 0, 1/2, 1, -1, 2, -2, inf, and each of the 7 resulting polynomials of size N/4 is split again,
 * recursively, down to leaves of size <= TOOM4_LEAF. A polynomial in the transform domain is
 * thus 7^depth leaves of N/4^depth words, stored depth-first (child j of a node is the j-th
 * contiguous block of 1/7 of the node).
 *
 * - forward transform  (IntPolynomial_ifft, TorusPolynomial_ifft) = evaluation down the tree
 * - "pointwise" product (LagrangeHalfCPolynomialMul/AddMul/SubMul) = one negacyclic schoolbook
 *   product per leaf (the striding split keeps every sub-product negacyclic and size-preserving)
 * - inverse transform  (TorusPolynomial_fft) = Bodrato interpolation + recomposition up the tree
 *
 * The transform domain has two formats, tagged in each polynomial:
 * - OPERAND: the evaluation of the degree-3 split (point 1/2 scaled by 8, leading coefficient at inf);
 *   this is what the forward transform produces, and what the leaf products take as inputs;
 * - PRODUCT: the layout of a degree-6 product (point 1/2 scaled by 64, x^6 coefficient at inf);
 *   this is what the leaf products produce, and what the interpolation (inverse transform) takes.
 * The bootstrapping only ever multiplies OPERANDs (key, decomposed accumulator), accumulates
 * PRODUCTs, and interpolates PRODUCTs, so no conversion happens there. Any other combination
 * allowed by the FFT interface is still supported, through exact conversions:
 * OPERAND -> PRODUCT is a leaf product with the transform of 1, PRODUCT -> OPERAND is an inverse
 * then a forward transform.
 *
 * All the arithmetic is exact, in 64-bit words (mod 2^64). Each level of interpolation divides
 * by at most 2^4, so the result is exact mod 2^(64 - 4*depth), which must be >= 2^32 (depth <= 8).
 * Since every step is linear, products can be accumulated in the transform domain and
 * interpolated only once (this is what tGswFFTExternMulToTLwe does), and the bootstrapping key
 * is transformed once at key setup (tGswToFFTConvert), like with the FFT.
 */

#include <cstdint>
#include <cstddef>
#include "tfhe.h"
#include "polynomials.h"

#ifndef TOOM4_LEAF
#define TOOM4_LEAF 16   // recursion cutoff (best base of the striding Toom-4 for N=1024, see best_bases.csv)
#endif

class Toom4StridedProcessor {
public:
    const int32_t N;        // polynomial size
    const int32_t depth;    // number of Toom-4 levels
    const int32_t leaf;     // leaf size  N / 4^depth
    const int32_t nleaves;  // 7^depth
    const size_t size;      // words per transformed polynomial = nleaves * leaf

    explicit Toom4StridedProcessor(int32_t N);

    void forward(uint64_t* tree, const uint64_t* a) const;    // a: N words
    void backward(uint64_t* r, const uint64_t* tree) const;   // r: N words
    void mulAcc(uint64_t* acc, const uint64_t* a, const uint64_t* b, int32_t sign) const; // acc += sign * a.b

    void toProduct(uint64_t* dst, const uint64_t* src) const;   // OPERAND -> PRODUCT (dst != src)
    void toOperand(uint64_t* dst, const uint64_t* src) const;   // PRODUCT -> OPERAND (dst may be src)

    static const Toom4StridedProcessor* get(int32_t N);       // one cached processor per N

    /** thread-local scratch: [0,N) coefficient buffer, [N,4N) recursion, then 3 transformed polynomials */
    uint64_t* workspace() const;
    uint64_t* tmpTree(int32_t i) const { return workspace() + 4 * (size_t) N + (size_t) i * size; }
private:
    size_t wsWords;
    uint64_t* one;    // OPERAND transform of the constant polynomial 1
};

enum Toom4Format { TOOM4_ZERO = 0, TOOM4_OPERAND = 1, TOOM4_PRODUCT = 2 };

/**
 * Must be layout-compatible with LagrangeHalfCPolynomial (at most 2 pointers)
 */
struct LagrangeHalfCPolynomial_IMPL {
    uint64_t* coefs;                     // the evaluation tree, proc->size words, then the format tag
    const Toom4StridedProcessor* proc;

    uint64_t& format() { return coefs[proc->size]; }
    uint64_t format() const { return coefs[proc->size]; }

    explicit LagrangeHalfCPolynomial_IMPL(int32_t N);
    ~LagrangeHalfCPolynomial_IMPL();
};

#endif // TOOM4_LAGRANGEHALFC_IMPL_H
