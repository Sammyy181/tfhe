/*
 * Two-level Nussbaumer(16, 8) transform (see nussbaumer.h for the overview).
 *
 * Both levels are a DFT whose "numbers" are polynomials mod y^n+1 and whose
 * root of unity is a power of y, so that a multiplication by a root is just
 * a negacyclic rotation of the n coefficients (no real multiplication).
 *
 * The forward DFT is decimation-in-frequency and leaves its output in
 * bit-reversed order; the inverse DFT is decimation-in-time and takes its
 * input in that same order. The pointwise products in between do not care
 * about the order, so no bit-reversal permutation is ever needed.
 */

#include <cstring>
#include <cassert>
#include "nussbaumer.h"

using namespace std;

// dst = src * y^k in Z[y]/(y^n+1), for 0 <= k < 2n. dst and src must not overlap.
static inline void negacyclicRotate(NussWord *dst, const NussWord *src, const int32_t n, int32_t k) {
    if (k < n) {
        for (int32_t i = 0; i < n - k; i++) dst[i + k] = src[i];
        for (int32_t i = n - k; i < n; i++) dst[i + k - n] = -src[i];
    } else {
        k -= n;
        for (int32_t i = 0; i < n - k; i++) dst[i + k] = -src[i];
        for (int32_t i = n - k; i < n; i++) dst[i + k - n] = src[i];
    }
}

// In-place forward DFT of L elements of Z[y]/(y^n+1) (L <= 2n, both powers of 2),
// with the L-th root of unity y^(2n/L). Output is in bit-reversed order.
static void ringForwardDFT(NussWord *x, const int32_t L, const int32_t n) {
    NussWord diff[NUSS_L1_SIZE];
    assert(n <= NUSS_L1_SIZE);
    for (int32_t h = L / 2; h >= 1; h /= 2) {
        for (int32_t s = 0; s < L; s += 2 * h) {
            for (int32_t j = 0; j < h; j++) {
                NussWord *u = x + (s + j) * n;
                NussWord *v = x + (s + j + h) * n;
                // (u, v) -> (u + v, (u - v) * w^j) where w = y^(n/h) has order 2h
                for (int32_t t = 0; t < n; t++) {
                    diff[t] = u[t] - v[t];
                    u[t] += v[t];
                }
                negacyclicRotate(v, diff, n, j * (n / h));
            }
        }
    }
}

// In-place inverse of ringForwardDFT, up to a factor L: takes the bit-reversed
// output of ringForwardDFT and returns L times the original input.
static void ringInverseDFT(NussWord *x, const int32_t L, const int32_t n) {
    NussWord rot[NUSS_L1_SIZE];
    assert(n <= NUSS_L1_SIZE);
    for (int32_t h = 1; h < L; h *= 2) {
        for (int32_t s = 0; s < L; s += 2 * h) {
            for (int32_t j = 0; j < h; j++) {
                NussWord *u = x + (s + j) * n;
                NussWord *v = x + (s + j + h) * n;
                // (u, v) -> (u + v * w^-j, u - v * w^-j)
                const int32_t k = j * (n / h);
                negacyclicRotate(rot, v, n, k == 0 ? 0 : 2 * n - k);
                for (int32_t t = 0; t < n; t++) {
                    v[t] = u[t] - rot[t];
                    u[t] += rot[t];
                }
            }
        }
    }
}


EXPORT void nussbaumerForward(NussWord *result, const int32_t *poly) {
    NussWord lvl1[NUSS_L1_POINTS * NUSS_L1_SIZE];

    // level 1: poly(X) = sum_{i<16} A_i(Y) X^i with Y = X^16.
    // A_i has 64 coefs; slots 16..31 of the length-32 transform are zero.
    for (int32_t i = 0; i < NUSS_L1_PIECES; i++)
        for (int32_t t = 0; t < NUSS_L1_SIZE; t++)
            lvl1[i * NUSS_L1_SIZE + t] = (NussWord) (int64_t) poly[NUSS_L1_PIECES * t + i];
    memset(lvl1 + NUSS_L1_PIECES * NUSS_L1_SIZE, 0,
           (NUSS_L1_POINTS - NUSS_L1_PIECES) * NUSS_L1_SIZE * sizeof(NussWord));
    ringForwardDFT(lvl1, NUSS_L1_POINTS, NUSS_L1_SIZE);

    // level 2, on each of the 32 polynomials E(Y) mod Y^64+1:
    // E(Y) = sum_{j<8} e_j(Z) Y^j with Z = Y^8.
    // e_j has 8 coefs; slots 8..15 of the length-16 transform are zero.
    for (int32_t s = 0; s < NUSS_L1_POINTS; s++) {
        const NussWord *e = lvl1 + s * NUSS_L1_SIZE;
        NussWord *out = result + s * NUSS_L2_POINTS * NUSS_L2_SIZE;
        for (int32_t j = 0; j < NUSS_L2_PIECES; j++)
            for (int32_t t = 0; t < NUSS_L2_SIZE; t++)
                out[j * NUSS_L2_SIZE + t] = e[NUSS_L2_PIECES * t + j];
        memset(out + NUSS_L2_PIECES * NUSS_L2_SIZE, 0,
               (NUSS_L2_POINTS - NUSS_L2_PIECES) * NUSS_L2_SIZE * sizeof(NussWord));
        ringForwardDFT(out, NUSS_L2_POINTS, NUSS_L2_SIZE);
    }
}


EXPORT void nussbaumerInverse(Torus32 *result, const NussWord *source) {
    NussWord lvl1[NUSS_L1_POINTS * NUSS_L1_SIZE];
    NussWord lvl2[NUSS_L2_POINTS * NUSS_L2_SIZE];
    NussWord out[NUSS_N];

    // undo level 2: 16 polynomials of size 8 -> one polynomial mod Y^64+1
    memset(lvl1, 0, sizeof(lvl1));
    for (int32_t s = 0; s < NUSS_L1_POINTS; s++) {
        memcpy(lvl2, source + s * NUSS_L2_POINTS * NUSS_L2_SIZE, sizeof(lvl2));
        ringInverseDFT(lvl2, NUSS_L2_POINTS, NUSS_L2_SIZE);
        // E(Y) = sum_j c_j(Y^8) Y^j, reduced with Y^64 = -1
        NussWord *e = lvl1 + s * NUSS_L1_SIZE;
        for (int32_t j = 0; j < NUSS_L2_POINTS; j++) {
            for (int32_t t = 0; t < NUSS_L2_SIZE; t++) {
                const int32_t idx = NUSS_L2_PIECES * t + j;
                if (idx < NUSS_L1_SIZE) e[idx] += lvl2[j * NUSS_L2_SIZE + t];
                else e[idx - NUSS_L1_SIZE] -= lvl2[j * NUSS_L2_SIZE + t];
            }
        }
    }

    // undo level 1: 32 polynomials of size 64 -> one polynomial mod X^1024+1
    ringInverseDFT(lvl1, NUSS_L1_POINTS, NUSS_L1_SIZE);
    memset(out, 0, sizeof(out));
    // result(X) = sum_i C_i(X^16) X^i, reduced with X^1024 = -1
    for (int32_t i = 0; i < NUSS_L1_POINTS; i++) {
        for (int32_t t = 0; t < NUSS_L1_SIZE; t++) {
            const int32_t idx = NUSS_L1_PIECES * t + i;
            if (idx < NUSS_N) out[idx] += lvl1[i * NUSS_L1_SIZE + t];
            else out[idx - NUSS_N] -= lvl1[i * NUSS_L1_SIZE + t];
        }
    }

    // the two inverse DFTs left a factor 32 * 16 = 2^9: this is an exact division
    for (int32_t i = 0; i < NUSS_N; i++)
        result[i] = (Torus32) (uint32_t) (out[i] >> NUSS_LOG_SCALE);
}
