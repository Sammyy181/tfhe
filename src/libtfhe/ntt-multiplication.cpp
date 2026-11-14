#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <map>
#include <mutex>
#include "polynomials_arithmetic.h"

using namespace std;

/**
 * Exact multiplication of an IntPolynomial by a TorusPolynomial mod X^N+1 with the
 * Number Theoretic Transform (NTT).
 *
 * The product has to be computed over the integers (the torus is Z/2^32, which has no
 * suitable roots of unity), so we compute it modulo 3 NTT-friendly primes q_0,q_1,q_2
 * (each q_i = 1 mod 2N, so that Z/q_i contains a primitive 2N-th root of unity psi, which
 * makes the transform negacyclic: no twisting or reduction mod X^N+1 is needed). The three
 * residues are then recombined with the CRT (Garner), centered around 0, and reduced mod 2^32.
 * Unlike the FFT, the result is exact as long as the true integer coefficients stay in
 * (-P/2, P/2), with P = q_0*q_1*q_2 ~ 2^90. The coefficients of poly2 are taken in
 * [-2^31, 2^31), so this holds as soon as tau * N * max|poly1| < 2^58.
 */

typedef uint64_t u64;
typedef unsigned __int128 u128;

static const int32_t NTT_NB_PRIMES = 3;
static const u64 NTT_Q[NTT_NB_PRIMES] = {2013265921ULL, 998244353ULL, 754974721ULL}; // 15*2^27+1, 119*2^23+1, 45*2^24+1
static const int32_t NTT_MAX_LOGN = 22;   // 2N must divide q_i-1 for the smallest 2-adic valuation (2^23)

static u64 mulmod(u64 a, u64 b, u64 q) { return (a * b) % q; }   // a,b < q < 2^31

static u64 powmod(u64 a, u64 e, u64 q) {
    u64 r = 1;
    a %= q;
    while (e) {
        if (e & 1) r = mulmod(r, a, q);
        a = mulmod(a, a, q);
        e >>= 1;
    }
    return r;
}

static int32_t bitrev(int32_t x, int32_t logn) {
    int32_t r = 0;
    for (int32_t i = 0; i < logn; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

struct NttTables {
    int32_t N;
    u64* fwd[NTT_NB_PRIMES];            // psi^bitrev(i)
    u64* inv[NTT_NB_PRIMES];            // psi^-bitrev(i)
    u64 nInv[NTT_NB_PRIMES];            // N^-1 mod q
};

//tables are built once per N (and never freed: there are only a few different N)
static const NttTables* getTables(const int32_t N) {
    static map<int32_t, NttTables*> cache;
    static mutex mtx;
    lock_guard<mutex> guard(mtx);
    map<int32_t, NttTables*>::iterator it = cache.find(N);
    if (it != cache.end()) return it->second;

    int32_t logn = 0;
    while ((1 << logn) < N) logn++;
    if ((1 << logn) != N || logn > NTT_MAX_LOGN) {   // N must be a power of 2, at most 2^22
        fprintf(stderr, "torusPolynomialMultNTT: unsupported N=%d\n", N);
        abort();
    }

    NttTables* t = new NttTables;
    t->N = N;
    u64* pw = new u64[N];
    for (int32_t p = 0; p < NTT_NB_PRIMES; p++) {
        const u64 q = NTT_Q[p];
        // primitive 2N-th root of unity: x^N = -1
        u64 psi = 0;
        for (u64 g = 2; ; g++) {
            const u64 x = powmod(g, (q - 1) / (2 * (u64) N), q);
            if (powmod(x, N, q) == q - 1) { psi = x; break; }
        }
        const u64 psiInv = powmod(psi, q - 2, q);
        t->fwd[p] = new u64[N];
        t->inv[p] = new u64[N];
        pw[0] = 1;
        for (int32_t i = 1; i < N; i++) pw[i] = mulmod(pw[i - 1], psi, q);
        for (int32_t i = 0; i < N; i++) t->fwd[p][i] = pw[bitrev(i, logn)];
        pw[0] = 1;
        for (int32_t i = 1; i < N; i++) pw[i] = mulmod(pw[i - 1], psiInv, q);
        for (int32_t i = 0; i < N; i++) t->inv[p][i] = pw[bitrev(i, logn)];
        t->nInv[p] = powmod(N, q - 2, q);
    }
    delete[] pw;
    cache[N] = t;
    return t;
}

// forward negacyclic NTT (Cooley-Tukey, input in natural order, output in bit-reversed order)
static void nttForward(u64* a, const int32_t N, const u64 q, const u64* psiRev) {
    int32_t t = N;
    for (int32_t m = 1; m < N; m <<= 1) {
        t >>= 1;
        for (int32_t i = 0; i < m; i++) {
            const int32_t j1 = 2 * i * t;
            const u64 S = psiRev[m + i];
            for (int32_t j = j1; j < j1 + t; j++) {
                const u64 U = a[j];
                const u64 V = mulmod(a[j + t], S, q);
                const u64 s = U + V;
                a[j] = (s >= q) ? s - q : s;
                a[j + t] = (U >= V) ? U - V : U + q - V;
            }
        }
    }
}

// inverse negacyclic NTT (Gentleman-Sande, input in bit-reversed order, output in natural order).
// This is the interpolation.
static void nttInverse(u64* a, const int32_t N, const u64 q, const u64* psiInvRev, const u64 nInv) {
    int32_t t = 1;
    for (int32_t m = N; m > 1; m >>= 1) {
        const int32_t h = m >> 1;
        int32_t j1 = 0;
        for (int32_t i = 0; i < h; i++) {
            const u64 S = psiInvRev[h + i];
            for (int32_t j = j1; j < j1 + t; j++) {
                const u64 U = a[j];
                const u64 V = a[j + t];
                const u64 s = U + V;
                a[j] = (s >= q) ? s - q : s;
                a[j + t] = mulmod((U >= V) ? U - V : U + q - V, S, q);
            }
            j1 += 2 * t;
        }
        t <<= 1;
    }
    for (int32_t i = 0; i < N; i++) a[i] = mulmod(a[i], nInv, q);
}

// signed value -> [0,q)
static u64 toResidue(int64_t v, u64 q) {
    int64_t r = v % (int64_t) q;
    return (u64) (r < 0 ? r + (int64_t) q : r);
}

// CRT (Garner) of the 3 residues, centered around 0, returned mod 2^32
static Torus32 crtToTorus(u64 r0, u64 r1, u64 r2) {
    static const u64 p0 = NTT_Q[0], p1 = NTT_Q[1], p2 = NTT_Q[2];
    static const u64 inv01 = powmod(p0 % p1, p1 - 2, p1);               // p0^-1 mod p1
    static const u64 inv012 = powmod((p0 * p1) % p2, p2 - 2, p2);       // (p0*p1)^-1 mod p2
    static const u128 P = (u128) p0 * p1 * p2;

    const u64 k1 = mulmod((r1 + p1 - r0 % p1) % p1, inv01, p1);
    const u64 partial = (r0 + k1 * p0) % p2;
    const u64 k2 = mulmod((r2 + p2 - partial) % p2, inv012, p2);
    u128 x = (u128) r0 + (u128) k1 * p0 + (u128) k2 * (p0 * p1);        // in [0,P)
    if (x > P / 2) x -= P;                                               // wraps mod 2^128: low 32 bits stay right
    return (Torus32) (uint32_t) x;
}

/**
 * Sum of several products with the NTT, with optional lazy interpolation.
 *
 * result = sum_{j<tau} poly1[j] * poly2[j]   (mod X^N+1, result is overwritten)
 *
 * lazy=false: each product goes through the 2 forward NTTs (per prime), the pointwise product,
 * the inverse NTT and the CRT, and the tau results are added.
 * lazy=true: the pointwise products are accumulated mod q in the evaluation space, and the
 * inverse NTT (the interpolation) and the CRT run only once on the sum, which saves
 * 3*(tau-1) inverse NTTs and tau-1 CRTs. Both modes give exactly the same result.
 */
EXPORT void torusPolynomialMultSumNTT(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy) {
    const int32_t N = result->N;
    const NttTables* tab = getTables(N);
    for (int32_t i = 0; i < N; i++) result->coefsT[i] = 0;

    const int32_t group = lazy ? tau : 1;
    u64* acc[NTT_NB_PRIMES];
    for (int32_t p = 0; p < NTT_NB_PRIMES; p++) acc[p] = new u64[N];
    u64* ta = new u64[N];
    u64* tb = new u64[N];

    for (int32_t start = 0; start < tau; start += group) {
        const int32_t g = (tau - start < group) ? tau - start : group;
        for (int32_t p = 0; p < NTT_NB_PRIMES; p++)
            for (int32_t i = 0; i < N; i++) acc[p][i] = 0;

        for (int32_t j = start; j < start + g; j++) {
            for (int32_t p = 0; p < NTT_NB_PRIMES; p++) {
                const u64 q = NTT_Q[p];
                for (int32_t i = 0; i < N; i++) {
                    ta[i] = toResidue(poly1[j]->coefs[i], q);
                    tb[i] = toResidue(poly2[j]->coefsT[i], q);   // Torus32 read as a signed integer
                }
                nttForward(ta, N, q, tab->fwd[p]);
                nttForward(tb, N, q, tab->fwd[p]);
                for (int32_t i = 0; i < N; i++) {                // accumulate in the evaluation space
                    const u64 s = acc[p][i] + mulmod(ta[i], tb[i], q);
                    acc[p][i] = (s >= q) ? s - q : s;
                }
            }
        }

        for (int32_t p = 0; p < NTT_NB_PRIMES; p++)
            nttInverse(acc[p], N, NTT_Q[p], tab->inv[p], tab->nInv[p]);
        for (int32_t i = 0; i < N; i++)
            result->coefsT[i] += crtToTorus(acc[0][i], acc[1][i], acc[2][i]);
    }

    for (int32_t p = 0; p < NTT_NB_PRIMES; p++) delete[] acc[p];
    delete[] ta;
    delete[] tb;
}

/** result = poly1 * poly2 mod X^N+1, exactly, via NTT */
EXPORT void torusPolynomialMultNTT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    torusPolynomialMultSumNTT(result, &poly1, &poly2, 1, false);
}
