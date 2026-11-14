/**
 * Toom-Cook interpolation variants, for a controlled comparison of interpolation strategies.
 *
 * All variants multiply an IntPolynomial by a TorusPolynomial mod X^N+1 exactly, in 64-bit
 * words (see multiplication.cpp for the exact-division / bit-loss argument). For a given K
 * the variants share the same evaluation code, workspace layout, base case and recomposition;
 * they differ only in the interpolation strategy (and, for the striding variant, in the data
 * layout):
 *
 *   TV_T3_OWN       classical Toom-3, own interpolation sequence (as in torusPolynomialMultToomCook3)
 *   TV_T3_BODRATO   classical Toom-3, Bodrato's interpolation sequence (WAIFI 2007), same points
 *   TV_T3_GENERIC   classical Toom-3, generic matrix interpolation (inverse Vandermonde), same points
 *   TV_T4_BODRATO   classical Toom-4, Bodrato's sequence (as in torusPolynomialMultToomCook4)
 *   TV_T4_GENERIC   classical Toom-4, generic matrix interpolation, same points
 *   TV_T4_STRIDED   striding Toom-4 (Bernstein; Bermudo Mera et al.): interleaved split, negacyclic
 *                   size-preserving sub-products, Bodrato's sequence, sequential recomposition
 *   TV_T5_GENERIC   classical Toom-5, generic matrix interpolation (as in torusPolynomialMultToomCook5)
 *
 * "Classical" means the operands are split into contiguous blocks, the sub-products have
 * 2h-1 coefficients, and the product is folded mod X^N+1 at the end. All variants evaluate and
 * interpolate coefficient by coefficient (vertical scanning).
 */
#include <cstdint>
#include <cstdlib>
#include "polynomials_arithmetic.h"

namespace {

typedef uint64_t u64;
typedef __int128 i128;

const u64 INV3 = 0xAAAAAAAAAAAAAAABULL, INV9 = 0x8E38E38E38E38E39ULL, INV15 = 0xEEEEEEEEEEEEEEEFULL;
static_assert(3ULL * 0xAAAAAAAAAAAAAAABULL == 1, "inverse of 3");
static_assert(9ULL * 0x8E38E38E38E38E39ULL == 1, "inverse of 9");
static_assert(15ULL * 0xEEEEEEEEEEEEEEEFULL == 1, "inverse of 15");

enum Interp { I_OWN3, I_BOD3, I_BOD4, I_GEN };

// ------------------------------------------------------------------ generic interpolation table

struct Frac { i128 n, d; };
i128 gcd128(i128 a, i128 b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { i128 t = a % b; a = b; b = t; } return a; }
Frac mk(i128 n, i128 d) { if (d < 0) { n = -n; d = -d; } i128 g = gcd128(n, d); if (g > 1) { n /= g; d /= g; } Frac f = { n, d }; return f; }
Frac fsub(Frac a, Frac b) { return mk(a.n * b.d - b.n * a.d, a.d * b.d); }
Frac fmul(Frac a, Frac b) { return mk(a.n * b.n, a.d * b.d); }
Frac fdiv(Frac a, Frac b) { return mk(a.n * b.d, a.d * b.n); }

struct GenTable {
    int32_t n;
    u64 M[9][9];
    u64 oddinv[9];
    int32_t shift[9];
    int32_t maxShift;
};

// points are homogeneous (p,q), i.e. the value p/q; the row of a product at p/q is p^j q^(n-1-j).
// The last row is the point at infinity.
GenTable buildGen(const int32_t K, const int32_t* p, const int32_t* q) {
    GenTable T;
    const int32_t n = 2 * K - 1;
    T.n = n;
    Frac a[9][18];
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < 2 * n; j++) a[i][j] = mk(0, 1);
    for (int32_t i = 0; i < n - 1; i++)
        for (int32_t j = 0; j < n; j++) {
            i128 v = 1;
            for (int32_t e = 0; e < j; e++) v *= p[i];
            for (int32_t e = 0; e < n - 1 - j; e++) v *= q[i];
            a[i][j] = mk(v, 1);
        }
    a[n - 1][n - 1] = mk(1, 1);
    for (int32_t i = 0; i < n; i++) a[i][n + i] = mk(1, 1);
    for (int32_t c = 0; c < n; c++) {               // Gauss-Jordan
        int32_t piv = c;
        while (a[piv][c].n == 0) piv++;
        for (int32_t j = 0; j < 2 * n; j++) { Frac t = a[c][j]; a[c][j] = a[piv][j]; a[piv][j] = t; }
        Frac pv = a[c][c];
        for (int32_t j = 0; j < 2 * n; j++) a[c][j] = fdiv(a[c][j], pv);
        for (int32_t i = 0; i < n; i++) {
            if (i == c || a[i][c].n == 0) continue;
            Frac f = a[i][c];
            for (int32_t j = 0; j < 2 * n; j++) a[i][j] = fsub(a[i][j], fmul(f, a[c][j]));
        }
    }
    T.maxShift = 0;
    for (int32_t i = 0; i < n; i++) {
        i128 den = 1;
        for (int32_t j = 0; j < n; j++) { const i128 d = a[i][n + j].d; den = den / gcd128(den, d) * d; }
        for (int32_t j = 0; j < n; j++) T.M[i][j] = (u64)(int64_t)(a[i][n + j].n * (den / a[i][n + j].d));
        int32_t s = 0;
        while (((den >> s) & 1) == 0) s++;
        const u64 odd = (u64)(den >> s);
        u64 x = odd;
        for (int32_t it = 0; it < 6; it++) x *= 2 - odd * x;   // odd^-1 mod 2^64
        T.oddinv[i] = x; T.shift[i] = s;
        if (s > T.maxShift) T.maxShift = s;
    }
    return T;
}

const int32_t P3[4] = { 0, 1, -1, -2 },          Q3[4] = { 1, 1, 1, 1 };
const int32_t P4[6] = { 0, 1, 1, -1, 2, -2 },    Q4[6] = { 1, 2, 1, 1, 1, 1 };   // 0, 1/2, 1, -1, 2, -2
const int32_t P5[8] = { 0, 1, -1, 2, -2, 3, -3, 4 }, Q5[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };

template<int K> const GenTable& gen() {
    static const GenTable T = buildGen(K, K == 3 ? P3 : (K == 4 ? P4 : P5), K == 3 ? Q3 : (K == 4 ? Q4 : Q5));
    return T;
}

// ------------------------------------------------------------------ evaluation (shared per K)

// e[j][t] = value of the pieces at point j, for one coefficient index t. a[i] = i-th piece coefficient.
inline void eval3(u64* const* e, int32_t t, u64 a0, u64 a1, u64 a2) {
    const u64 s = a0 + a2;
    const u64 m = s - a1;
    e[0][t] = a0; e[1][t] = s + a1; e[2][t] = m; e[3][t] = ((m + a2) << 1) - a0; e[4][t] = a2;
}
inline void eval4(u64* const* e, int32_t t, u64 a0, u64 a1, u64 a2, u64 a3) {
    const u64 s02 = a0 + a2, s13 = a1 + a3, ev = a0 + (a2 << 2), od = (a1 << 1) + (a3 << 3);
    e[0][t] = a0; e[1][t] = (a0 << 3) + (a1 << 2) + (a2 << 1) + a3;
    e[2][t] = s02 + s13; e[3][t] = s02 - s13; e[4][t] = ev + od; e[5][t] = ev - od; e[6][t] = a3;
}
inline void eval5(u64* const* e, int32_t t, const u64* a) {      // Horner at 0, +-1, +-2, +-3, 4, inf
    for (int32_t j = 0; j < 8; j++) {
        const u64 x = (u64)(int64_t) P5[j];
        u64 v = a[4];
        for (int32_t i = 3; i >= 0; i--) v = v * x + a[i];
        e[j][t] = v;
    }
    e[8][t] = a[4];
}

// ------------------------------------------------------------------ interpolation strategies

// own Toom-3 sequence (identical to multiplication.cpp); r = {r(0), r(1), r(-1), r(-2), r(inf)} -> c0..c4
inline void interpOwn3(u64* const* r, int32_t t) {
    const u64 c0 = r[0][t], c4 = r[4][t], r1 = r[1][t], rm = r[2][t], r2 = r[3][t];
    const u64 c2 = ((r1 + rm) >> 1) - c0 - c4;
    const u64 s = (r1 - rm) >> 1;                       // c1 + c3
    const u64 u = (c0 + (c2 << 2) + (c4 << 4) - r2) >> 1; // c1 + 4 c3
    const u64 c3 = (u - s) * INV3;
    r[1][t] = s - c3; r[2][t] = c2; r[3][t] = c3;
}
// Bodrato's Toom-3 sequence for the points 0, 1, -1, -2, inf (Bodrato, WAIFI 2007)
inline void interpBod3(u64* const* r, int32_t t) {
    const u64 r0 = r[0][t], ri = r[4][t];
    u64 r1 = r[1][t], rm1 = r[2][t], r3 = r[3][t], r2;
    r3 = (r3 - r1) * INV3;              // (r(-2) - r(1)) / 3
    r1 = (r1 - rm1) >> 1;               // (r(1) - r(-1)) / 2
    r2 = rm1 - r0;                      // r(-1) - r(0)
    r3 = ((r2 - r3) >> 1) + (ri << 1);  // (r2 - r3) / 2 + 2 r(inf)
    r2 = r2 + r1 - ri;
    r1 = r1 - r3;
    r[1][t] = r1; r[2][t] = r2; r[3][t] = r3;
}
// Bodrato's Toom-4 sequence (identical to toom4_interp in multiplication.cpp and to Algorithm 11 of
// Bermudo Mera's thesis); r = {0, 1/2 (x64), 1, -1, 2, -2, inf} -> c0..c6
inline void interpBod4(u64* const* rr, int32_t t) {
    u64 r0 = rr[0][t], r1 = rr[1][t], r2 = rr[2][t], r3 = rr[3][t], r4 = rr[4][t], r5 = rr[5][t], r6 = rr[6][t];
    r1 = r1 + r4;
    r5 = r5 - r4;
    r3 = (r3 - r2) >> 1;
    r4 = r4 - r0;
    r4 = r4 - (r6 << 6);
    r4 = (r4 << 1) + r5;
    r2 = r2 + r3;
    r1 = r1 - 65 * r2;
    r2 = r2 - r6;
    r2 = r2 - r0;
    r1 = r1 + 45 * r2;
    r4 = ((r4 - (r2 << 3)) * INV3) >> 3;    // /24
    r5 = r5 + r1;
    r1 = ((r1 + (r3 << 4)) * INV9) >> 1;    // /18
    r3 = -(r3 + r1);
    r5 = ((r1 * 30 - r5) * INV15) >> 2;     // /60
    r2 = r2 - r4;
    r1 = r1 - r5;
    // with the inputs in the order 0, 1/2, 1, -1, 2, -2, inf the sequence leaves c_i in r_i
    // (Algorithm 11 of Bermudo Mera's thesis uses the mirrored order, infinity first)
    rr[1][t] = r1; rr[2][t] = r2; rr[3][t] = r3; rr[4][t] = r4; rr[5][t] = r5;
}
// generic: c_i = ((sum_j M_ij r_j) * odd_i^-1) >> s_i
template<int K> inline void interpGen(u64* const* r, int32_t t, const GenTable& T) {
    const int32_t n = 2 * K - 1;
    u64 v[9];
    for (int32_t j = 0; j < n; j++) v[j] = r[j][t];
    for (int32_t i = 0; i < n; i++) {
        u64 acc = 0;
        for (int32_t j = 0; j < n; j++) acc += T.M[i][j] * v[j];
        r[i][t] = (acc * T.oddinv[i]) >> T.shift[i];
    }
}

// bits of the 64-bit word lost per recursion level
int32_t lossOf(const int32_t K, const Interp I) {
    if (I == I_OWN3 || I == I_BOD3) return 1;
    if (I == I_BOD4) return 4;
    return K == 3 ? gen<3>().maxShift : (K == 4 ? gen<4>().maxShift : gen<5>().maxShift);
}

// ------------------------------------------------------------------ classical recursion

void naivePlain(u64* R, const u64* A, const u64* B, const int32_t n) {
    for (int32_t i = 0; i < 2 * n - 1; i++) R[i] = 0;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++) R[i + j] += A[i] * B[j];
}

template<int K>
size_t wsClassical(const int32_t size, const int32_t base) {
    if (size <= base) return 0;
    const int32_t h = (size + K - 1) / K, n = 2 * K - 1;
    return (size_t)(2 * K + 4 * n) * h + wsClassical<K>(h, base);
}

template<int K, Interp I>
void recClassical(u64* R, const u64* A, const u64* B, const int32_t size, u64* ws, const int32_t base) {
    if (size <= base) { naivePlain(R, A, B, size); return; }
    const int32_t n = 2 * K - 1, h = (size + K - 1) / K, rl = 2 * h - 1;
    u64* Ap = ws; ws += K * h;
    u64* Bp = ws; ws += K * h;
    u64 *eA[9], *eB[9], *r[9];
    for (int32_t j = 0; j < n; j++) { eA[j] = ws; ws += h; }
    for (int32_t j = 0; j < n; j++) { eB[j] = ws; ws += h; }
    for (int32_t j = 0; j < n; j++) { r[j] = ws; ws += 2 * h; }
    for (int32_t i = 0; i < K * h; i++) { Ap[i] = (i < size) ? A[i] : 0; Bp[i] = (i < size) ? B[i] : 0; }

    for (int32_t t = 0; t < h; t++) {
        if (K == 3) { eval3(eA, t, Ap[t], Ap[h + t], Ap[2 * h + t]); eval3(eB, t, Bp[t], Bp[h + t], Bp[2 * h + t]); }
        else if (K == 4) { eval4(eA, t, Ap[t], Ap[h + t], Ap[2 * h + t], Ap[3 * h + t]); eval4(eB, t, Bp[t], Bp[h + t], Bp[2 * h + t], Bp[3 * h + t]); }
        else {
            u64 pa[5], pb[5];
            for (int32_t i = 0; i < 5; i++) { pa[i] = Ap[i * h + t]; pb[i] = Bp[i * h + t]; }
            eval5(eA, t, pa); eval5(eB, t, pb);
        }
    }
    for (int32_t j = 0; j < n; j++) recClassical<K, I>(r[j], eA[j], eB[j], h, ws, base);

    if (I == I_OWN3)      for (int32_t t = 0; t < rl; t++) interpOwn3(r, t);
    else if (I == I_BOD3) for (int32_t t = 0; t < rl; t++) interpBod3(r, t);
    else if (I == I_BOD4) for (int32_t t = 0; t < rl; t++) interpBod4(r, t);
    else { const GenTable& T = gen<K>(); for (int32_t t = 0; t < rl; t++) interpGen<K>(r, t, T); }

    const int32_t Rlen = 2 * size - 1;
    for (int32_t i = 0; i < Rlen; i++) R[i] = 0;
    for (int32_t k = 0; k < n; k++)
        for (int32_t i = 0; i < rl && k * h + i < Rlen; i++) R[k * h + i] += r[k][i];
}

template<int K>
int32_t depthClassical(int32_t size, const int32_t base) {
    int32_t d = 0;
    while (size > base) { size = (size + K - 1) / K; d++; }
    return d;
}

template<int K, Interp I>
int32_t multClassical(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2, const int32_t base) {
    const int32_t N = poly1->N;
    if (depthClassical<K>(N, base) * lossOf(K, I) > 32) return -1;   // precision budget (64 - 32 bits)
    u64* buf = new u64[4 * (size_t)N + wsClassical<K>(N, base)];
    u64* A = buf; u64* B = A + N; u64* R = B + N; u64* ws = R + 2 * N;
    for (int32_t i = 0; i < N; i++) { A[i] = (u64)(int64_t) poly1->coefs[i]; B[i] = (u64)(int64_t) poly2->coefsT[i]; }
    recClassical<K, I>(R, A, B, N, ws, base);
    for (int32_t i = 0; i < N - 1; i++) result->coefsT[i] = (Torus32)(uint32_t)(R[i] - R[N + i]);
    result->coefsT[N - 1] = (Torus32)(uint32_t) R[N - 1];
    delete[] buf;
    return 0;
}

// ------------------------------------------------------------------ striding Toom-4 (negacyclic)

// R = A * B mod (y^n + 1), all of size n
void naiveNega(u64* R, const u64* A, const u64* B, const int32_t n) {
    for (int32_t i = 0; i < n; i++) {
        u64 s = 0;
        for (int32_t j = 0; j <= i; j++) s += A[j] * B[i - j];
        for (int32_t j = i + 1; j < n; j++) s -= A[j] * B[n + i - j];
        R[i] = s;
    }
}

size_t wsStrided(const int32_t n, const int32_t base) {
    if (n <= base || (n & 3)) return 0;
    const int32_t r = n / 4;
    return (size_t) 21 * r + wsStrided(r, base);
}

int32_t depthStrided(int32_t n, const int32_t base) {
    int32_t d = 0;
    while (n > base && !(n & 3)) { n /= 4; d++; }
    return d;
}

// a(x) = sum_{i<4} x^i A_i(y), y = x^4, A_i(y) = sum_t a[4t+i] y^t: every sub-product is a negacyclic
// product of size n/4, and recomposition writes the result already reduced mod x^n + 1
void recStrided(u64* R, const u64* A, const u64* B, const int32_t n, u64* ws, const int32_t base) {
    if (n <= base || (n & 3)) { naiveNega(R, A, B, n); return; }
    const int32_t r = n / 4;
    u64 *eA[7], *eB[7], *w[7];
    for (int32_t j = 0; j < 7; j++) { eA[j] = ws; ws += r; }
    for (int32_t j = 0; j < 7; j++) { eB[j] = ws; ws += r; }
    for (int32_t j = 0; j < 7; j++) { w[j] = ws; ws += r; }
    for (int32_t t = 0; t < r; t++) {           // sequential access: a[4t], a[4t+1], a[4t+2], a[4t+3]
        eval4(eA, t, A[4 * t], A[4 * t + 1], A[4 * t + 2], A[4 * t + 3]);
        eval4(eB, t, B[4 * t], B[4 * t + 1], B[4 * t + 2], B[4 * t + 3]);
    }
    for (int32_t j = 0; j < 7; j++) recStrided(w[j], eA[j], eB[j], r, ws, base);
    for (int32_t t = 0; t < r; t++) interpBod4(w, t);
    // C(x) = sum_{j<7} C_j(y) x^j, and x^j = y x^(j-4) for j >= 4 (y C(y) wraps negatively)
    for (int32_t t = 0; t < r; t++) {
        for (int32_t j = 0; j < 3; j++) R[4 * t + j] = w[j][t] + (t > 0 ? w[j + 4][t - 1] : (u64)0 - w[j + 4][r - 1]);
        R[4 * t + 3] = w[3][t];
    }
}

int32_t multStrided(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2, const int32_t base) {
    const int32_t N = poly1->N;
    if (depthStrided(N, base) * lossOf(4, I_BOD4) > 32) return -1;
    u64* buf = new u64[3 * (size_t)N + wsStrided(N, base)];
    u64* A = buf; u64* B = A + N; u64* R = B + N; u64* ws = R + N;
    for (int32_t i = 0; i < N; i++) { A[i] = (u64)(int64_t) poly1->coefs[i]; B[i] = (u64)(int64_t) poly2->coefsT[i]; }
    recStrided(R, A, B, N, ws, base);
    for (int32_t i = 0; i < N; i++) result->coefsT[i] = (Torus32)(uint32_t) R[i];
    delete[] buf;
    return 0;
}

} // namespace

EXPORT int32_t torusPolynomialMultToomVariant(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2, const int32_t variant, const int32_t base) {
    switch (variant) {
        case TV_T3_OWN:     return multClassical<3, I_OWN3>(result, poly1, poly2, base);
        case TV_T3_BODRATO: return multClassical<3, I_BOD3>(result, poly1, poly2, base);
        case TV_T3_GENERIC: return multClassical<3, I_GEN>(result, poly1, poly2, base);
        case TV_T4_BODRATO: return multClassical<4, I_BOD4>(result, poly1, poly2, base);
        case TV_T4_GENERIC: return multClassical<4, I_GEN>(result, poly1, poly2, base);
        case TV_T4_STRIDED: return multStrided(result, poly1, poly2, base);
        case TV_T5_GENERIC: return multClassical<5, I_GEN>(result, poly1, poly2, base);
        default: return -2;
    }
}

EXPORT int32_t toomVariantDepth(const int32_t variant, const int32_t N, const int32_t base) {
    switch (variant) {
        case TV_T3_OWN: case TV_T3_BODRATO: case TV_T3_GENERIC: return depthClassical<3>(N, base);
        case TV_T4_BODRATO: case TV_T4_GENERIC: return depthClassical<4>(N, base);
        case TV_T4_STRIDED: return depthStrided(N, base);
        case TV_T5_GENERIC: return depthClassical<5>(N, base);
        default: return -1;
    }
}
