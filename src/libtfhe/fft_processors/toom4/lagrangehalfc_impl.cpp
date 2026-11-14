/**
 * Toom-4 striding transform processor (see lagrangehalfc_impl.h): replaces the FFT in the
 * LagrangeHalfCPolynomial interface, hence in the whole bootstrapping (key conversion,
 * external product, blind rotation), without touching the bootstrapping code.
 */
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <memory>
#include <vector>
#include "lagrangehalfc_impl.h"
#include "lagrangehalfc_arithmetic.h"

typedef uint64_t u64;

namespace {

const u64 INV3 = 0xAAAAAAAAAAAAAAABULL, INV9 = 0x8E38E38E38E38E39ULL, INV15 = 0xEEEEEEEEEEEEEEEFULL;
static_assert(3ULL * 0xAAAAAAAAAAAAAAABULL == 1, "inverse of 3");
static_assert(9ULL * 0x8E38E38E38E38E39ULL == 1, "inverse of 9");
static_assert(15ULL * 0xEEEEEEEEEEEEEEEFULL == 1, "inverse of 15");

// evaluation of a0 + a1 x + a2 x^2 + a3 x^3 at 0, 1/2 (times 8), 1, -1, 2, -2, inf
// (identical to eval4 in toom-variants.cpp)
inline void eval4(u64* const* e, int32_t t, u64 a0, u64 a1, u64 a2, u64 a3) {
    const u64 s02 = a0 + a2, s13 = a1 + a3, ev = a0 + (a2 << 2), od = (a1 << 1) + (a3 << 3);
    e[0][t] = a0; e[1][t] = (a0 << 3) + (a1 << 2) + (a2 << 1) + a3;
    e[2][t] = s02 + s13; e[3][t] = s02 - s13; e[4][t] = ev + od; e[5][t] = ev - od; e[6][t] = a3;
}

// Bodrato's Toom-4 interpolation sequence (identical to interpBod4 in toom-variants.cpp);
// rr = {r(0), 64 r(1/2), r(1), r(-1), r(2), r(-2), r(inf)} -> c0..c6 in place
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
    rr[1][t] = r1; rr[2][t] = r2; rr[3][t] = r3; rr[4][t] = r4; rr[5][t] = r5;
}

// size of one transformed polynomial of size n at a given depth
inline size_t treeWords(int32_t n, int32_t depth) {
    size_t w = (size_t) n;
    for (int32_t d = 0; d < depth; d++) w = w / 4 * 7;
    return w;
}

// evaluation down the tree: in (n words) -> out (treeWords(n, d) words)
void fwdRec(u64* out, const u64* in, const int32_t n, const int32_t d, u64* ws) {
    if (d == 0) { memcpy(out, in, n * sizeof(u64)); return; }
    const int32_t r = n / 4;
    u64* e[7];
    for (int32_t j = 0; j < 7; j++) e[j] = ws + j * r;
    for (int32_t t = 0; t < r; t++)
        eval4(e, t, in[4 * t], in[4 * t + 1], in[4 * t + 2], in[4 * t + 3]);
    const size_t cs = treeWords(r, d - 1);
    for (int32_t j = 0; j < 7; j++) fwdRec(out + j * cs, e[j], r, d - 1, ws + 7 * r);
}

// interpolation + recomposition up the tree: tree -> R (n words, already reduced mod x^n + 1)
void bwdRec(u64* R, const u64* tree, const int32_t n, const int32_t d, u64* ws) {
    if (d == 0) { memcpy(R, tree, n * sizeof(u64)); return; }
    const int32_t r = n / 4;
    u64* w[7];
    for (int32_t j = 0; j < 7; j++) w[j] = ws + j * r;
    const size_t cs = treeWords(r, d - 1);
    for (int32_t j = 0; j < 7; j++) bwdRec(w[j], tree + j * cs, r, d - 1, ws + 7 * r);
    for (int32_t t = 0; t < r; t++) interpBod4(w, t);
    // C(x) = sum_{j<7} C_j(y) x^j with y = x^4, and x^j = y x^(j-4) for j >= 4 (y C(y) wraps negatively)
    for (int32_t t = 0; t < r; t++) {
        for (int32_t j = 0; j < 3; j++) R[4 * t + j] = w[j][t] + (t > 0 ? w[j + 4][t - 1] : (u64) 0 - w[j + 4][r - 1]);
        R[4 * t + 3] = w[3][t];
    }
}

// acc += A * B mod (y^m + 1)   (schoolbook, size m)
inline void negaMacAdd(u64* __restrict acc, const u64* __restrict A, const u64* __restrict B, const int32_t m) {
    for (int32_t i = 0; i < m; i++) {
        u64 s = 0;
        for (int32_t j = 0; j <= i; j++) s += A[j] * B[i - j];
        for (int32_t j = i + 1; j < m; j++) s -= A[j] * B[m + i - j];
        acc[i] += s;
    }
}
inline void negaMacSub(u64* __restrict acc, const u64* __restrict A, const u64* __restrict B, const int32_t m) {
    for (int32_t i = 0; i < m; i++) {
        u64 s = 0;
        for (int32_t j = 0; j <= i; j++) s += A[j] * B[i - j];
        for (int32_t j = i + 1; j < m; j++) s -= A[j] * B[m + i - j];
        acc[i] -= s;
    }
}

} // namespace


Toom4StridedProcessor::Toom4StridedProcessor(const int32_t N) :
    N(N),
    depth([N]() { int32_t d = 0, n = N; while (n > TOOM4_LEAF && (n & 3) == 0 && d < 8) { n /= 4; d++; } return d; }()),
    leaf(N >> (2 * depth)),
    nleaves([this]() { int32_t p = 1; for (int32_t d = 0; d < depth; d++) p *= 7; return p; }()),
    size((size_t) nleaves * (size_t) leaf)
{
    // 4 bits lost per interpolation level: the 32 low bits (the torus) must stay exact
    assert(4 * depth <= 32);
    // scratch: N words for a coefficient-domain operand, 3N for the recursion (7n/4 per level < 2.4 N),
    // then 3 transformed polynomials for the format conversions and the constant helpers
    wsWords = 4 * (size_t) N + 3 * size;
    one = new u64[size];
    u64* a = workspace();
    memset(a, 0, N * sizeof(u64));
    a[0] = 1;
    forward(one, a);
}

const Toom4StridedProcessor* Toom4StridedProcessor::get(const int32_t N) {
    static std::mutex mtx;
    static std::map<int32_t, std::unique_ptr<Toom4StridedProcessor> > procs;
    std::lock_guard<std::mutex> lock(mtx);
    std::unique_ptr<Toom4StridedProcessor>& p = procs[N];
    if (!p) p.reset(new Toom4StridedProcessor(N));
    return p.get();
}

u64* Toom4StridedProcessor::workspace() const {
    thread_local std::vector<u64> buf;
    if (buf.size() < wsWords) buf.resize(wsWords);
    return buf.data();
}

void Toom4StridedProcessor::forward(u64* tree, const u64* a) const {
    fwdRec(tree, a, N, depth, workspace() + N);
}

void Toom4StridedProcessor::backward(u64* r, const u64* tree) const {
    bwdRec(r, tree, N, depth, workspace() + N);
}

void Toom4StridedProcessor::mulAcc(u64* acc, const u64* a, const u64* b, const int32_t sign) const {
    if (sign > 0) for (int32_t l = 0; l < nleaves; l++) negaMacAdd(acc + l * leaf, a + l * leaf, b + l * leaf, leaf);
    else          for (int32_t l = 0; l < nleaves; l++) negaMacSub(acc + l * leaf, a + l * leaf, b + l * leaf, leaf);
}

// a = a * 1: the leaf products with the transform of 1 give the PRODUCT layout of a
void Toom4StridedProcessor::toProduct(u64* dst, const u64* src) const {
    assert(dst != src);
    memset(dst, 0, size * sizeof(u64));
    mulAcc(dst, src, one, +1);
}

void Toom4StridedProcessor::toOperand(u64* dst, const u64* src) const {
    u64* a = workspace();          // coefficient buffer
    backward(a, src);
    forward(dst, a);
}


LagrangeHalfCPolynomial_IMPL::LagrangeHalfCPolynomial_IMPL(const int32_t N) {
    proc = Toom4StridedProcessor::get(N);
    coefs = new u64[proc->size + 1];
    memset(coefs, 0, proc->size * sizeof(u64));
    format() = TOOM4_ZERO;
}

LagrangeHalfCPolynomial_IMPL::~LagrangeHalfCPolynomial_IMPL() {
    delete[] coefs;
}

static inline LagrangeHalfCPolynomial_IMPL* impl(LagrangeHalfCPolynomial* p) { return (LagrangeHalfCPolynomial_IMPL*) p; }
static inline const LagrangeHalfCPolynomial_IMPL* impl(const LagrangeHalfCPolynomial* p) { return (const LagrangeHalfCPolynomial_IMPL*) p; }

// the OPERAND layout of p: p itself, or a converted copy in the scratch tree "slot"
// (also copied when p aliases the destination "avoid", which is about to be overwritten)
static const u64* operandOf(const LagrangeHalfCPolynomial_IMPL* p, const int32_t slot, const u64* avoid) {
    const Toom4StridedProcessor* proc = p->proc;
    if (p->format() == TOOM4_PRODUCT) {
        u64* t = proc->tmpTree(slot);
        proc->toOperand(t, p->coefs);
        return t;
    }
    if (p->coefs == avoid) {
        u64* t = proc->tmpTree(slot);
        memcpy(t, p->coefs, proc->size * sizeof(u64));
        return t;
    }
    return p->coefs;
}

// brings an accumulator to the PRODUCT layout (in place)
static void makeProduct(LagrangeHalfCPolynomial_IMPL* acc) {
    if (acc->format() == TOOM4_OPERAND) {
        const Toom4StridedProcessor* proc = acc->proc;
        u64* t = proc->tmpTree(2);
        proc->toProduct(t, acc->coefs);
        memcpy(acc->coefs, t, proc->size * sizeof(u64));
    }
    acc->format() = TOOM4_PRODUCT;
}


//initialize the LagrangeHalfCPolynomial structure
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj, const int32_t N) {
    new(obj) LagrangeHalfCPolynomial_IMPL(N);
}
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N) {
    for (int32_t i = 0; i < nbelts; i++) new(obj + i) LagrangeHalfCPolynomial_IMPL(N);
}

//destroys the LagrangeHalfCPolynomial structure
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    impl(obj)->~LagrangeHalfCPolynomial_IMPL();
}
EXPORT void destroy_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    for (int32_t i = 0; i < nbelts; i++) impl(obj + i)->~LagrangeHalfCPolynomial_IMPL();
}


/**
 * transforms (names kept from the FFT interface: "ifft" = to the transform domain)
 */
EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const int32_t N = r->proc->N;
    assert(p->N == N);
    u64* a = r->proc->workspace();
    for (int32_t i = 0; i < N; i++) a[i] = (u64)(int64_t) p->coefs[i];
    r->proc->forward(r->coefs, a);
    r->format() = TOOM4_OPERAND;
}

EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const int32_t N = r->proc->N;
    assert(p->N == N);
    u64* a = r->proc->workspace();
    for (int32_t i = 0; i < N; i++) a[i] = (u64)(int64_t) p->coefsT[i];
    r->proc->forward(r->coefs, a);
    r->format() = TOOM4_OPERAND;
}

EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    const LagrangeHalfCPolynomial_IMPL* s = impl(p);
    const Toom4StridedProcessor* proc = s->proc;
    const int32_t N = proc->N;
    assert(result->N == N);
    if (s->format() == TOOM4_ZERO) {
        for (int32_t i = 0; i < N; i++) result->coefsT[i] = 0;
        return;
    }
    const u64* tree = s->coefs;
    if (s->format() == TOOM4_OPERAND) {
        u64* t = proc->tmpTree(0);
        proc->toProduct(t, tree);
        tree = t;
    }
    u64* R = proc->workspace();
    proc->backward(R, tree);
    for (int32_t i = 0; i < N; i++) result->coefsT[i] = (Torus32)(uint32_t) R[i];   // the torus = 32 low bits
}


//MISC OPERATIONS

/** sets to zero */
EXPORT void LagrangeHalfCPolynomialClear(LagrangeHalfCPolynomial* reps) {
    LagrangeHalfCPolynomial_IMPL* r = impl(reps);
    memset(r->coefs, 0, r->proc->size * sizeof(u64));
    r->format() = TOOM4_ZERO;
}

// OPERAND transform of the polynomial whose only nonzero coefficients are given
static void transformSparse(const Toom4StridedProcessor* proc, u64* tree, const int32_t nb, const int32_t* idx, const u64* val) {
    u64* a = proc->workspace();
    memset(a, 0, proc->N * sizeof(u64));
    for (int32_t i = 0; i < nb; i++) a[idx[i]] += val[i];
    proc->forward(tree, a);
}

/** sets to this torus32 constant */
EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const int32_t idx = 0;
    const u64 val = (u64)(int64_t) mu;
    transformSparse(r->proc, r->coefs, 1, &idx, &val);
    r->format() = TOOM4_OPERAND;
}

EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const Toom4StridedProcessor* proc = r->proc;
    u64* c = proc->tmpTree(1);
    const int32_t idx = 0;
    const u64 val = (u64)(int64_t) mu;
    transformSparse(proc, c, 1, &idx, &val);
    if (r->format() == TOOM4_PRODUCT) {
        u64* t = proc->tmpTree(0);
        proc->toProduct(t, c);
        c = t;
    } else {
        r->format() = TOOM4_OPERAND;   // ZERO + OPERAND
    }
    for (size_t i = 0; i < proc->size; i++) r->coefs[i] += c[i];
}

/** sets to X^ai-1 (0 <= ai < 2N) */
EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int32_t ai) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const int32_t N = r->proc->N;
    assert(ai >= 0 && ai < 2 * N);
    const int32_t idx[2] = { ai < N ? ai : ai - N, 0 };
    const u64 val[2] = { ai < N ? (u64) 1 : (u64) -1, (u64) -1 };
    transformSparse(r->proc, r->coefs, 2, idx, val);
    r->format() = TOOM4_OPERAND;
}


/** termwise multiplication in the transform domain (one negacyclic product per leaf) */
EXPORT void LagrangeHalfCPolynomialMul(
        LagrangeHalfCPolynomial* result,
        const LagrangeHalfCPolynomial* a,
        const LagrangeHalfCPolynomial* b) {
    LagrangeHalfCPolynomial_IMPL* r = impl(result);
    const Toom4StridedProcessor* proc = r->proc;
    const u64* pa = operandOf(impl(a), 0, r->coefs);
    const u64* pb = operandOf(impl(b), 1, r->coefs);
    memset(r->coefs, 0, proc->size * sizeof(u64));
    proc->mulAcc(r->coefs, pa, pb, +1);
    r->format() = TOOM4_PRODUCT;
}

EXPORT void LagrangeHalfCPolynomialAddTo(
        LagrangeHalfCPolynomial* accum,
        const LagrangeHalfCPolynomial* a) {
    LagrangeHalfCPolynomial_IMPL* r = impl(accum);
    const LagrangeHalfCPolynomial_IMPL* s = impl(a);
    const Toom4StridedProcessor* proc = r->proc;
    if (s->format() == TOOM4_ZERO) return;
    const u64* aa = s->coefs;
    if (r->format() == TOOM4_ZERO) {
        r->format() = s->format();
    } else if (r->format() != s->format()) {   // mixed layouts: add in the PRODUCT layout
        if (s->format() == TOOM4_OPERAND) {
            u64* t = proc->tmpTree(0);
            proc->toProduct(t, aa);
            aa = t;
        }
        makeProduct(r);
    }
    for (size_t i = 0; i < proc->size; i++) r->coefs[i] += aa[i];
}

/** termwise multiplication and addTo in the transform domain */
EXPORT void LagrangeHalfCPolynomialAddMul(
        LagrangeHalfCPolynomial* accum,
        const LagrangeHalfCPolynomial* a,
        const LagrangeHalfCPolynomial* b) {
    LagrangeHalfCPolynomial_IMPL* r = impl(accum);
    const u64* pa = operandOf(impl(a), 0, r->coefs);
    const u64* pb = operandOf(impl(b), 1, r->coefs);
    makeProduct(r);
    r->proc->mulAcc(r->coefs, pa, pb, +1);
}

EXPORT void LagrangeHalfCPolynomialSubMul(
        LagrangeHalfCPolynomial* accum,
        const LagrangeHalfCPolynomial* a,
        const LagrangeHalfCPolynomial* b) {
    LagrangeHalfCPolynomial_IMPL* r = impl(accum);
    const u64* pa = operandOf(impl(a), 0, r->coefs);
    const u64* pb = operandOf(impl(b), 1, r->coefs);
    makeProduct(r);
    r->proc->mulAcc(r->coefs, pa, pb, -1);
}
