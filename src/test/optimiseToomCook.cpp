// Finds the best Toom-Cook-3 recursion cutoff (the size at or below which the
// naive product is used instead of recursing). In multiplication.cpp this is the
// compile-time constant TOOM_BASE (currently 12).
//
// The Toom-3 code below is a copy of toom3_rec / ToomCook3_aux from multiplication.cpp,
// except that the cutoff is a runtime parameter, so we can sweep it in one run.
// The wrapper does the same work as torusPolynomialMultToomCook3 (alloc, 64-bit lift,
// product, truncation to 32 bits, reduction mod X^N+1), so the timings are end-to-end.
// Every cutoff is checked against the naive product mod X^N+1 (exact equality).
//
// Build (from the tfhe directory):
//   g++ -O3 -march=native -std=c++11 src/test/optimiseToomCook.cpp -o optimiseToomCook
// Run:
//   ./optimiseToomCook                 # N=1024, cutoffs 1..200
//   ./optimiseToomCook 2048 100 2000   # N, max cutoff, repetitions per cutoff
//
// Keep -O3 (and ideally the flags the library is built with); the optimum depends on it.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <vector>
#include <algorithm>

using namespace std;

typedef int32_t Torus32;
typedef uint64_t u64;

// ---------------------------------------------------------------- Toom-3 (parameterised)

static void toom3_naive64(u64* R, const u64* A, const u64* B, const int32_t n) {
    for (int32_t i = 0; i < 2*n-1; i++) R[i] = 0;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            R[i+j] += A[i]*B[j];
}

static size_t toom3_ws(const int32_t size, const int32_t base) {
    if (size <= base) return 0;
    const int32_t h = (size+2)/3;
    return 22*(size_t)h + toom3_ws(h, base);
}

static void toom3_rec(u64* R, const u64* A, const u64* B, const int32_t size, u64* ws, const int32_t base) {
    if (size <= base) {
        toom3_naive64(R, A, B, size);
        return;
    }
    const int32_t h = (size+2)/3;
    const int32_t rl = 2*h-1;
    const u64 inv3 = 0xAAAAAAAAAAAAAAABULL;

    u64* Ap  = ws; ws += 3*h;
    u64* Bp  = ws; ws += 3*h;
    u64* eA1 = ws; ws += h;  u64* eAm = ws; ws += h;  u64* eA2 = ws; ws += h;
    u64* eB1 = ws; ws += h;  u64* eBm = ws; ws += h;  u64* eB2 = ws; ws += h;
    u64* r0  = ws; ws += 2*h;
    u64* r1  = ws; ws += 2*h;
    u64* rm  = ws; ws += 2*h;
    u64* r2  = ws; ws += 2*h;
    u64* rinf= ws; ws += 2*h;

    for (int32_t i = 0; i < 3*h; i++) {
        Ap[i] = (i < size) ? A[i] : 0;
        Bp[i] = (i < size) ? B[i] : 0;
    }
    const u64 *a0 = Ap, *a1 = Ap+h, *a2 = Ap+2*h;
    const u64 *b0 = Bp, *b1 = Bp+h, *b2 = Bp+2*h;

    for (int32_t i = 0; i < h; i++) {
        u64 e = a0[i] + a2[i];
        eA1[i] = e + a1[i];
        eAm[i] = e - a1[i];
        eA2[i] = ((eAm[i] + a2[i]) << 1) - a0[i];
        e = b0[i] + b2[i];
        eB1[i] = e + b1[i];
        eBm[i] = e - b1[i];
        eB2[i] = ((eBm[i] + b2[i]) << 1) - b0[i];
    }

    toom3_rec(r0,   a0,  b0,  h, ws, base);
    toom3_rec(r1,   eA1, eB1, h, ws, base);
    toom3_rec(rm,   eAm, eBm, h, ws, base);
    toom3_rec(r2,   eA2, eB2, h, ws, base);
    toom3_rec(rinf, a2,  b2,  h, ws, base);

    for (int32_t i = 0; i < rl; i++) {
        u64 c0 = r0[i], c4 = rinf[i];
        u64 c2 = ((r1[i] + rm[i]) >> 1) - c0 - c4;
        u64 t  = (r1[i] - rm[i]) >> 1;
        u64 u  = (c0 + (c2 << 2) + (c4 << 4) - r2[i]) >> 1;
        u64 c3 = (u - t) * inv3;
        u64 c1 = t - c3;
        r1[i] = c2; rm[i] = c1; r2[i] = c3;
    }

    const int32_t Rlen = 2*size-1;
    for (int32_t i = 0; i < Rlen; i++) R[i] = 0;
    const u64* c[5] = { r0, rm, r1, r2, rinf };
    for (int32_t k = 0; k < 5; k++)
        for (int32_t i = 0; i < rl && k*h+i < Rlen; i++)
            R[k*h+i] += c[k][i];
}

// same work as torusPolynomialMultToomCook3, with a runtime cutoff
static void toomMult(Torus32* result, const int32_t* a, const Torus32* b, const int32_t N, const int32_t base) {
    Torus32* R = new Torus32[2*N-1];
    u64* ws = new u64[4*(size_t)N + toom3_ws(N, base)];
    u64* A64 = ws; u64* B64 = A64 + N; u64* R64 = B64 + N; u64* sub = R64 + 2*N;

    for (int32_t i = 0; i < N; i++) {
        A64[i] = (u64)(int64_t) a[i];
        B64[i] = (u64)(int64_t) b[i];
    }
    toom3_rec(R64, A64, B64, N, sub, base);
    for (int32_t i = 0; i < 2*N-1; i++)
        R[i] = (Torus32)(uint32_t) R64[i];

    for (int32_t i = 0; i < N-1; ++i)
        result[i] = R[i] - R[N+i];
    result[N-1] = R[N-1];

    delete[] R;
    delete[] ws;
}

// ---------------------------------------------------------------- reference

// naive product mod X^N+1 (as torusPolynomialMultNaive_aux), unsigned arithmetic to avoid UB
static void naiveMult(Torus32* result, const int32_t* a, const Torus32* b, const int32_t N) {
    for (int32_t i = 0; i < N; i++) {
        uint32_t ri = 0;
        for (int32_t j = 0; j <= i; j++) ri += (uint32_t)a[j] * (uint32_t)b[i-j];
        for (int32_t j = i+1; j < N; j++) ri -= (uint32_t)a[j] * (uint32_t)b[N+i-j];
        result[i] = (Torus32) ri;
    }
}

static double nowSec() {
    return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
    const int32_t N       = argc > 1 ? atoi(argv[1]) : 1024;
    const int32_t maxBase = argc > 2 ? atoi(argv[2]) : 200;
    const int32_t reps    = argc > 3 ? atoi(argv[3]) : 1000;
    const int32_t rounds  = 5;   // each cutoff is timed in several rounds; we keep the best (least noisy)

    vector<int32_t> a(N);
    vector<Torus32> b(N), expected(N), got(N);
    srand(12345);
    for (int32_t i = 0; i < N; i++) {
        a[i] = (rand() % (4*N)) - (2*N);   // same ranges as test-multiplication.cpp
        b[i] = rand();
    }
    naiveMult(expected.data(), a.data(), b.data(), N);

    printf("N=%d, cutoffs 1..%d, %d reps x %d rounds per cutoff\n\n", N, maxBase, reps, rounds);
    printf("%8s %14s %10s\n", "cutoff", "us/mult", "check");

    // sizes at which the recursion actually changes: cutoffs between two consecutive
    // sizes of the recursion tree (N, ceil(N/3), ...) give identical code paths.
    vector<double> us(maxBase+1, 0);
    int32_t bestBase = 1;
    double bestUs = 1e300;
    bool allOk = true;

    for (int32_t base = 1; base <= maxBase; base++) {
        toomMult(got.data(), a.data(), b.data(), N, base);
        bool ok = (got == expected);
        allOk = allOk && ok;

        double best = 1e300;
        for (int32_t r = 0; r < rounds; r++) {
            double t0 = nowSec();
            for (int32_t k = 0; k < reps; k++) {
                toomMult(got.data(), a.data(), b.data(), N, base);
                asm volatile("" : : "r"(got.data()) : "memory");  // keep the compiler from dropping the call
            }
            best = min(best, (nowSec() - t0) / reps);
        }
        us[base] = best * 1e6;
        if (us[base] < bestUs) { bestUs = us[base]; bestBase = base; }
        printf("%8d %14.3f %10s\n", base, us[base], ok ? "ok" : "MISMATCH");
        fflush(stdout);
    }

    // naive-only baseline = a cutoff >= N
    double t0 = nowSec();
    const int32_t nreps = max(1, reps/10);
    for (int32_t k = 0; k < nreps; k++) {
        toomMult(got.data(), a.data(), b.data(), N, N);
        asm volatile("" : : "r"(got.data()) : "memory");
    }
    double naiveUs = (nowSec() - t0) / nreps * 1e6;

    printf("\n=== summary ===\n");
    printf("all cutoffs correct:     %s\n", allOk ? "yes" : "NO");
    printf("best cutoff:             %d  (%.3f us/mult)\n", bestBase, bestUs);
    printf("current cutoff (12):     %.3f us/mult\n", maxBase >= 12 ? us[12] : -1.0);
    printf("pure naive (cutoff=N):   %.3f us/mult\n", naiveUs);
    printf("\nfive fastest cutoffs:\n");
    vector<int32_t> idx;
    for (int32_t i = 1; i <= maxBase; i++) idx.push_back(i);
    sort(idx.begin(), idx.end(), [&](int32_t x, int32_t y){ return us[x] < us[y]; });
    for (int i = 0; i < 5 && i < (int)idx.size(); i++)
        printf("  cutoff %3d: %.3f us\n", idx[i], us[idx[i]]);
    printf("\nIf several cutoffs are within noise of each other, prefer the one that is also\n"
           "good for other N (rerun with e.g. 512, 1024, 2048) and set TOOM_BASE to it.\n");
    return allOk ? 0 : 1;
}
