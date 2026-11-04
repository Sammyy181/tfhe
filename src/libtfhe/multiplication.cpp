#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"



using namespace std;



IntPolynomial::IntPolynomial(const int32_t N): N(N)
{
    this->coefs = new int32_t[N]; 
}

IntPolynomial::~IntPolynomial() {
    delete[] coefs;
}



TorusPolynomial::TorusPolynomial(const int32_t N): N(N)
{
    this->coefsT = new Torus32[N]; 
}

TorusPolynomial::~TorusPolynomial() {
    delete[] coefsT;
}


void torusPolynomialMultNaive_plain_aux(Torus32* __restrict result, const int32_t* __restrict poly1, const Torus32* __restrict poly2, const int32_t N) {
    const int32_t _2Nm1 = 2*N-1;
    Torus32 ri;
    for (int32_t i=0; i<N; i++) {
	ri=0;
	for (int32_t j=0; j<=i; j++) {
	    ri += poly1[j]*poly2[i-j];
	}
	result[i]=ri;
    }
    for (int32_t i=N; i<_2Nm1; i++) {
	ri=0;
	for (int32_t j=i-N+1; j<N; j++) {
	    ri += poly1[j]*poly2[i-j];
	}
	result[i]=ri;
    }
}

void torusPolynomialMultNaive_aux(Torus32* __restrict result, const int32_t* __restrict poly1, const Torus32* __restrict poly2, const int32_t N) {
    Torus32 ri;
    for (int32_t i=0; i<N; i++) {
		ri=0;
			for (int32_t j=0; j<=i; j++) {
		    	ri += poly1[j]*poly2[i-j];
			}
			for (int32_t j=i+1; j<N; j++) {
		    	ri -= poly1[j]*poly2[N+i-j];
			}
		result[i]=ri;
    }
}

/**
 * This is the naive external multiplication of an integer polynomial
 * with a torus polynomial. (this function should yield exactly the same
 * result as the karatsuba or fft version) 
 */
EXPORT void torusPolynomialMultNaive(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    assert(result!=poly2);
    assert(poly2->N==N && result->N==N);
    torusPolynomialMultNaive_aux(result->coefsT, poly1->coefs, poly2->coefsT, N);
}



/**
 * This function multiplies 2 polynomials (an integer poly and a torus poly) by using Karatsuba
 * The karatsuba function is torusPolynomialMultKaratsuba: it takes in input two polynomials and multiplies them 
 * To do that, it uses the auxiliary function Karatsuba_aux, which is recursive ad which works with 
 * the vectors containing the coefficients of the polynomials (primitive types)
 */

// A and B of size = size
// R of size = 2*size-1
EXPORT void Karatsuba_aux(Torus32* R, const int32_t* A, const Torus32* B, const int32_t size, const char* buf){
    const int32_t h = size / 2;
	const int32_t sm1 = size-1;

	//we stop the karatsuba recursion at h=4, because on my machine,
	//it seems to be optimal
	if (h<=4)
	{
	    torusPolynomialMultNaive_plain_aux(R, A, B, size);
	    return;
	}

	//we split the polynomials in 2
	int32_t* Atemp = (int32_t*) buf; buf += h*sizeof(int32_t);
	Torus32* Btemp = (Torus32*) buf; buf+= h*sizeof(Torus32);
	Torus32* Rtemp = (Torus32*) buf; buf+= size*sizeof(Torus32); 
	//Note: in the above line, I have put size instead of sm1 so that buf remains aligned on a power of 2

	for (int32_t i = 0; i < h; ++i)
	    Atemp[i] = A[i] + A[h+i];
	for (int32_t i = 0; i < h; ++i)
	    Btemp[i] = B[i] + B[h+i];

	// Karatsuba recursively
	Karatsuba_aux(R, A, B, h, buf); // (R[0],R[2*h-2]), (A[0],A[h-1]), (B[0],B[h-1])
	Karatsuba_aux(R+size, A+h, B+h, h, buf); // (R[2*h],R[4*h-2]), (A[h],A[2*h-1]), (B[h],B[2*h-1])
	Karatsuba_aux(Rtemp, Atemp, Btemp, h, buf);
	R[sm1]=0; //this one needs to be set manually
	for (int32_t i = 0; i < sm1; ++i) 
	    Rtemp[i] -= R[i] + R[size+i];
	for (int32_t i = 0; i < sm1; ++i) 
	    R[h+i] += Rtemp[i];

}

// poly1, poly2 and result are polynomials mod X^N+1
EXPORT void torusPolynomialMultKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
	const int32_t N = poly1->N;
	Torus32* R = new Torus32[2*N-1];
	char* buf = new char[16*N]; //that's large enough to store every tmp variables (2*2*N*4)
	
	// Karatsuba 
	Karatsuba_aux(R, poly1->coefs, poly2->coefsT, N, buf);

	// reduction mod X^N+1
	for (int32_t i = 0; i < N-1; ++i) 
	    result->coefsT[i] = R[i] - R[N+i];
	result->coefsT[N-1] = R[N-1];
	
	delete[] R;
	delete[] buf;
}

// poly1, poly2 and result are polynomials mod X^N+1
EXPORT void torusPolynomialAddMulRKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
	const int32_t N = poly1->N;
	Torus32* R = new Torus32[2*N-1];
	char* buf = new char[16*N]; //that's large enough to store every tmp variables (2*2*N*4)
	
	// Karatsuba 
	Karatsuba_aux(R, poly1->coefs, poly2->coefsT, N, buf);

	// reduction mod X^N+1
	for (int32_t i = 0; i < N-1; ++i) 
	    result->coefsT[i] += R[i] - R[N+i];
	result->coefsT[N-1] += R[N-1];
	
	delete[] R;
	delete[] buf;
}

// poly1, poly2 and result are polynomials mod X^N+1
EXPORT void torusPolynomialSubMulRKaratsuba(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
	const int32_t N = poly1->N;
	Torus32* R = new Torus32[2*N-1];
	char* buf = new char[16*N]; //that's large enough to store every tmp variables (2*2*N*4)
	
	// Karatsuba 
	Karatsuba_aux(R, poly1->coefs, poly2->coefsT, N, buf);

	// reduction mod X^N+1
	for (int32_t i = 0; i < N-1; ++i) 
	    result->coefsT[i] -= R[i] - R[N+i];
	result->coefsT[N-1] -= R[N-1];
	
	delete[] R;
	delete[] buf;
}

/**
 * Toom-Cook 3 multiplication (same contract as Karatsuba_aux: A and B of size = size,
 * R of size 2*size-1, plain product with no reduction mod X^N+1).
 *
 * Toom-3 splits A and B in 3 pieces (A = a0 + a1.Y + a2.Y^2, with Y=X^h) and evaluates
 * them at 0, 1, -1, -2, infinity, which gives 5 products of size h instead of 9.
 * The interpolation needs exact divisions by 2 and 3. Division by 3 is fine mod 2^k
 * (3 is invertible), but division by 2 is not: it destroys the top bit. So all the
 * recursion is done on 64-bit words; each division by 2 costs one correct top bit, and
 * we need only the lowest 32 bits at the end, so we have plenty of margin.
 *
 * The workspace "buf" must be at least ToomCook3_buf_size(size) bytes.
 */

typedef uint64_t u64;

//below this size, we use the naive multiplication
static const int32_t TOOM_BASE = 16;

//maximum number of products accumulated together by one (lazy) interpolation
static const int32_t TOOM_MAX_TAU = 16;

//R = sum_p A[p].B[p] (plain products of tau pairs of size n)
static void toom_naive64_sum(u64* R, const u64* const* A, const u64* const* B, const int32_t tau, const int32_t n) {
    for (int32_t i = 0; i < 2*n-1; i++) R[i] = 0;
    for (int32_t p = 0; p < tau; p++)
        for (int32_t i = 0; i < n; i++)
            for (int32_t j = 0; j < n; j++)
                R[i+j] += A[p][i]*B[p][j];
}

//workspace needed by toom3_rec (in u64 words)
static size_t toom3_ws(const int32_t size, const int32_t tau) {
    if (size <= TOOM_BASE) return 0;
    const int32_t h = (size+2)/3;
    return (12*(size_t)tau + 10)*h + toom3_ws(h, tau);
}

//workspace needed by ToomCook3_aux (in bytes)
static size_t ToomCook3_buf_size(const int32_t size) {
    return (4*(size_t)size + toom3_ws(size, 1)) * sizeof(u64);
}

//R = sum_{p<tau} A[p].B[p]. The pairs are evaluated one by one, the products at each point
//are accumulated (recursively: each of the 5 recursive calls is itself a sum of tau products),
//and the interpolation is done only once, on the sum. With tau=1, it is a plain product.
static void toom3_rec(u64* R, const u64* const* A, const u64* const* B, const int32_t tau, const int32_t size, u64* ws) {
    if (size <= TOOM_BASE) {
        toom_naive64_sum(R, A, B, tau, size);
        return;
    }
    const int32_t h = (size+2)/3;       // size of each piece (the last one is zero-padded)
    const int32_t rl = 2*h-1;           // size of a product of two pieces
    const u64 inv3 = 0xAAAAAAAAAAAAAAABULL; // 3^-1 mod 2^64

    //carve the workspace
    u64 *Ap[TOOM_MAX_TAU], *Bp[TOOM_MAX_TAU], *eA[TOOM_MAX_TAU], *eB[TOOM_MAX_TAU];
    for (int32_t p = 0; p < tau; p++) {
        Ap[p] = ws; ws += 3*h;          // zero-padded copies, so that the 3 pieces all have size h
        Bp[p] = ws; ws += 3*h;
        eA[p] = ws; ws += 3*h;          // evaluations at 1, -1, -2
        eB[p] = ws; ws += 3*h;
    }
    u64* r0  = ws; ws += 2*h;
    u64* r1  = ws; ws += 2*h;
    u64* rm  = ws; ws += 2*h;
    u64* r2  = ws; ws += 2*h;
    u64* rinf= ws; ws += 2*h;

    for (int32_t p = 0; p < tau; p++) {
        for (int32_t i = 0; i < 3*h; i++) {
            Ap[p][i] = (i < size) ? A[p][i] : 0;
            Bp[p][i] = (i < size) ? B[p][i] : 0;
        }
        const u64 *a0 = Ap[p], *a1 = Ap[p]+h, *a2 = Ap[p]+2*h;
        const u64 *b0 = Bp[p], *b1 = Bp[p]+h, *b2 = Bp[p]+2*h;
        u64 *eA1 = eA[p], *eAm = eA[p]+h, *eA2 = eA[p]+2*h;
        u64 *eB1 = eB[p], *eBm = eB[p]+h, *eB2 = eB[p]+2*h;

        //evaluations at 1, -1, -2 (p(-2) = 2.(p(-1)+a2) - a0, which reuses p(-1))
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
    }

    //5 recursive products at 0, 1, -1, -2, infinity (they run one after the other so they share the workspace)
    u64* rr[5] = { r0, r1, rm, r2, rinf };
    for (int32_t k = 0; k < 5; k++) {
        const u64 *pa[TOOM_MAX_TAU], *pb[TOOM_MAX_TAU];
        for (int32_t p = 0; p < tau; p++) {
            if (k == 0)      { pa[p] = Ap[p];            pb[p] = Bp[p]; }
            else if (k < 4)  { pa[p] = eA[p]+(k-1)*h;    pb[p] = eB[p]+(k-1)*h; }
            else             { pa[p] = Ap[p]+2*h;        pb[p] = Bp[p]+2*h; }
        }
        toom3_rec(rr[k], pa, pb, tau, h, ws);
    }

    //interpolation: the product is c0 + c1.Y + c2.Y^2 + c3.Y^3 + c4.Y^4, with
    //  r0=c0, rinf=c4, r1=c0+c1+c2+c3+c4, rm=c0-c1+c2-c3+c4, r2=c0-2c1+4c2-8c3+16c4
    //c0 and c4 are already known. We overwrite r1, rm and r2 by c2, c1, c3.
    for (int32_t i = 0; i < rl; i++) {
        u64 c0 = r0[i], c4 = rinf[i];
        u64 c2 = ((r1[i] + rm[i]) >> 1) - c0 - c4;     // exact /2
        u64 t  = (r1[i] - rm[i]) >> 1;                 // c1 + c3
        u64 u  = (c0 + (c2 << 2) + (c4 << 4) - r2[i]) >> 1; // c1 + 4.c3
        u64 c3 = (u - t) * inv3;                       // exact /3
        u64 c1 = t - c3;
        r1[i] = c2; rm[i] = c1; r2[i] = c3;
    }

    //recomposition: R = sum_k c_k . X^(k.h) (terms beyond 2*size-2 are zero padding)
    const int32_t Rlen = 2*size-1;
    for (int32_t i = 0; i < Rlen; i++) R[i] = 0;
    const u64* c[5] = { r0, rm, r1, r2, rinf };
    for (int32_t k = 0; k < 5; k++)
        for (int32_t i = 0; i < rl && k*h+i < Rlen; i++)
            R[k*h+i] += c[k][i];
}

// A and B of size = size
// R of size = 2*size-1
// buf of size >= ToomCook3_buf_size(size) bytes
EXPORT void ToomCook3_aux(Torus32* R, const int32_t* A, const Torus32* B, const int32_t size, const char* buf){
    u64* ws = (u64*) buf;
    u64* A64 = ws; ws += size;
    u64* B64 = ws; ws += size;
    u64* R64 = ws; ws += 2*size;

    for (int32_t i = 0; i < size; i++) {
        A64[i] = (u64)(int64_t) A[i];
        B64[i] = (u64)(int64_t) B[i];
    }
    const u64 *pa[1] = { A64 }, *pb[1] = { B64 };
    toom3_rec(R64, pa, pb, 1, size, ws);
    for (int32_t i = 0; i < 2*size-1; i++)
        R[i] = (Torus32) (uint32_t) R64[i];  // keep the 32 low bits: this is the torus
}

// poly1, poly2 and result are polynomials mod X^N+1
EXPORT void torusPolynomialMultToomCook3(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
	const int32_t N = poly1->N;
	Torus32* R = new Torus32[2*N - 1];
	char* buf = new char[ToomCook3_buf_size(N)];

	//Toom-Cook
	ToomCook3_aux(R, poly1->coefs, poly2->coefsT, N, buf);

	// reduction mod X^N + 1
	for (int32_t i = 0; i < N-1; ++i)
		result->coefsT[i] = R[i] - R[N+i];
	result->coefsT[N-1] = R[N-1];

	delete[] R;
	delete[] buf;
}

/**
 * Toom-Cook 4 and Toom-Cook 5 multiplication (same contract as ToomCook3_aux).
 *
 * Toom-K splits A and B in K pieces (A = sum a_i.Y^i, with Y=X^h) and evaluates them at
 * 2K-1 points, which gives 2K-1 products of size h instead of K^2 (7 for K=4, 9 for K=5).
 * Toom-5 uses the finite points 0, 1, -1, 2, -2, 3, ... plus infinity. Toom-4 uses the points
 * 0, 1/2, 1, -1, 2, -2 and infinity, with Bodrato's hand-written interpolation sequence
 * (toom4_interp).
 *
 * For Toom-5 the interpolation matrix (the inverse of the Vandermonde matrix) is computed once, with
 * exact rational arithmetic, the first time it is needed. Row i of the inverse is written
 * as (integer row) / (2^s . odd). The odd part is inverted mod 2^64, and the 2^s is a
 * shift: as in Toom-3, all the recursion is done on 64-bit words, and each level of
 * recursion costs s correct top bits. Since only the 32 low bits are needed at the end,
 * the recursion is allowed to lose 32 bits in total (depth . s <= 32). The wrappers check
 * it, and fall back on Karatsuba for the (very large N) cases where it would not hold.
 */

typedef __int128 i128;

struct Frac { i128 n, d; }; // d > 0, always reduced

static i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i128 t = a % b; a = b; b = t; }
    return a;
}

static Frac fracMake(i128 n, i128 d) {
    if (d < 0) { n = -n; d = -d; }
    i128 g = gcd128(n, d);
    if (g > 1) { n /= g; d /= g; }
    Frac f = { n, d };
    return f;
}
static Frac fracSub(Frac a, Frac b) { return fracMake(a.n*b.d - b.n*a.d, a.d*b.d); }
static Frac fracMul(Frac a, Frac b) { return fracMake(a.n*b.n, a.d*b.d); }
static Frac fracDiv(Frac a, Frac b) { return fracMake(a.n*b.d, a.d*b.n); }

struct ToomInterp {
    int32_t n;          // number of evaluation points = 2K-1 (the last one is infinity)
    int64_t pts[9];     // the n-1 finite points
    u64 M[9][9];        // integer interpolation matrix (mod 2^64)
    u64 oddinv[9];      // inverse mod 2^64 of the odd part of the denominator of each row
    int32_t shift[9];   // power of 2 in the denominator of each row
    int32_t maxShift;
};

static ToomInterp toomBuildInterp(const int32_t K) {
    ToomInterp T;
    const int32_t n = 2*K-1;
    T.n = n;
    for (int32_t j = 0; j < n-1; j++)
        T.pts[j] = (j == 0) ? 0 : ((j % 2) ? (j+1)/2 : -(j/2));  // 0, 1, -1, 2, -2, 3, ...

    //Vandermonde matrix V (r = V.c) next to the identity, then Gauss-Jordan
    Frac a[9][18];
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < 2*n; j++)
            a[i][j] = fracMake(0, 1);
    for (int32_t i = 0; i < n-1; i++) {
        i128 p = 1;
        for (int32_t j = 0; j < n; j++) { a[i][j] = fracMake(p, 1); p *= T.pts[i]; }
    }
    a[n-1][n-1] = fracMake(1, 1);   // infinity: the leading coefficient
    for (int32_t i = 0; i < n; i++) a[i][n+i] = fracMake(1, 1);

    for (int32_t c = 0; c < n; c++) {
        int32_t piv = c;
        while (a[piv][c].n == 0) piv++;
        for (int32_t j = 0; j < 2*n; j++) { Frac t = a[c][j]; a[c][j] = a[piv][j]; a[piv][j] = t; }
        Frac pv = a[c][c];
        for (int32_t j = 0; j < 2*n; j++) a[c][j] = fracDiv(a[c][j], pv);
        for (int32_t i = 0; i < n; i++) {
            if (i == c || a[i][c].n == 0) continue;
            Frac f = a[i][c];
            for (int32_t j = 0; j < 2*n; j++) a[i][j] = fracSub(a[i][j], fracMul(f, a[c][j]));
        }
    }

    //row i of the inverse is a[i][n..2n-1]: put it on a common denominator
    T.maxShift = 0;
    for (int32_t i = 0; i < n; i++) {
        i128 den = 1;
        for (int32_t j = 0; j < n; j++) {
            const i128 d = a[i][n+j].d;
            den = den / gcd128(den, d) * d;
        }
        for (int32_t j = 0; j < n; j++)
            T.M[i][j] = (u64)(int64_t)(a[i][n+j].n * (den / a[i][n+j].d));
        int32_t s = 0;
        while (((den >> s) & 1) == 0) s++;
        const u64 odd = (u64)(den >> s);
        u64 x = odd;                                        // Newton: x = odd^-1 mod 2^64
        for (int32_t it = 0; it < 6; it++) x *= 2 - odd*x;
        T.oddinv[i] = x;
        T.shift[i] = s;
        if (s > T.maxShift) T.maxShift = s;
    }
    return T;
}

template<int K>
static const ToomInterp& toomInterp() {
    static const ToomInterp T = toomBuildInterp(K);
    return T;
}

//workspace needed by toomk_rec (in u64 words)
template<int K>
static size_t toomk_ws(const int32_t size, const int32_t base, const int32_t tau) {
    if (size <= base) return 0;
    const int32_t h = (size+K-1)/K;
    return ((size_t)tau*(2*K + 2*(2*K-1)) + 2*(2*K-1))*h + toomk_ws<K>(h, base, tau);
}

//number of correct top bits lost by one level of recursion
//(Toom-4 uses the hand-written interpolation sequence below; see the tracking there)
template<int K>
static int32_t toomk_loss() {
    return (K == 4) ? 4 : toomInterp<K>().maxShift;
}

//true if the recursion cannot lose more than 32 correct bits
template<int K>
static bool toomk_safe(int32_t size, const int32_t base) {
    int32_t depth = 0;
    while (size > base) { size = (size+K-1)/K; depth++; }
    return depth*toomk_loss<K>() <= 32;
}

//exact division by (odd . 2^s), given inv = odd^-1 mod 2^64: costs s correct top bits
static inline u64 exdiv(const u64 x, const u64 inv, const int32_t s) { return (x*inv) >> s; }
static_assert(3ULL*0xAAAAAAAAAAAAAAABULL == 1, "inverse of 3");
static_assert(9ULL*0x8E38E38E38E38E39ULL == 1, "inverse of 9");
static_assert(15ULL*0xEEEEEEEEEEEEEEEFULL == 1, "inverse of 15");

//Toom-4 evaluation of one coefficient index t of the pieces a0..a3, at the 7 points
//0, 1/2 (scaled by 8: 8.a0+4.a1+2.a2+a3 = 8.p(1/2)), 1, -1, 2, -2, infinity
static inline void toom4_eval(u64** e, const int32_t t, const u64 a0, const u64 a1, const u64 a2, const u64 a3) {
    const u64 s02 = a0 + a2, s13 = a1 + a3;
    const u64 ev = a0 + (a2 << 2), od = (a1 << 1) + (a3 << 3);
    e[0][t] = a0;
    e[1][t] = (a0 << 3) + (a1 << 2) + (a2 << 1) + a3;
    e[2][t] = s02 + s13;
    e[3][t] = s02 - s13;
    e[4][t] = ev + od;
    e[5][t] = ev - od;
    e[6][t] = a3;
}

//Toom-4 interpolation (Bodrato's sequence). Input: the 7 products at 0, 1/2 (scaled by 64),
//1, -1, 2, -2, infinity; output: the 7 coefficients c0..c6 (in place).
//Number of correct top bits lost (tracked through the sequence): 1 for the /2, then
//min(64,63)-3 = 60 for the /24, 62 for the /18, 60 for the /60, i.e. 4 bits in total.
static inline void toom4_interp(u64& r0, u64& r1, u64& r2, u64& r3, u64& r4, u64& r5, u64& r6) {
    const u64 inv3 = 0xAAAAAAAAAAAAAAABULL, inv9 = 0x8E38E38E38E38E39ULL, inv15 = 0xEEEEEEEEEEEEEEEFULL;
    r1 = r1 + r4;
    r5 = r5 - r4;
    r3 = (r3 - r2) >> 1;                 // exact /2
    r4 = r4 - r0;
    r4 = r4 - (r6 << 6);
    r4 = (r4 << 1) + r5;
    r2 = r2 + r3;
    r1 = r1 - 65*r2;
    r2 = r2 - r6;
    r2 = r2 - r0;
    r1 = r1 + 45*r2;
    r4 = exdiv(r4 - (r2 << 3), inv3, 3); // exact /24
    r5 = r5 + r1;
    r1 = exdiv(r1 + (r3 << 4), inv9, 1); // exact /18
    r3 = -(r3 + r1);
    r5 = exdiv((r1*30) - r5, inv15, 2);  // exact /60
    r2 = r2 - r4;
    r1 = r1 - r5;
}

template<int K>
static void toomk_rec(u64* R, const u64* const* A, const u64* const* B, const int32_t tau, const int32_t size, u64* ws, const int32_t base, const ToomInterp& T) {
    if (size <= base) {
        toom_naive64_sum(R, A, B, tau, size);
        return;
    }
    const int32_t n = 2*K-1;
    const int32_t h = (size+K-1)/K;     // size of each piece (the last one is zero-padded)
    const int32_t rl = 2*h-1;           // size of a product of two pieces

    //carve the workspace
    u64 *Ap[TOOM_MAX_TAU], *Bp[TOOM_MAX_TAU], *eA[TOOM_MAX_TAU][9], *eB[TOOM_MAX_TAU][9], *r[9];
    for (int32_t p = 0; p < tau; p++) {
        Ap[p] = ws; ws += K*h;
        Bp[p] = ws; ws += K*h;
        for (int32_t j = 0; j < n; j++) { eA[p][j] = ws; ws += h; }
        for (int32_t j = 0; j < n; j++) { eB[p][j] = ws; ws += h; }
    }
    for (int32_t j = 0; j < n; j++) { r[j] = ws; ws += 2*h; }

    for (int32_t p = 0; p < tau; p++) {
        const u64 *a = Ap[p], *b = Bp[p];
        for (int32_t i = 0; i < K*h; i++) {
            Ap[p][i] = (i < size) ? A[p][i] : 0;
            Bp[p][i] = (i < size) ? B[p][i] : 0;
        }

        //evaluations
        if (K == 4) {
            for (int32_t t = 0; t < h; t++) {
                toom4_eval(eA[p], t, a[t], a[h+t], a[2*h+t], a[3*h+t]);
                toom4_eval(eB[p], t, b[t], b[h+t], b[2*h+t], b[3*h+t]);
            }
        } else
        for (int32_t j = 0; j < n; j++) {   //Horner
            if (j == n-1) {                                  // infinity
                for (int32_t t = 0; t < h; t++) { eA[p][j][t] = a[(K-1)*h+t]; eB[p][j][t] = b[(K-1)*h+t]; }
            } else if (T.pts[j] == 0) {
                for (int32_t t = 0; t < h; t++) { eA[p][j][t] = a[t]; eB[p][j][t] = b[t]; }
            } else {
                const u64 x = (u64) T.pts[j];
                for (int32_t t = 0; t < h; t++) {
                    u64 va = a[(K-1)*h+t], vb = b[(K-1)*h+t];
                    for (int32_t i = K-2; i >= 0; i--) {
                        va = va*x + a[i*h+t];
                        vb = vb*x + b[i*h+t];
                    }
                    eA[p][j][t] = va; eB[p][j][t] = vb;
                }
            }
        }
    }

    //2K-1 recursive products, each one is the sum over the tau pairs
    //(one after the other so they share the workspace)
    for (int32_t j = 0; j < n; j++) {
        const u64 *pa[TOOM_MAX_TAU], *pb[TOOM_MAX_TAU];
        for (int32_t p = 0; p < tau; p++) { pa[p] = eA[p][j]; pb[p] = eB[p][j]; }
        toomk_rec<K>(r[j], pa, pb, tau, h, ws, base, T);
    }

    //interpolation: r[i] is overwritten by the i-th coefficient c_i
    if (K == 4) {
        for (int32_t t = 0; t < rl; t++)
            toom4_interp(r[0][t], r[1][t], r[2][t], r[3][t], r[4][t], r[5][t], r[6][t]);
    } else
    for (int32_t t = 0; t < rl; t++) {
        u64 v[9];
        for (int32_t j = 0; j < n; j++) v[j] = r[j][t];
        for (int32_t i = 0; i < n; i++) {
            u64 acc = 0;
            for (int32_t j = 0; j < n; j++) acc += T.M[i][j]*v[j];
            r[i][t] = (acc*T.oddinv[i]) >> T.shift[i];   // exact division by the denominator
        }
    }

    //recomposition: R = sum_k c_k . X^(k.h) (terms beyond 2*size-2 are zero padding)
    const int32_t Rlen = 2*size-1;
    for (int32_t i = 0; i < Rlen; i++) R[i] = 0;
    for (int32_t k = 0; k < n; k++)
        for (int32_t i = 0; i < rl && k*h+i < Rlen; i++)
            R[k*h+i] += r[k][i];
}

//below this size, we use the naive multiplication (to be tuned on each machine)
static const int32_t TOOM4_BASE = 16;
static const int32_t TOOM5_BASE = 20;

template<int K>
static size_t toomk_buf_size(const int32_t size, const int32_t base) {
    return (4*(size_t)size + toomk_ws<K>(size, base, 1)) * sizeof(u64);
}

template<int K>
static void toomk_aux(Torus32* R, const int32_t* A, const Torus32* B, const int32_t size, const char* buf, const int32_t base) {
    u64* ws = (u64*) buf;
    u64* A64 = ws; ws += size;
    u64* B64 = ws; ws += size;
    u64* R64 = ws; ws += 2*size;

    for (int32_t i = 0; i < size; i++) {
        A64[i] = (u64)(int64_t) A[i];
        B64[i] = (u64)(int64_t) B[i];
    }
    const u64 *pa[1] = { A64 }, *pb[1] = { B64 };
    toomk_rec<K>(R64, pa, pb, 1, size, ws, base, toomInterp<K>());
    for (int32_t i = 0; i < 2*size-1; i++)
        R[i] = (Torus32) (uint32_t) R64[i];  // keep the 32 low bits: this is the torus
}

template<int K>
static void toomk_mult(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2, const int32_t base) {
    const int32_t N = poly1->N;
    if (!toomk_safe<K>(N, base)) { // too many levels for the 64-bit words
        torusPolynomialMultKaratsuba(result, poly1, poly2);
        return;
    }
    Torus32* R = new Torus32[2*N - 1];
    char* buf = new char[toomk_buf_size<K>(N, base)];

    toomk_aux<K>(R, poly1->coefs, poly2->coefsT, N, buf, base);

    // reduction mod X^N + 1
    for (int32_t i = 0; i < N-1; ++i)
        result->coefsT[i] = R[i] - R[N+i];
    result->coefsT[N-1] = R[N-1];

    delete[] R;
    delete[] buf;
}

// A and B of size = size
// R of size = 2*size-1
// buf of size >= toomk_buf_size<4 or 5>(size, TOOM4_BASE or TOOM5_BASE) bytes
EXPORT void ToomCook4_aux(Torus32* R, const int32_t* A, const Torus32* B, const int32_t size, const char* buf){
    toomk_aux<4>(R, A, B, size, buf, TOOM4_BASE);
}
EXPORT void ToomCook5_aux(Torus32* R, const int32_t* A, const Torus32* B, const int32_t size, const char* buf){
    toomk_aux<5>(R, A, B, size, buf, TOOM5_BASE);
}

// poly1, poly2 and result are polynomials mod X^N+1
EXPORT void torusPolynomialMultToomCook4(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
    toomk_mult<4>(result, poly1, poly2, TOOM4_BASE);
}
EXPORT void torusPolynomialMultToomCook5(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
    toomk_mult<5>(result, poly1, poly2, TOOM5_BASE);
}

/**
 * Sum of several products, with optional lazy interpolation (Toom-3 and Toom-4).
 *
 * result = sum_{j<tau} poly1[j] * poly2[j]   (mod X^N+1, result is overwritten)
 *
 * lazy=false: the tau products are computed one by one (evaluation, pointwise products and
 * interpolation each time) and added.
 * lazy=true: each pair is evaluated, the pointwise products of all the pairs are added
 * in the evaluation space, and the interpolation runs only once on the sum. The same
 * is done at every level of the recursion (the sub-products are themselves sums), so the
 * interpolation cost does not depend on tau. Since evaluation and interpolation are linear,
 * and the 64-bit arithmetic is exact mod 2^64, the result is the same as with lazy=false.
 * The pairs are processed by groups of at most TOOM_MAX_TAU.
 */
static size_t toom_sum_words(const int32_t K, const int32_t N, const int32_t g) {
    return (size_t)(2*g+2)*N + ((K == 3) ? toom3_ws(N, g) : toomk_ws<4>(N, TOOM4_BASE, g));
}

static void toom_sum_impl(const int32_t K, TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy) {
    const int32_t N = result->N;
    for (int32_t i = 0; i < N; i++) result->coefsT[i] = 0;

    if (K == 4 && !toomk_safe<4>(N, TOOM4_BASE)) { // too many levels for the 64-bit words
        for (int32_t j = 0; j < tau; j++)
            torusPolynomialAddMulRKaratsuba(result, poly1[j], poly2[j]);
        return;
    }

    const int32_t group = lazy ? TOOM_MAX_TAU : 1;
    Torus32* R = new Torus32[2*N - 1];
    u64* buf = new u64[toom_sum_words(K, N, group)];

    for (int32_t start = 0; start < tau; start += group) {
        const int32_t g = (tau-start < group) ? tau-start : group;
        u64* ws = buf;
        const u64 *pa[TOOM_MAX_TAU], *pb[TOOM_MAX_TAU];
        for (int32_t p = 0; p < g; p++) {
            u64* A64 = ws; ws += N;
            u64* B64 = ws; ws += N;
            for (int32_t i = 0; i < N; i++) {
                A64[i] = (u64)(int64_t) poly1[start+p]->coefs[i];
                B64[i] = (u64)(int64_t) poly2[start+p]->coefsT[i];
            }
            pa[p] = A64; pb[p] = B64;
        }
        u64* R64 = ws; ws += 2*N;

        if (K == 3) toom3_rec(R64, pa, pb, g, N, ws);
        else        toomk_rec<4>(R64, pa, pb, g, N, ws, TOOM4_BASE, toomInterp<4>());

        for (int32_t i = 0; i < 2*N-1; i++)
            R[i] = (Torus32) (uint32_t) R64[i];  // keep the 32 low bits: this is the torus

        // reduction mod X^N + 1
        for (int32_t i = 0; i < N-1; ++i)
            result->coefsT[i] += R[i] - R[N+i];
        result->coefsT[N-1] += R[N-1];
    }

    delete[] R;
    delete[] buf;
}

EXPORT void torusPolynomialMultSumToomCook3(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy){
    toom_sum_impl(3, result, poly1, poly2, tau, lazy);
}
EXPORT void torusPolynomialMultSumToomCook4(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy){
    toom_sum_impl(4, result, poly1, poly2, tau, lazy);
}

/**
 * Sum of several products with the FFT, with optional lazy interpolation.
 *
 * result = sum_{j<tau} poly1[j] * poly2[j]   (mod X^N+1, result is overwritten)
 *
 * lazy=false: each product goes through the 2 forward FFTs, the pointwise product and the
 * inverse FFT, and the tau results are added (this is what torusPolynomialAddMulRFFT does).
 * lazy=true: the pointwise products are accumulated in the Lagrange (evaluation) space, and
 * the inverse FFT, which is the interpolation, runs only once on the sum. Compared with
 * lazy=false it saves tau-1 inverse FFTs. Like the FFT itself, the result is approximate
 * (the rounding errors of the tau products add up, but the rounding happens only once).
 */
EXPORT void torusPolynomialMultSumFFT(TorusPolynomial* result, const IntPolynomial* const* poly1, const TorusPolynomial* const* poly2, const int32_t tau, const bool lazy){
    const int32_t N = result->N;
    if (!lazy) {
        for (int32_t i = 0; i < N; i++) result->coefsT[i] = 0;
        for (int32_t j = 0; j < tau; j++)
            torusPolynomialAddMulRFFT(result, poly1[j], poly2[j]);
        return;
    }
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3, N);
    LagrangeHalfCPolynomialClear(tmp+2);
    for (int32_t j = 0; j < tau; j++) {
        IntPolynomial_ifft(tmp+0, poly1[j]);
        TorusPolynomial_ifft(tmp+1, poly2[j]);
        LagrangeHalfCPolynomialAddMul(tmp+2, tmp+0, tmp+1);   // accumulate in the evaluation space
    }
    TorusPolynomial_fft(result, tmp+2);                       // one single interpolation
    delete_LagrangeHalfCPolynomial_array(3, tmp);
}
