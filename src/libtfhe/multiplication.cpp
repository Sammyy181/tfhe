#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "polynomials_arithmetic.h"
#include "nussbaumer.h"



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

	// Karatsuba recursivly
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

/* **********************************************************************
 * Nussbaumer(16, 8) + Winograd-Schoolbook multiplication
 * (see nussbaumer.h for the overview)
 * ******************************************************************** */

/**
 * r = a * b mod Z^8+1, with the Winograd-Schoolbook algorithm (Algorithm 1 of
 * Wang et al., "A High-Throughput Toom-Cook-4 Polynomial Multiplier ...").
 *
 * The schoolbook product is cut in F(2,3) blocks. One block takes a tile of
 * 4 coefficients of a and a filter of 3 coefficients of b, and produces
 *    c0 = a0.f0 + a1.f1 + a2.f2
 *    c1 = a1.f0 + a2.f1 + a3.f2
 * with 4 multiplications instead of 6:
 *    m0 = (a0 - a2).f0              m1 = (a1 + a2).(f0 + f1 + f2)
 *    m3 = (a1 - a3).f2              m2 = (a2 - a1).(f0 - f1 + f2)
 *    c0 = m0 + (m1 + m2)/2          c1 = (m1 - m2)/2 - m3
 * m1 + m2 and m1 - m2 are always even, so the two halvings are exact.
 *
 * For size 8: a is padded to A' = [0,0,a0..a7,0,0] (5 tiles, stride 2) and b
 * is reversed and padded to B' = [0,b7..b0] (3 filters, stride 3), which
 * gives 5 * 3 * 4 = 60 multiplications instead of 64.
 */
static inline void winogradSchoolbook8(NussWord* __restrict r, const NussWord* __restrict a, const NussWord* __restrict b) {
    // input transform of the 5 tiles (shared by the 3 filters)
    const NussWord Ap[12] = {0, 0, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], 0, 0};
    NussWord ta[5][4];
    for (int32_t q = 0; q < 5; q++) {
	const NussWord* T = Ap + 2*q;
	ta[q][0] = T[0] - T[2];
	ta[q][1] = T[1] + T[2];
	ta[q][2] = T[2] - T[1];
	ta[q][3] = T[1] - T[3];
    }

    // C = full product (15 coefficients, C[15] is always 0)
    NussWord C[16] = {0};
    for (int32_t i = 6; i >= 0; i -= 3) {
	// filter (B'[i], B'[i+1], B'[i+2]) = (b[8-i], b[7-i], b[6-i]), where b[8] = 0
	const NussWord f0 = (i == 0) ? 0 : b[8-i];
	const NussWord f1 = b[7-i];
	const NussWord f2 = b[6-i];
	// input transform of the filter
	const NussWord s = f0 + f2;
	const NussWord fp = s + f1;
	const NussWord fm = s - f1;
	for (int32_t q = 0; q < 5; q++) {
	    // element-wise multiplication
	    const NussWord m0 = ta[q][0] * f0;
	    const NussWord m1 = ta[q][1] * fp;
	    const NussWord m2 = ta[q][2] * fm;
	    const NussWord m3 = ta[q][3] * f2;
	    // output transform, accumulated at the position of this block
	    C[2*q + 6 - i] += m0 + ((m1 + m2) >> 1);
	    C[2*q + 7 - i] += ((m1 - m2) >> 1) - m3;
	}
    }

    // reduction mod Z^8+1
    for (int32_t n = 0; n < 7; n++) r[n] = C[n] - C[n+8];
    r[7] = C[7];
}

/**
 * result = a * b in the Nussbaumer domain:
 * 512 independent products of size 8 mod Z^8+1, with Winograd-Schoolbook.
 */
EXPORT void nussbaumerPointwiseMulWS(NussWord* result, const NussWord* a, const NussWord* b) {
    for (int32_t p = 0; p < NUSS_NB_PRODUCTS; p++)
	winogradSchoolbook8(result + p*NUSS_L2_SIZE, a + p*NUSS_L2_SIZE, b + p*NUSS_L2_SIZE);
}

/**
 * accum += poly1 * poly2 in the Nussbaumer domain (no inverse transform here:
 * the caller applies nussbaumerInverse once all the products are accumulated).
 * poly1 is transformed once and multiplied with ncols transformed polynomials.
 */
EXPORT void torusPolynomialAddMulRNussbaumer(NussWord* accum, const IntPolynomial* poly1, const NussWord* poly2Nuss, int32_t ncols) {
    assert(poly1->N == NUSS_N);
    NussWord poly1Nuss[NUSS_SIZE];
    NussWord prod[NUSS_SIZE];

    nussbaumerForward(poly1Nuss, poly1->coefs);
    for (int32_t c = 0; c < ncols; c++) {
	nussbaumerPointwiseMulWS(prod, poly1Nuss, poly2Nuss + c*NUSS_SIZE);
	NussWord* acc = accum + c*NUSS_SIZE;
	for (int32_t i = 0; i < NUSS_SIZE; i++) acc[i] += prod[i];
    }
}

// poly1, poly2 and result are polynomials mod X^N+1
// standalone version: forward transform of both, product, inverse transform
EXPORT void torusPolynomialMulNussWS(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2){
    assert(poly1->N == NUSS_N && poly2->N == NUSS_N && result->N == NUSS_N);
    NussWord poly1Nuss[NUSS_SIZE];
    NussWord poly2Nuss[NUSS_SIZE];
    NussWord prod[NUSS_SIZE];

    nussbaumerForward(poly1Nuss, poly1->coefs);
    nussbaumerForward(poly2Nuss, poly2->coefsT);
    nussbaumerPointwiseMulWS(prod, poly1Nuss, poly2Nuss);
    nussbaumerInverse(result->coefsT, prod);
}
