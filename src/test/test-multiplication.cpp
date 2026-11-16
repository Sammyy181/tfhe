#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstring>
#include "tfhe.h"
#include "polynomials.h"

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}


int32_t main(int32_t argc, char **argv) {
    const int32_t count = 1000; //number of tests to compare the 4 types of multiplication

    const int32_t N = 1024;
    IntPolynomial *a = new_IntPolynomial(N);
    TorusPolynomial *b = new_TorusPolynomial(N);
    TorusPolynomial *resNaive = new_TorusPolynomial(N);
    TorusPolynomial *resFFT = new_TorusPolynomial(N);
    TorusPolynomial *resKaratsuba = new_TorusPolynomial(N);
    TorusPolynomial *resNussbaumer = new_TorusPolynomial(N);
    //LagrangeHalfCPolynomial* test_fft = new_LagrangeHalfCPolynomial(N);

    double cycles_naive[count];
    double cycles_karatsuba[count];
    double cycles_fft[count];
    double cycles_nussbaumer[count];

    for (int32_t i = 0; i < count; ++i) {
        for (int32_t i = 0; i < N; i++) {
            a->coefs[i] = (rand() % (4 * N)) - (2 * N);
            b->coefsT[i] = rand();
        }

        //measure the execution time
        clock_t cstart, cend;

        cstart = clock();
        torusPolynomialMultNaive(resNaive, a, b);
        cend = clock();
        cycles_naive[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultKaratsuba(resKaratsuba, a, b);
        cend = clock();
        cycles_karatsuba[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultFFT(resFFT, a, b);
        //TorusPolynomial_ifft(test_fft,resNaive);
        //TorusPolynomial_fft(resFFT,test_fft);
        cend = clock();
        cycles_fft[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMulNussWS(resNussbaumer, a, b);
        cend = clock();
        cycles_nussbaumer[i] = cend - cstart;

        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resFFT->coefsT[i])) > 1) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resFFT->coefsT[i] << endl;
                dieDramatically("Naive != FFT\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resKaratsuba->coefsT[i])) > 1) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resKaratsuba->coefsT[i] << endl;
                dieDramatically("Naive != Karatsuba\n");
            }
        }
        //Nussbaumer + Winograd-Schoolbook is an exact computation: no tolerance
        for (int32_t i = 0; i < N; i++) {
            if (resNaive->coefsT[i] != resNussbaumer->coefsT[i]) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resNussbaumer->coefsT[i] << endl;
                dieDramatically("Naive != Nussbaumer\n");
            }
        }
    }

    // computing the average number of cycles per type of multiplication
    double cyc_naive = 0;
    double cyc_karatsuba = 0;
    double cyc_fft = 0;
    double cyc_nussbaumer = 0;
    for (int32_t i = 0; i < count; ++i) {
        cyc_naive += cycles_naive[i];
        cyc_karatsuba += cycles_karatsuba[i];
        cyc_fft += cycles_fft[i];
        cyc_nussbaumer += cycles_nussbaumer[i];
    }
    cyc_naive = cyc_naive / count;
    cyc_karatsuba = cyc_karatsuba / count;
    cyc_fft = cyc_fft / count;
    cyc_nussbaumer = cyc_nussbaumer / count;

    cout << "torusPolynomialMultNaive: " << cyc_naive << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultNaive time: " << (cyc_naive / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultKaratsuba: " << cyc_karatsuba << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultKaratsuba time: " << (cyc_karatsuba / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultFFT: " << cyc_fft << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultFFT time: " << (cyc_fft / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMulNussWS: " << cyc_nussbaumer << " clock cycles (average)" << endl;
    cout << "torusPolynomialMulNussWS time: " << (cyc_nussbaumer / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;


    // ******************************************************************************
    // Sum of products, as in one external product (CMux) of the bootstrapping:
    //     result = sum_{p < kpl} a_p * b_p
    // where the a_p are the kpl = (k+1)*l digit polynomials of the decomposed
    // accumulator and the b_p are polynomials of the bootstrapping key.
    // kpl and the range of the digits come from the default parameter set.
    //
    // The b_p are converted to the FFT / Nussbaumer domain beforehand (this is
    // not timed: in the bootstrapping, the key is converted once and for all).
    //  - Naive:      kpl full multiplications, summed
    //  - FFT:        kpl forward transforms of a_p, products accumulated in the
    //                FFT domain, 1 inverse transform
    //  - Nussbaumer: kpl forward transforms of a_p, Winograd-Schoolbook products
    //                accumulated in the Nussbaumer domain, 1 inverse transform
    //                (lazy interpolation)
    // ******************************************************************************
    TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(100);
    const int32_t kpl = params->tgsw_params->kpl;
    const int32_t halfBg = params->tgsw_params->halfBg;
    const int32_t countSum = 200; //number of tests for the sums of products

    IntPolynomial *as = new_IntPolynomial_array(kpl, N);
    TorusPolynomial *bs = new_TorusPolynomial_array(kpl, N);
    TorusPolynomial *tmp = new_TorusPolynomial(N);
    TorusPolynomial *sumNaive = new_TorusPolynomial(N);
    TorusPolynomial *sumFFT = new_TorusPolynomial(N);
    TorusPolynomial *sumNussbaumer = new_TorusPolynomial(N);
    LagrangeHalfCPolynomial *bsFFT = new_LagrangeHalfCPolynomial_array(kpl, N);
    LagrangeHalfCPolynomial *aFFT = new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomial *accumFFT = new_LagrangeHalfCPolynomial(N);
    NussWord *bsNussbaumer = new NussWord[kpl * NUSS_SIZE];
    NussWord *accumNussbaumer = new NussWord[NUSS_SIZE];

    double cyc_sum_naive = 0;
    double cyc_sum_fft = 0;
    double cyc_sum_nussbaumer = 0;
    int32_t maxErrFFT = 0;

    for (int32_t trial = 0; trial < countSum; ++trial) {
        for (int32_t p = 0; p < kpl; p++) {
            for (int32_t i = 0; i < N; i++) {
                as[p].coefs[i] = (rand() % (2 * halfBg)) - halfBg;
                bs[p].coefsT[i] = rand();
            }
            //conversion of the "key" polynomials (not timed)
            TorusPolynomial_ifft(bsFFT + p, bs + p);
            nussbaumerForward(bsNussbaumer + p * NUSS_SIZE, bs[p].coefsT);
        }

        clock_t cstart, cend;

        cstart = clock();
        torusPolynomialClear(sumNaive);
        for (int32_t p = 0; p < kpl; p++) {
            torusPolynomialMultNaive(tmp, as + p, bs + p);
            torusPolynomialAddTo(sumNaive, tmp);
        }
        cend = clock();
        cyc_sum_naive += cend - cstart;

        cstart = clock();
        LagrangeHalfCPolynomialClear(accumFFT);
        for (int32_t p = 0; p < kpl; p++) {
            IntPolynomial_ifft(aFFT, as + p);
            LagrangeHalfCPolynomialAddMul(accumFFT, aFFT, bsFFT + p);
        }
        TorusPolynomial_fft(sumFFT, accumFFT);
        cend = clock();
        cyc_sum_fft += cend - cstart;

        cstart = clock();
        memset(accumNussbaumer, 0, NUSS_SIZE * sizeof(NussWord));
        for (int32_t p = 0; p < kpl; p++) {
            torusPolynomialAddMulRNussbaumer(accumNussbaumer, as + p, bsNussbaumer + p * NUSS_SIZE, 1);
        }
        nussbaumerInverse(sumNussbaumer->coefsT, accumNussbaumer);
        cend = clock();
        cyc_sum_nussbaumer += cend - cstart;

        //FFT: approximate (floating point), we record the largest error
        for (int32_t i = 0; i < N; i++) {
            int32_t err = abs(int32_t(sumNaive->coefsT[i] - sumFFT->coefsT[i]));
            if (err > maxErrFFT) maxErrFFT = err;
            if (err > 1) {
                cerr << i << " " << sumNaive->coefsT[i] << " vs. " << sumFFT->coefsT[i] << endl;
                dieDramatically("sum: Naive != FFT\n");
            }
        }
        //Nussbaumer: exact
        for (int32_t i = 0; i < N; i++) {
            if (sumNaive->coefsT[i] != sumNussbaumer->coefsT[i]) {
                cerr << i << " " << sumNaive->coefsT[i] << " vs. " << sumNussbaumer->coefsT[i] << endl;
                dieDramatically("sum: Naive != Nussbaumer\n");
            }
        }
    }
    cyc_sum_naive = cyc_sum_naive / countSum;
    cyc_sum_fft = cyc_sum_fft / countSum;
    cyc_sum_nussbaumer = cyc_sum_nussbaumer / countSum;

    cout << "----- sum of " << kpl << " products (digits in [-" << halfBg << "," << halfBg << "]) -----" << endl;
    cout << "Naive sum: " << cyc_sum_naive << " clock cycles (average)" << endl;
    cout << "Naive sum time: " << (cyc_sum_naive / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "FFT sum (1 inverse transform): " << cyc_sum_fft << " clock cycles (average)" << endl;
    cout << "FFT sum time: " << (cyc_sum_fft / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << "FFT sum largest error: " << maxErrFFT << " (in units of 2^-32)" << endl;
    cout << endl;
    cout << "Nussbaumer sum (lazy interpolation, 1 inverse transform): " << cyc_sum_nussbaumer << " clock cycles (average)" << endl;
    cout << "Nussbaumer sum time: " << (cyc_sum_nussbaumer / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << "Nussbaumer sum largest error: 0 (exact)" << endl;
    cout << endl;

    delete[] accumNussbaumer;
    delete[] bsNussbaumer;
    delete_LagrangeHalfCPolynomial(accumFFT);
    delete_LagrangeHalfCPolynomial(aFFT);
    delete_LagrangeHalfCPolynomial_array(kpl, bsFFT);
    delete_TorusPolynomial(sumNussbaumer);
    delete_TorusPolynomial(sumFFT);
    delete_TorusPolynomial(sumNaive);
    delete_TorusPolynomial(tmp);
    delete_TorusPolynomial_array(kpl, bs);
    delete_IntPolynomial_array(kpl, as);
    delete_gate_bootstrapping_parameters(params);

    delete_IntPolynomial(a);
    delete_TorusPolynomial(b);
    delete_TorusPolynomial(resNaive);
    delete_TorusPolynomial(resKaratsuba);
    delete_TorusPolynomial(resNussbaumer);
    delete_TorusPolynomial(resFFT);
    return 0;
}
