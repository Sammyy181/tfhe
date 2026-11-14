#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
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
    const int32_t count = 1000; //number of tests to compare the 3 types of multiplication

    const int32_t N = 1024;
    IntPolynomial *a = new_IntPolynomial(N);
    TorusPolynomial *b = new_TorusPolynomial(N);
    TorusPolynomial *resNaive = new_TorusPolynomial(N);
    TorusPolynomial *resFFT = new_TorusPolynomial(N);
    TorusPolynomial *resKaratsuba = new_TorusPolynomial(N);
    TorusPolynomial *resToomCook3 = new_TorusPolynomial(N);
    TorusPolynomial *resToomCook4 = new_TorusPolynomial(N);
    TorusPolynomial *resToomCook5 = new_TorusPolynomial(N);
    TorusPolynomial *resToomCook4S = new_TorusPolynomial(N);
    TorusPolynomial *resNTT = new_TorusPolynomial(N);
    //LagrangeHalfCPolynomial* test_fft = new_LagrangeHalfCPolynomial(N);

    double cycles_naive[count];
    double cycles_karatsuba[count];
    double cycles_fft[count];
    double cycles_toomcook3[count];
    double cycles_toomcook4[count];
    double cycles_toomcook5[count];
    double cycles_toomcook4s[count];
    double cycles_ntt[count];

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
        torusPolynomialMultToomCook3(resToomCook3, a, b);
        cend = clock();
        cycles_toomcook3[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultToomCook4(resToomCook4, a, b);
        cend = clock();
        cycles_toomcook4[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultToomCook5(resToomCook5, a, b);
        cend = clock();
        cycles_toomcook5[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultToomCook4Strided(resToomCook4S, a, b);
        cend = clock();
        cycles_toomcook4s[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultFFT(resFFT, a, b);
        //TorusPolynomial_ifft(test_fft,resNaive);
        //TorusPolynomial_fft(resFFT,test_fft);
        cend = clock();
        cycles_fft[i] = cend - cstart;

        cstart = clock();
        torusPolynomialMultNTT(resNTT, a, b);
        cend = clock();
        cycles_ntt[i] = cend - cstart;

        for (int32_t i = 0; i < N; i++) {
            if (resNaive->coefsT[i] != resNTT->coefsT[i]) {   // the NTT is exact
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resNTT->coefsT[i] << endl;
                dieDramatically("Naive != NTT\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resFFT->coefsT[i])) > 1) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resFFT->coefsT[i] << endl;
                dieDramatically("Naive != FFT\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resKaratsuba->coefsT[i])) > 0) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resKaratsuba->coefsT[i] << endl;
                dieDramatically("Naive != Karatsuba\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resToomCook3->coefsT[i])) > 0) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resToomCook3->coefsT[i] << endl;
                dieDramatically("Naive != Toom Cook-3\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resToomCook4->coefsT[i])) > 0) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resToomCook4->coefsT[i] << endl;
                dieDramatically("Naive != Toom Cook-4\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (abs(int32_t(resNaive->coefsT[i] - resToomCook5->coefsT[i])) > 0) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resToomCook5->coefsT[i] << endl;
                dieDramatically("Naive != Toom Cook-5\n");
            }
        }
        for (int32_t i = 0; i < N; i++) {
            if (resNaive->coefsT[i] != resToomCook4S->coefsT[i]) {
                cerr << i << " " << resNaive->coefsT[i] << " vs. " << resToomCook4S->coefsT[i] << endl;
                dieDramatically("Naive != striding Toom Cook-4\n");
            }
        }
    }

    // computing the average number of cycles per type of multiplication
    double cyc_naive = 0;
    double cyc_karatsuba = 0;
    double cyc_fft = 0;
    double cyc_toomcook3 = 0;
    double cyc_toomcook4 = 0;
    double cyc_toomcook5 = 0;
    double cyc_toomcook4s = 0;
    double cyc_ntt = 0;
    for (int32_t i = 0; i < count; ++i) {
        cyc_naive += cycles_naive[i];
        cyc_karatsuba += cycles_karatsuba[i];
        cyc_fft += cycles_fft[i];
        cyc_toomcook3 += cycles_toomcook3[i];
        cyc_toomcook4 += cycles_toomcook4[i];
        cyc_toomcook5 += cycles_toomcook5[i];
        cyc_toomcook4s += cycles_toomcook4s[i];
        cyc_ntt += cycles_ntt[i];
    }
    cyc_naive = cyc_naive / count;
    cyc_karatsuba = cyc_karatsuba / count;
    cyc_fft = cyc_fft / count;
    cyc_toomcook3 = cyc_toomcook3 / count;
    cyc_toomcook4 = cyc_toomcook4 / count;
    cyc_toomcook5 = cyc_toomcook5 / count;
    cyc_toomcook4s = cyc_toomcook4s / count;
    cyc_ntt = cyc_ntt / count;

    cout << "torusPolynomialMultNaive: " << cyc_naive << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultNaive time: " << (cyc_naive / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultKaratsuba: " << cyc_karatsuba << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultKaratsuba time: " << (cyc_karatsuba / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultFFT: " << cyc_fft << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultFFT time: " << (cyc_fft / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultNTT: " << cyc_ntt << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultNTT time: " << (cyc_ntt / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultToomCook3: "<< cyc_toomcook3 << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultToomCook3 time: " << (cyc_toomcook3 / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultToomCook4: "<< cyc_toomcook4 << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultToomCook4 time: " << (cyc_toomcook4 / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultToomCook5: "<< cyc_toomcook5 << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultToomCook5 time: " << (cyc_toomcook5 / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultToomCook4Strided: "<< cyc_toomcook4s << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultToomCook4Strided time: " << (cyc_toomcook4s / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;

    const int32_t tau = 20;        
    const int32_t countSum = 200;
    IntPolynomial *as[tau];
    TorusPolynomial *bs[tau];
    for (int32_t j = 0; j < tau; j++) {
        as[j] = new_IntPolynomial(N);
        bs[j] = new_TorusPolynomial(N);
    }
    TorusPolynomial *tmp = new_TorusPolynomial(N);
    TorusPolynomial *sumNaive = new_TorusPolynomial(N);
    TorusPolynomial *sumRes = new_TorusPolynomial(N);
    const int32_t nSum = 9;
    double cyc_sum[nSum] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    const char* sumNames[nSum] = {"Naive", "FFT", "FFT lazy", "Toom-3 non-lazy", "Toom-3 lazy", "Toom-4 non-lazy", "Toom-4 lazy", "NTT", "NTT lazy"};
    int32_t maxErrFFT = 0, maxErrFFTLazy = 0;
    for (int32_t it = 0; it < countSum; ++it) {
        for (int32_t j = 0; j < tau; j++)
            for (int32_t i = 0; i < N; i++) {
                as[j]->coefs[i] = (rand() % (4 * N)) - (2 * N);
                bs[j]->coefsT[i] = rand();
            }
        for (int32_t i = 0; i < N; i++) sumNaive->coefsT[i] = 0;
        for (int32_t j = 0; j < tau; j++) {
            torusPolynomialMultNaive(tmp, as[j], bs[j]);
            for (int32_t i = 0; i < N; i++) sumNaive->coefsT[i] += tmp->coefsT[i];
        }
        for (int32_t m = 0; m < nSum; m++) {
            clock_t cstart = clock();
            if (m < 2) {   // Naive and FFT: product by product, then added
                for (int32_t i = 0; i < N; i++) sumRes->coefsT[i] = 0;
                for (int32_t j = 0; j < tau; j++) {
                    if (m == 0) torusPolynomialMultNaive(tmp, as[j], bs[j]);
                    else        torusPolynomialMultFFT(tmp, as[j], bs[j]);
                    for (int32_t i = 0; i < N; i++) sumRes->coefsT[i] += tmp->coefsT[i];
                }
            }
            else if (m == 2) torusPolynomialMultSumFFT(sumRes, as, bs, tau, true);
            else if (m < 5)  torusPolynomialMultSumToomCook3(sumRes, as, bs, tau, m == 4);
            else if (m < 7)  torusPolynomialMultSumToomCook4(sumRes, as, bs, tau, m == 6);
            else             torusPolynomialMultSumNTT(sumRes, as, bs, tau, m == 8);
            clock_t cend = clock();
            cyc_sum[m] += cend - cstart;
            const int32_t tol = (m == 1 || m == 2) ? tau : 0;   // the FFT is approximate (error <= 1 per product)
            for (int32_t i = 0; i < N; i++) {
                const int32_t err = abs(int32_t(sumNaive->coefsT[i] - sumRes->coefsT[i]));
                if (m == 1 && err > maxErrFFT) maxErrFFT = err;
                if (m == 2 && err > maxErrFFTLazy) maxErrFFTLazy = err;
                if (err > tol) {
                    cerr << i << " " << sumNaive->coefsT[i] << " vs. " << sumRes->coefsT[i] << endl;
                    dieDramatically(string("Naive sum != ") + sumNames[m] + "\n");
                }
            }
        }
    }
    // cout << "sum of " << tau << " products:" << endl;
    // for (int32_t m = 0; m < nSum; m++)
    //     cout << "  " << sumNames[m] << ": " << (cyc_sum[m] / countSum) << " clock cycles (average)" << endl;
    // cout << "  max error vs. naive: FFT " << maxErrFFT << ", FFT lazy " << maxErrFFTLazy << endl;
    // cout << endl;
    for (int32_t j = 0; j < tau; j++) {
        delete_IntPolynomial(as[j]);
        delete_TorusPolynomial(bs[j]);
    }
    delete_TorusPolynomial(tmp);
    delete_TorusPolynomial(sumNaive);
    delete_TorusPolynomial(sumRes);

    delete_IntPolynomial(a);
    delete_TorusPolynomial(b);
    delete_TorusPolynomial(resNaive);
    delete_TorusPolynomial(resKaratsuba);
    delete_TorusPolynomial(resNTT);
    delete_TorusPolynomial(resToomCook4S);
    return 0;
}
