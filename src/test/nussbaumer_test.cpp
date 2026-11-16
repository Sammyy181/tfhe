#include <gtest/gtest.h>
#include <cstring>
#include "tfhe.h"
#include "nussbaumer.h"

using namespace std;

// All the Nussbaumer + Winograd-Schoolbook computations are exact integer
// computations mod 2^32: unlike the FFT tests, every comparison below is an
// exact equality with the naive/Karatsuba reference.

namespace {

    const int32_t N = NUSS_N;

    void randomIntPolynomial(IntPolynomial *a, int32_t bound) {
        for (int32_t i = 0; i < N; i++) a->coefs[i] = uniformTorus32_distrib(generator) % bound;
    }

    void assertSamePolynomial(const TorusPolynomial *a, const TorusPolynomial *b) {
        for (int32_t i = 0; i < N; i++)
            ASSERT_EQ(a->coefsT[i], b->coefsT[i]) << "coefficient " << i;
    }

    // inverse(forward(a)) = a, including for full-range 32-bit coefficients
    TEST(NussbaumerTest, forwardThenInverseIsIdentity) {
        TorusPolynomial *a = new_TorusPolynomial(N);
        TorusPolynomial *b = new_TorusPolynomial(N);
        NussWord *aNuss = new NussWord[NUSS_SIZE];
        for (int32_t trials = 0; trials < 10; ++trials) {
            torusPolynomialUniform(a);
            nussbaumerForward(aNuss, a->coefsT);
            nussbaumerInverse(b->coefsT, aNuss);
            assertSamePolynomial(a, b);
        }
        delete[] aNuss;
        delete_TorusPolynomial(b);
        delete_TorusPolynomial(a);
    }

    // one product, against the naive multiplication
    TEST(NussbaumerTest, torusPolynomialMulNussWS) {
        IntPolynomial *a = new_IntPolynomial(N);
        TorusPolynomial *b = new_TorusPolynomial(N);
        TorusPolynomial *aB = new_TorusPolynomial(N);
        TorusPolynomial *aBref = new_TorusPolynomial(N);
        for (int32_t trials = 0; trials < 20; ++trials) {
            // small digits (like a gadget decomposition), then larger and
            // larger integers, up to full-range 32-bit coefficients
            if (trials < 10) randomIntPolynomial(a, 512);
            else if (trials < 15) randomIntPolynomial(a, 1 << 20);
            else for (int32_t i = 0; i < N; i++) a->coefs[i] = uniformTorus32_distrib(generator);
            torusPolynomialUniform(b);

            torusPolynomialMultNaive(aBref, a, b);
            torusPolynomialMulNussWS(aB, a, b);
            assertSamePolynomial(aBref, aB);
        }
        delete_TorusPolynomial(aBref);
        delete_TorusPolynomial(aB);
        delete_TorusPolynomial(b);
        delete_IntPolynomial(a);
    }

    // extreme coefficients (largest possible intermediate values)
    TEST(NussbaumerTest, torusPolynomialMulNussWSExtremeValues) {
        IntPolynomial *a = new_IntPolynomial(N);
        TorusPolynomial *b = new_TorusPolynomial(N);
        TorusPolynomial *aB = new_TorusPolynomial(N);
        TorusPolynomial *aBref = new_TorusPolynomial(N);
        const int32_t extremes[3] = {INT32_MIN, INT32_MAX, -1};
        for (int32_t va = 0; va < 3; va++) {
            for (int32_t vb = 0; vb < 3; vb++) {
                for (int32_t i = 0; i < N; i++) {
                    a->coefs[i] = extremes[va];
                    b->coefsT[i] = extremes[vb];
                }
                torusPolynomialMultNaive(aBref, a, b);
                torusPolynomialMulNussWS(aB, a, b);
                assertSamePolynomial(aBref, aB);
            }
        }
        delete_TorusPolynomial(aBref);
        delete_TorusPolynomial(aB);
        delete_TorusPolynomial(b);
        delete_IntPolynomial(a);
    }

    // lazy interpolation: several products accumulated in the Nussbaumer
    // domain, a single inverse transform at the end
    TEST(NussbaumerTest, torusPolynomialAddMulRNussbaumer) {
        const int32_t NBPRODS = 6;
        IntPolynomial *a = new_IntPolynomial(N);
        TorusPolynomial *b = new_TorusPolynomial(N);
        TorusPolynomial *res = new_TorusPolynomial(N);
        TorusPolynomial *ref = new_TorusPolynomial(N);
        NussWord *bNuss = new NussWord[NUSS_SIZE];
        NussWord *accum = new NussWord[NUSS_SIZE];
        for (int32_t trials = 0; trials < 5; ++trials) {
            torusPolynomialClear(ref);
            memset(accum, 0, NUSS_SIZE * sizeof(NussWord));
            for (int32_t p = 0; p < NBPRODS; p++) {
                randomIntPolynomial(a, 512);
                torusPolynomialUniform(b);
                torusPolynomialAddMulRKaratsuba(ref, a, b);

                nussbaumerForward(bNuss, b->coefsT);
                torusPolynomialAddMulRNussbaumer(accum, a, bNuss, 1);
            }
            nussbaumerInverse(res->coefsT, accum);
            assertSamePolynomial(ref, res);
        }
        delete[] accum;
        delete[] bNuss;
        delete_TorusPolynomial(ref);
        delete_TorusPolynomial(res);
        delete_TorusPolynomial(b);
        delete_IntPolynomial(a);
    }

    // the external product, against the same computation with Karatsuba:
    // result_j = sum_p deca_p * gsw[p][j]
    void testExternMul(int32_t k, int32_t l, int32_t Bgbit) {
        TLweParams *tlwe_params = new_TLweParams(N, k, 0., 1.);
        TGswParams *params = new_TGswParams(l, Bgbit, tlwe_params);
        const int32_t kpl = params->kpl;

        TGswSample *gsw = new_TGswSample(params);
        TGswSampleNussbaumer *gswNuss = new_TGswSampleNussbaumer(params);
        TLweSample *accum = new_TLweSample(tlwe_params);
        TorusPolynomial *expected = new_TorusPolynomial_array(k + 1, N);
        IntPolynomial *deca = new_IntPolynomial_array(kpl, N);

        for (int32_t trials = 0; trials < 3; ++trials) {
            // the multiplication does not care whether the samples are valid
            // encryptions: uniformly random polynomials are the hardest case
            for (int32_t p = 0; p < kpl; p++)
                for (int32_t j = 0; j <= k; j++)
                    torusPolynomialUniform(&gsw->all_sample[p].a[j]);
            for (int32_t j = 0; j <= k; j++)
                torusPolynomialUniform(&accum->a[j]);
            tGswToNussbaumerConvert(gswNuss, gsw, params);

            for (int32_t i = 0; i <= k; i++)
                tGswTorus32PolynomialDecompH(deca + i * l, accum->a + i, params);
            for (int32_t j = 0; j <= k; j++) {
                torusPolynomialClear(expected + j);
                for (int32_t p = 0; p < kpl; p++)
                    torusPolynomialAddMulRKaratsuba(expected + j, deca + p, &gsw->all_sample[p].a[j]);
            }

            tGswNussbaumerExternMulToTLwe(accum, gswNuss, params);

            for (int32_t j = 0; j <= k; j++)
                assertSamePolynomial(expected + j, &accum->a[j]);
        }

        delete_IntPolynomial_array(kpl, deca);
        delete_TorusPolynomial_array(k + 1, expected);
        delete_TLweSample(accum);
        delete_TGswSampleNussbaumer(gswNuss);
        delete_TGswSample(gsw);
        delete_TGswParams(params);
        delete_TLweParams(tlwe_params);
    }

    TEST(NussbaumerTest, tGswNussbaumerExternMulToTLwe) {
        testExternMul(1, 3, 10); // the bootstrapping parameters (k=1, l=3, Bg=1024)
        testExternMul(2, 4, 8);
    }

}
