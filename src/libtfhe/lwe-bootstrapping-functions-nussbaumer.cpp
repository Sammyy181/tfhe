/*
 * Bootstrapping Nussbaumer functions
 *
 * Line-by-line counterparts of the functions of lwe-bootstrapping-functions-fft.cpp,
 * where the FFT external product is replaced by the Nussbaumer(16, 8) +
 * Winograd-Schoolbook external product (tGswNussbaumerExternMulToTLwe).
 */

#include <iostream>
#include <cassert>
#include "tfhe.h"
#include "nussbaumer.h"

using namespace std;


// converts the bootstrapping key to the Nussbaumer domain (done once)
// the keyswitch key is shared with bk, not copied
EXPORT LweBootstrappingKeyNussbaumer *new_LweBootstrappingKeyNussbaumer(const LweBootstrappingKey *bk) {
    const LweParams *in_out_params = bk->in_out_params;
    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t n = in_out_params->n;

    // Bootstrapping Key Nussbaumer
    TGswSampleNussbaumer *bkNussbaumer = new_TGswSampleNussbaumer_array(n, bk_params);
    for (int32_t i = 0; i < n; ++i) {
        tGswToNussbaumerConvert(&bkNussbaumer[i], &bk->bk[i], bk_params);
    }

    LweBootstrappingKeyNussbaumer *obj = new LweBootstrappingKeyNussbaumer;
    obj->in_out_params = in_out_params;
    obj->bk_params = bk_params;
    obj->accum_params = accum_params;
    obj->extract_params = extract_params;
    obj->bkNussbaumer = bkNussbaumer;
    obj->ks = bk->ks;
    return obj;
}

EXPORT void delete_LweBootstrappingKeyNussbaumer(LweBootstrappingKeyNussbaumer *obj) {
    delete_TGswSampleNussbaumer_array(obj->in_out_params->n, (TGswSampleNussbaumer *) obj->bkNussbaumer);
    delete obj;
}


void tfhe_MuxRotate_Nussbaumer(TLweSample *result, const TLweSample *accum, const TGswSampleNussbaumer *bki,
                               const int32_t barai, const TGswParams *bk_params) {
    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    tLweMulByXaiMinusOne(result, barai, accum, bk_params->tlwe_params);
    // temp *= BKi
    tGswNussbaumerExternMulToTLwe(result, bki, bk_params);
    // ACC += temp
    tLweAddTo(result, accum, bk_params->tlwe_params);
}


/**
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW Nussbaumer samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
EXPORT void tfhe_blindRotate_Nussbaumer(TLweSample *accum,
                                        const TGswSampleNussbaumer *bk,
                                        const int32_t *bara,
                                        const int32_t n,
                                        const TGswParams *bk_params) {

    TLweSample *temp = new_TLweSample(bk_params->tlwe_params);
    TLweSample *temp2 = temp;
    TLweSample *temp3 = accum;

    for (int32_t i = 0; i < n; i++) {
        const int32_t barai = bara[i];
        if (barai == 0) continue; //indeed, this is an easy case!

        tfhe_MuxRotate_Nussbaumer(temp2, temp3, bk + i, barai, bk_params);
        swap(temp2, temp3);
    }
    if (temp3 != accum) {
        tLweCopy(accum, temp3, bk_params->tlwe_params);
    }

    delete_TLweSample(temp);
}


/**
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW Nussbaumer samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
EXPORT void tfhe_blindRotateAndExtract_Nussbaumer(LweSample *result,
                                                  const TorusPolynomial *v,
                                                  const TGswSampleNussbaumer *bk,
                                                  const int32_t barb,
                                                  const int32_t *bara,
                                                  const int32_t n,
                                                  const TGswParams *bk_params) {

    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t N = accum_params->N;
    const int32_t _2N = 2 * N;

    // Test polynomial
    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    // Accumulator
    TLweSample *acc = new_TLweSample(accum_params);

    // testvector = X^{2N-barb}*v
    if (barb != 0) torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    else torusPolynomialCopy(testvectbis, v);
    tLweNoiselessTrivial(acc, testvectbis, accum_params);
    // Blind rotation
    tfhe_blindRotate_Nussbaumer(acc, bk, bara, n, bk_params);
    // Extraction
    tLweExtractLweSample(result, acc, extract_params, accum_params);

    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}


/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void tfhe_bootstrap_woKS_Nussbaumer(LweSample *result,
                                           const LweBootstrappingKeyNussbaumer *bk,
                                           Torus32 mu,
                                           const LweSample *x) {

    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk->accum_params;
    const LweParams *in_params = bk->in_out_params;
    const int32_t N = accum_params->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = in_params->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[N];

    // Modulus switching
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    for (int32_t i = 0; i < n; i++) {
        bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
    }

    // the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) testvect->coefsT[i] = mu;

    // Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_Nussbaumer(result, testvect, bk->bkNussbaumer, barb, bara, n, bk_params);

    delete[] bara;
    delete_TorusPolynomial(testvect);
}


/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void tfhe_bootstrap_Nussbaumer(LweSample *result,
                                      const LweBootstrappingKeyNussbaumer *bk,
                                      Torus32 mu,
                                      const LweSample *x) {

    LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

    tfhe_bootstrap_woKS_Nussbaumer(u, bk, mu, x);
    // Key switching
    lweKeySwitch(result, bk->ks, u);

    delete_LweSample(u);
}
