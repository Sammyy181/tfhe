#include "tfhe_gate_bootstrapping_structures.h"
#include "lwekey.h"
#include "lweparams.h"
#include "lwesamples.h"
#include "tlwe.h"
#include "tgsw.h"

TFheGateBootstrappingParameterSet::TFheGateBootstrappingParameterSet(const int32_t ks_t, const int32_t ks_basebit, const LweParams* const in_out_params, const TGswParams* const tgsw_params):
    ks_t(ks_t),
    ks_basebit(ks_basebit),
    in_out_params(in_out_params),
    tgsw_params(tgsw_params)
{}

TFheGateBootstrappingCloudKeySet::TFheGateBootstrappingCloudKeySet(
        const TFheGateBootstrappingParameterSet* const params, 
        const LweBootstrappingKey* const bk,
        const LweBootstrappingKeyFFT* const bkFFT,
        const LweBootstrappingKeyNussbaumer* const bkNuss):
    params(params),bk(bk),bkFFT(bkFFT),bkNuss(bkNuss)
{}

TFheGateBootstrappingSecretKeySet::TFheGateBootstrappingSecretKeySet(
        const TFheGateBootstrappingParameterSet* const params, 
        const LweBootstrappingKey* const bk,
        const LweBootstrappingKeyFFT* const bkFFT,
        const LweBootstrappingKeyNussbaumer* const bkNuss,
        const LweKey* lwe_key,
        const TGswKey* tgsw_key):
    params(params),
    lwe_key(lwe_key),
    tgsw_key(tgsw_key),
    cloud(params,bk,bkFFT,bkNuss) 
{}
