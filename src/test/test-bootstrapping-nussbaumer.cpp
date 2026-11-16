#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}

#ifndef NDEBUG
extern const TLweKey *debug_accum_key;
extern const LweKey *debug_extract_key;
extern const LweKey *debug_in_key;
#endif

int32_t main(int32_t argc, char **argv){
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif
    const int32_t nb_samples = 64;
    const Torus32 mu_boot = modSwitchToTorus32(1, 8);

    // generate params
    int32_t minimum_lambda = 100;
    TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *in_out_params = params->in_out_params;
    // generate the secret keyset
    TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);

    // generate input samples
    LweSample *test_in = new_LweSample_array(nb_samples, in_out_params);
    for (int32_t i = 0; i < nb_samples; ++i){
        lweSymEncrypt(test_in + i, modSwitchToTorus32(i, nb_samples), 0.01, keyset->lwe_key);
    }
    // output samples
    LweSample *test_out = new_LweSample_array(nb_samples, in_out_params);

    // bootstrap input samples
    cout << "starting bootstrapping..." << endl;
    clock_t begin = clock();
    for (int32_t i = 0; i < nb_samples; ++i){
        tfhe_bootstrap_Nussbaumer(test_out + i, keyset->cloud.bkNuss, mu_boot, test_in + i);
    }
    clock_t end = clock();
    cout << "finished " << nb_samples << "bootstrappings" << endl;
    cout << "time per bootstrapping (microsecs)... " << (end - begin) / double(nb_samples) << endl;
    
    for (int32_t i = 0; i < nb_samples; ++i) {
        Torus32 phase = lwePhase(test_out + i, keyset->lwe_key);
        cout << "phase " << i << " = " << t32tod(phase) << endl;
    }

    delete_LweSample_array(nb_samples, test_out);
    delete_LweSample_array(nb_samples, test_in);
    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}