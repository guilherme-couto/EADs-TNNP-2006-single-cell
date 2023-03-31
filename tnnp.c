/*-----------------------------------------------------
ten Tusscher model 2006
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

/*------------------------------------------------------
Electrophysiological parameters for ten Tusscher model 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D
--------------------------------------------------------*/
// Constants
double R = 8314.472;      // Gas constant -> (???) [8.314472 J/(K*mol)]
double T = 310.0;         // Temperature -> K
double F = 96485.3415;    // Faraday constant -> (???) [96.4867 C/mmol]
double RTONF = 26.713761; // R*T/F -> (???)

// Tissue properties
double beta = 1400.0;       // Surface area-to-volume ratio -> cm^-1
double Cm = 0.185;          // Cell capacitance per unit surface area -> uF/ (???)^2 (ten Tusscher)
double sigma_long = 1.334;  // Conductivity -> mS/cm
double sigma_trans = 0.176; // Conductivity -> mS/cm
// double sigma = 0.1;         // Conductivity (isotropic) -> mS/cm (Sachetto)
// double beta = 0.14;            // Surface area-to-volume ratio -> um^-1
// double Cm = 0.185;             // Cell capacitance per unit surface area -> uF/ (???)^2 (um???) (ten Tusscher) (Sachetto)
// double sigma = 0.00001;        // Conductivity (isotropic) -> mS/um

// Intracellular volumes
double V_C = 0.016404;    // Cellular volume -> (???) [16404 um^3]
double V_SR = 0.001094;   // Sarcoplasmic reticulum volume -> (???) [1094 um^3]
double V_SS = 0.00005468; // Subsarcolemmal space volume -> (???) [54.68 um^3]

// External concentrations
double K_o = 5.4;  // Extracellular potassium (K+) concentration -> mM
double Na_o = 140; // Extracellular sodium (Na+) concentration -> mM
double Ca_o = 2.0; // Extracellular calcium (Ca++) concentration -> mM

// Parameters for currents
double G_Na = 14.838; // Maximal I_Na (sodium current) conductance -> nS/pF
double G_K1 = 5.405;  // Maximal I_K1 (late rectifier potassium current) conductance -> nS/pF
double G_to = 0.294;  // Maximal I_to (transient outward potassium current) conductance -> nS/pF (epi and M cells)
// double G_to = 0.073;        // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)
double G_Kr = 0.2*0.153; // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF
// double G_Ks = 1.0*0.392; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)
double G_Ks = 1*0.098;        // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells) (Sachetto)
double p_KNa = 0.03;  // Relative I_Ks permeability to Na+ over K+ -> dimensionless
// double G_CaL = 2.0*3.98e-5;     // Maximal I_CaL (L-type calcium current) conductance -> cm/ms/uF
double G_CaL = 6*7.96e-5;     // Maximal I_CaL (L-type calcium current) conductance -> mm^3/ms/uF (mTP06b) (M cells)
double k_NaCa = 1000.0;     // Maximal I_NaCa (Na+/Ca++ exchanger current) -> pA/pF
double gamma_I_NaCa = 0.35; // Voltage dependence parameter of I_NaCa -> dimensionless
double K_mCa = 1.38;        // Half-saturation constant of I_NaCa for intracellular Ca++ -> mM
double K_mNa_i = 87.5;      // Half-saturation constant of I_NaCa for intracellular Na+ -> mM
double k_sat = 0.1;         // Saturation factor for I_NaCa -> dimensionless
double alpha = 2.5;         // Factor enhancing outward nature of I_NaCa -> dimensionless
double P_NaK = 2.724;       // Maximal I_NaK (Na+/K+ pump current) -> pA/pF
double K_mK = 1.0;          // Half-saturation constant of I_NaK for Ko -> mM
double K_mNa = 40.0;        // Half-saturation constant of I_NaK for intracellular Na+ -> mM
double G_pK = 0.0146;       // Maximal I_pK (plateau potassium current) conductance -> nS/pF
double G_pCa = 0.1238;      // Maximal I_pCa (plateau calcium current) conductance -> nS/pF
double K_pCa = 0.0005;      // Half-saturation constant of I_pCa for intracellular Ca++ -> mM
double G_bNa = 0.00029;     // Maximal I_bNa (sodium background current) conductance -> nS/pF
double G_bCa = 0.000592;    // Maximal I_bCa (calcium background current) conductance -> nS/pF

// Intracellular calcium flux dynamics
double V_maxup = 0.006375; // Maximal I_up -> mM/ms
double K_up = 0.00025;     // Half-saturation constant of I_up -> mM
double V_rel = 0.102;      // Maximal I_rel conductance -> mM/ms
double k1_prime = 0.15;    // R to O and RI to I I_rel transition rate -> mM^-2*ms^-1
double k2_prime = 0.045;   // O to I  and R to RI I_rel transition rate -> mM^-1*ms^-1
double k3 = 0.06;          // O to R and I to RI I_rel transition rate -> ms^-1
double k4 = 0.005;         // I to O and RI to I I_rel transition rate -> ms^-1
double EC = 1.5;           // Half-saturation constant of k_Ca_SR -> mM
double max_SR = 2.5;       // Maximum value of k_Ca_SR -> dimensionless
double min_SR = 1.0;       // Minimum value of k_Ca_SR -> dimensionless
double V_leak = 0.00036;   // Maximal I_leak conductance -> mM/ms
double V_xfer = 0.0038;    // Maximal I_xfer conductance -> mM/ms

// Calcium buffering dynamics
double Buf_C = 0.2;       // Total cytoplasmic buffer concentration -> mM
double K_bufc = 0.001;    // Half-saturation constant of cytoplasmic buffers -> mM
double Buf_SR = 10.0;     // Total sarcoplasmic reticulum buffer concentration -> mM
double K_bufsr = 0.3;     // Half-saturation constant of sarcoplasmic reticulum buffers -> mM
double Buf_SS = 0.4;      // Total subspace buffer concentration -> mM
double K_bufss = 0.00025; // Half-saturation constant of subspace buffer -> mM

/*----------------------------------------
Initial Conditions
-----------------------------------------*/
double V_init = -86.2;       // Initial membrane potential -> mV
double X_r1_init = 0.0;   // Initial rapid time-dependent potassium current Xr1 gate -> dimensionless
double X_r2_init = 1.0;    // Initial rapid time-dependent potassium current Xr2 gate -> dimensionless
double X_s_init = 0.0;     // Initial slow time-dependent potassium current Xs gate -> dimensionless
double m_init = 0.0;      // Initial fast sodium current m gate -> dimensionless
double h_init = 0.75;       // Initial fast sodium current h gate -> dimensionless
double j_init = 0.75;       // Initial fast sodium current j gate -> dimensionless
double d_init = 0.0;     // Initial L-type calcium current d gate -> dimensionless
double f_init = 1.0;       // Initial L-type calcium current f gate -> dimensionless
double f2_init = 1.0;      // Initial L-type calcium current f2 gate -> dimensionless
double fCass_init = 1.0;   // Initial L-type calcium current fCass gate -> dimensionless
double s_init = 1.0;     // Initial transient outward current s gate -> dimensionless
double r_init = 0.0;      // Initial transient outward current r gate -> dimensionless
double Ca_i_init = 0.00007;  // Initial intracellular Ca++ concentration -> mM
double Ca_SR_init = 1.3;     // Initial sarcoplasmic reticulum Ca++ concentration -> mM
double Ca_SS_init = 0.00007;  // Initial subspace Ca++ concentration -> mM
double R_prime_init = 1.0; // Initial ryanodine receptor -> dimensionless
double Na_i_init = 7.67;     // Initial intracellular Na+ concentration -> mM
double K_i_init = 138.3;     // Initial intracellular K+ concentration -> mM

void initialize_from_file(char *filename, double *V, double *X_r1, double *X_r2, double *X_s, double *m, double *h, double *j, double *d, double *f, double *f2, double *fCass, double *s, double *r, double *Ca_i, double *Ca_SR, double *Ca_SS, double *R_prime, double *Na_i, double *K_i)
{
    char *ptr;
    FILE *fp_read;
    fp_read = fopen(filename, "r");
    char line[20];
    fgets(line, 20, fp_read);
    double num = strtod(line, &ptr);
    *V = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *X_r1 = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *X_r2 = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *X_s = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *m = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *h = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *j = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *d = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *f = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *f2 = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *fCass = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *s = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *r = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *Ca_i = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *Ca_SR = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *Ca_SS = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *R_prime = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *Na_i = num;
    fgets(line, 20, fp_read);
    num = strtod(line, &ptr);
    *K_i = num;

}

/*----------------------------------------
Currents functions
------------------------------------------*/
// Reversal potentials for Na+, K+ and Ca++
double E_Na(double Na_i)
{
    return RTONF * log(Na_o / Na_i);
}
double E_K(double K_i)
{
    return RTONF * log(K_o / K_i);
}
double E_Ca(double Ca_i)
{
    return 0.5 * RTONF * log(Ca_o / Ca_i);
}

// Reversal potential for Ks
double E_Ks(double K_i, double Na_i)
{
    return RTONF * log((K_o + p_KNa * Na_o) / (K_i + p_KNa * Na_i));
}

// Fast sodium (Na+) current
double I_Na(double V, double m, double h, double j, double Na_i)
{
    return G_Na * pow(m, 3.0) * h * j * (V - E_Na(Na_i));
}
double m_inf(double V)
{
    return 1.0 / pow((1.0 + exp((-56.86 - V) / 9.03)), 2.0);
}
double alpha_m(double V)
{
    return 1.0 / (1.0 + exp((-60.0 - V) / 5.0));
}
double beta_m(double V)
{
    return (0.1 / (1.0 + exp((V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((V - 50.0) / 200.0)));
}
double tau_m(double V)
{
    return alpha_m(V) * beta_m(V);
}
double h_inf(double V)
{
    return 1.0 / pow((1.0 + exp((V + 71.55) / 7.43)), 2.0);
}
double alpha_h(double V)
{
    if (V >= -40.0)
    {
        return 0.0;
    }
    else
    {
        return 0.057 * exp(-(80.0 + V) / 6.8);
    }
}
double beta_h(double V)
{
    if (V >= -40)
    {
        return 0.77 / (0.13 * (1.0 + exp((V + 10.66) / (-11.1))));
    }
    else
    {
        return 2.7 * exp(0.079 * V) + 3.1e5 * exp(0.3485 * V);
    }
}
double tau_h(double V)
{
    return 1.0 / (alpha_h(V) + beta_h(V));
}
double j_inf(double V)
{
    return 1.0 / pow((1.0 + exp((V + 71.55) / 7.43)), 2.0);
}
double alpha_j(double V)
{
    if (V >= -40.0)
    {
        return 0.0;
    }
    else
    {
        return ((-25428.0 * exp(0.2444 * V) - (6.948e-6 * exp((-0.04391) * V))) * (V + 37.78)) / (1.0 + exp(0.311 * (V + 79.23)));
    }
}
double beta_j(double V)
{
    if (V >= -40.0)
    {
        return (0.6 * exp(0.057 * V)) / (1.0 + exp(-0.1 * (V + 32.0)));
    }
    else
    {
        return (0.02424 * exp(-0.01052 * V)) / (1.0 + exp(-0.1378 * (V + 40.14)));
    }
}
double tau_j(double V)
{
    return 1.0 / (alpha_j(V) + beta_j(V));
}

// L-type Ca2+ current
double I_CaL(double V, double d, double f, double f2, double fCass, double Ca_SS)
{
    return G_CaL * d * f * f2 * fCass * 4.0 * (V - 15.0) * (pow(F, 2.0) / (R * T)) * (0.25 * Ca_SS * exp(2.0 * (V - 15.0) * F / (R * T)) - Ca_o) / (exp(2.0 * (V - 15.0) * F / (R * T)) - 1.0);
}
double d_inf(double V)
{
    return 1.0 / (1.0 + exp((-8.0 - V) / 7.5));
}
double alpha_d(double V)
{
    return (1.4 / (1.0 + exp((-35.0 - V) / 13.0))) + 0.25;
}
double beta_d(double V)
{
    return 1.4 / (1.0 + exp((V + 5.0) / 5.0));
}
double gamma_d(double V)
{
    return 1.0 / (1.0 + exp((50.0 - V) / 20.0));
}
double tau_d(double V)
{
    return alpha_d(V) * beta_d(V) + gamma_d(V);
}
double f_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 20.0) / 7.0));
}
double alpha_f(double V)
{
    return 1102.5 * exp(-(pow((V + 27.0), 2.0)) / 225.0);
}
double beta_f(double V)
{
    return 200.0 / (1.0 + exp((13.0 - V) / 10.0));
}
double gamma_f(double V)
{
    return 180.0 / (1.0 + exp((V + 30.0) / 10.0)) + 20.0;
}
double tau_f(double V)
{
    return (alpha_f(V) + beta_f(V) + gamma_f(V))*0.5;
}
double f2_inf(double V)
{
    return 0.67 / (1.0 + exp((V + 35.0) / 7.0)) + 0.33;
}
double alpha_f2(double V)
{
    return 600.0 * exp(-(pow((V + 25.0), 2.0)) / 170.0);
}
double beta_f2(double V)
{
    return 31.0 / (1.0 + exp((25.0 - V) / 10.0));
}
double gamma_f2(double V)
{
    return 16.0 / (1.0 + exp((V + 30.0) / 10.0));
}
double tau_f2(double V)
{
    return alpha_f2(V) + beta_f2(V) + gamma_f2(V);
}
double fCass_inf(double Ca_SS)
{
    return 0.6 / (1.0 + pow((Ca_SS / 0.05), 2.0)) + 0.4;
}
double tau_fCass(double Ca_SS)
{
    return 80.0 / (1.0 + pow((Ca_SS / 0.05), 2.0)) + 2.0;
}

// Transient outward current
double I_to(double V, double r, double s, double K_i)
{
    return G_to * r * s * (V - E_K(K_i));
}
double r_inf(double V)
{
    return 1.0 / (1.0 + exp((20.0 - V) / 6.0));
}
double tau_r(double V)
{
    return 9.5 * exp(-(pow((V + 40.0), 2.0)) / 1800.0) + 0.8;
}
/* for epicardial and M cells */
double s_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 20.0) / 5.0));
}
double tau_s(double V)
{
    return 85.0 * exp(-(pow((V + 45.0), 2.0)) / 320.0) + 5.0 / (1.0 + exp((V - 20.0) / 5.0)) + 3.0;
}

/* for endocardial cells */
/* double s_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 28.0) / 5.0));
}
double tau_s(double V)
{
    return 1000.0 * exp(-(pow((V + 67.0), 2.0)) / 1000.0) + 8.0;
}
 */

// Slow delayed rectifier current
double I_Ks(double V, double X_s, double K_i, double Na_i)
{
    return G_Ks * pow(X_s, 2.0) * (V - E_Ks(K_i, Na_i));
}
double x_s_inf(double V)
{
    return 1.0 / (1.0 + exp((-5.0 - V) / 14.0));
}
double alpha_x_s(double V)
{
    return 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0));
}
double beta_x_s(double V)
{
    return 1.0 / (1.0 + exp((V - 35.0) / 15.0));
}
double tau_x_s(double V)
{
    return alpha_x_s(V) * beta_x_s(V) + 80.0;
}

// Rapid delayed rectifier current
double I_Kr(double V, double X_r1, double X_r2, double K_i)
{
    return G_Kr * sqrt(K_o / 5.4) * X_r1 * X_r2 * (V - E_K(K_i));
}
double x_r1_inf(double V)
{
    return 1.0 / (1.0 + exp((-26.0 - V) / 7.0));
}
double alpha_x_r1(double V)
{
    return 450.0 / (1.0 + exp((-45.0 - V) / 10.0));
}
double beta_x_r1(double V)
{
    return 6.0 / (1.0 + exp((V + 30.0) / 11.5));
}
double tau_x_r1(double V)
{
    return alpha_x_r1(V) * beta_x_r1(V);
}
double x_r2_inf(double V)
{
    return 1.0 / (1.0 + exp((V + 88.0) / 24.0));
}
double alpha_x_r2(double V)
{
    return 3.0 / (1.0 + exp((-60.0 - V) / 20.0));
}
double beta_x_r2(double V)
{
    return 1.12 / (1.0 + exp((V - 60.0) / 20.0));
}
double tau_x_r2(double V)
{
    return alpha_x_r2(V) * beta_x_r2(V);
}

// Inward rectifier K+ current
double alpha_K1(double V, double K_i)
{
    return 0.1 / (1.0 + exp(0.06 * (V - E_K(K_i) - 200.0)));
}
double beta_K1(double V, double K_i)
{
    return (3.0 * exp(0.0002 * (V - E_K(K_i) + 100.0)) + exp(0.1 * (V - E_K(K_i) - 10.0))) / (1.0 + exp(-0.5 * (V - E_K(K_i))));
}
double x_K1_inf(double V, double K_i)
{
    return alpha_K1(V, K_i) / (alpha_K1(V, K_i) + beta_K1(V, K_i));
}
double I_K1(double V, double K_i)
{
    return G_K1 * x_K1_inf(V, K_i) * (V - E_K(K_i));
}

// Na+/Ca++ exchanger current
double I_NaCa(double V, double Na_i, double Ca_i)
{
    return k_NaCa * (1.0 / (pow(K_mNa_i, 3.0) + pow(Na_o, 3.0))) * (1.0 / (K_mCa + Ca_o)) * (1.0 / (1.0 + k_sat * exp((gamma_I_NaCa - 1) * V * F / (R * T)))) * (exp(gamma_I_NaCa * V * F / (R * T)) * pow(Na_i, 3.0) * Ca_o - exp((gamma_I_NaCa - 1) * V * F / (R * T)) * pow(Na_o, 3.0) * Ca_i * 2.5);
}

// Na+/K+ pump current
double I_NaK(double V, double Na_i)
{
    return P_NaK * (K_o / (K_o + K_mK)) * (Na_i / (Na_i + K_mNa)) * (1.0 / (1.0 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0353 * exp(-V * F / (R * T))));
}

// I_pCa
double I_pCa(double V, double Ca_i)
{
    return (G_pCa * Ca_i) / (K_pCa + Ca_i);
}

// I_pK
double I_pK(double V, double K_i)
{
    return G_pK * (V - E_K(K_i)) * (1.0 / (1.0 + exp((25.0 - V) / 5.98)));
}

// Background currents
double I_bNa(double V, double Na_i)
{
    return G_bNa * (V - E_Na(Na_i));
}
double I_bCa(double V, double Ca_i)
{
    return G_bCa * (V - E_Ca(Ca_i));
}

// Calcium dynamics
double I_leak(double Ca_SR, double Ca_i)
{
    return V_leak * (Ca_SR - Ca_i);
}
double I_up(double Ca_i)
{
    return V_maxup / (1.0 + (pow(K_up, 2.0) / pow(Ca_i, 2.0)));
}
double k_casr(double Ca_SR)
{
    return max_SR - ((max_SR - min_SR) / (1.0 + pow((EC / Ca_SR), 2.0)));
}
double k1(double Ca_SR)
{
    return k1_prime / k_casr(Ca_SR);
}
double O(double Ca_SR, double Ca_SS, double R_prime)
{
    return (k1(Ca_SR) * pow(Ca_SS, 2.0) * R_prime) / (k3 + (k1(Ca_SR) * pow(Ca_SS, 2.0)));
}
double I_rel(double Ca_SR, double Ca_SS, double R_prime)
{
    return V_rel * O(Ca_SR, Ca_SS, R_prime) * (Ca_SR - Ca_SS);
}
double I_xfer(double Ca_SS, double Ca_i)
{
    return V_xfer * (Ca_SS - Ca_i);
}
double k2(double Ca_SR)
{
    return k2_prime * k_casr(Ca_SR);
}
double Ca_ibufc(double Ca_i)
{
    return Buf_C * Ca_i / (Ca_i + K_bufc);
}
double Ca_srbufsr(double Ca_SR)
{
    return Buf_SR * Ca_SR / (Ca_SR + K_bufsr);
}
double Ca_ssbufss(double Ca_SS)
{
    return Buf_SS * Ca_SS / (Ca_SS + K_bufss);
}
double bc(double Ca_i, double dCa_i)
{
    return Buf_C - Ca_ibufc(Ca_i) - dCa_i - Ca_i + K_bufc;
}
double cc(double Ca_i, double dCa_i)
{
    return K_bufc * (Ca_ibufc(Ca_i) + dCa_i + Ca_i);
}
double bjsr(double Ca_SR, double dCa_SR)
{
    return Buf_SR - Ca_srbufsr(Ca_SR) - dCa_SR - Ca_SR + K_bufsr;
}
double cjsr(double Ca_SR, double dCa_SR)
{
    return K_bufsr * (Ca_srbufsr(Ca_SR) + dCa_SR + Ca_SR);
}
double bcss(double Ca_SS, double dCa_SS)
{
    return Buf_SS - Ca_ssbufss(Ca_SS) - dCa_SS - Ca_SS + K_bufss;
}
double ccss(double Ca_SS, double dCa_SS)
{
    return K_bufss * (Ca_ssbufss(Ca_SS) + dCa_SS + Ca_SS);
}

/*----------------------------------------
Simulation parameters
-----------------------------------------*/

// double simulation_time = 1900000; // End time -> ms 30 min

// double simulation_time = 30000;
double dt = 0.02;              // Time step -> ms
int L = 2;                     // Length of the domain (square tissue) -> cm

/*----------------------------------------
Main function
-----------------------------------------*/
int main(int argc, char *argv[])
{

    double I_total = 0.0;
    double V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i;

    // Initialize variables
    V = V_init;
    X_r1 = X_r1_init;
    X_r2 = X_r2_init;
    X_s = X_s_init;
    m = m_init;
    h = h_init;
    j = j_init;
    d = d_init;
    f = f_init;
    f2 = f2_init;
    fCass = fCass_init;
    s = s_init;
    r = r_init;
    Ca_i = Ca_i_init;
    Ca_SR = Ca_SR_init;
    Ca_SS = Ca_SS_init;
    R_prime = R_prime_init;
    Na_i = Na_i_init;
    K_i = K_i_init;
    
    // initialize_from_file("last-teste-3.txt", &V, &X_r1, &X_r2, &X_s, &m, &h, &j, &d, &f, &f2, &fCass, &s, &r, &Ca_i, &Ca_SR, &Ca_SS, &R_prime, &Na_i, &K_i);
    
    double simulation_time = atof(argv[1]); // End time -> ms
    double time = 0.0;
    double last_time = 0.0;
    int step = 0;
    double Istim;
    double dR_prime_dt, dCa_SR, dCa_SS, dCa_i, dNa_i_dt, dK_i_dt;

    // Stimulation parameters
    double stimduration = 2.0;
    // double stimstrength=-38;
    double stimstrength = -20; // EAD
    double tbegin = 0;
    double tend = tbegin + stimduration;
    double period = atof(argv[2]);
    double period2 = atof(argv[3]);
    int stim_count = 1;
    int stim_count2 = 0;
    int stim_count_max = atof(argv[4]);
    int stim_count_max2 = atof(argv[5]);
    double time_until_next_stim = atof(argv[6]);

    // Files to save V and time
    FILE *fp_v, *fp_t, *fp_last, *fp_last_ap = fopen("last-ap.txt", "w"), *fp_last_ap_time = fopen("last-ap-time.txt", "w");
    fp_v = fopen("V.txt", "w");
    fp_t = fopen("t.txt", "w");
    

    // Forward Euler
    for (step = 0; time < simulation_time; step++)
    {
        // Save V and time
        fprintf(fp_v, "%lf\n", V);
        fprintf(fp_t, "%lf\n", time);

        // Stimulus
        if (time >= tbegin && time <= tend)
        {
            Istim = stimstrength;
        }
        else
        {
            Istim = 0.0;
        }
        if (time > tend && (stim_count < stim_count_max ||  stim_count2 < stim_count_max2))
        {
            if (stim_count < stim_count_max)
            {
                stim_count += 1;
                tbegin = tbegin + period;
            }
            else if (stim_count2 < stim_count_max2)
            {
                if (stim_count == stim_count_max && stim_count2 == 0)
                {
                    tbegin = tbegin + time_until_next_stim;
                }
                stim_count2 += 1;
                tbegin = tbegin + period2;
            }
            tend = tbegin + stimduration;
        }

        // Update total current
        I_total = Istim + I_Na(V, m, h, j, Na_i) + I_bNa(V, Na_i) + I_K1(V, K_i) + I_to(V, r, s, K_i) + I_Kr(V, X_r1, X_r2, K_i) + I_Ks(V, X_s, K_i, Na_i) + I_CaL(V, d, f, f2, fCass, Ca_SS) + I_NaK(V, Na_i) + I_NaCa(V, Na_i, Ca_i) + I_pCa(V, Ca_i) + I_pK(V, K_i) + I_bCa(V, Ca_i);
        

        // Update concentrations
        dR_prime_dt = ((-k2(Ca_SS)) * Ca_SS * R_prime) + (k4 * (1.0 - R_prime));
        dCa_i = dt * ((-(I_bCa(V, Ca_i) + I_pCa(V, Ca_i) - 2.0 * I_NaCa(V, Na_i, Ca_i)) * (1.0 / (2.0 * V_C * F)) * Cm) - (I_up(Ca_i) - I_leak(Ca_SR, Ca_i)) * (V_SR / V_C) + I_xfer(Ca_SS, Ca_i));
        dCa_SR = dt * (I_up(Ca_i) - I_leak(Ca_SR, Ca_i) - I_rel(Ca_SR, Ca_SS, R_prime));
        dCa_SS = dt * (-I_xfer(Ca_SS, Ca_i) * (V_C / V_SS) + I_rel(Ca_SR, Ca_SS, R_prime) * (V_SR / V_SS) + (-I_CaL(V, d, f, f2, fCass, Ca_SS) * (1.0 / (2.0 * V_C * F)) * Cm));
        dNa_i_dt = -(I_Na(V, m, h, j, Na_i) + I_bNa(V, Na_i) + (3.0 * I_NaK(V, Na_i)) + (3.0 * I_NaCa(V, Na_i, Ca_i))) * (1 / (V_C * F)) * Cm;
        dK_i_dt = -(Istim + I_K1(V, K_i) + I_to(V, r, s, K_i) + I_Kr(V, X_r1, X_r2, K_i) + I_Ks(V, X_s, K_i, Na_i) + I_pK(V, K_i) - (2.0 * I_NaK(V, Na_i))) * (1 / (V_C * F)) * Cm;

        R_prime = R_prime + dR_prime_dt * dt;
        Ca_SR = (sqrt(pow(bjsr(Ca_SR, dCa_SR), 2.0) + 4.0 * cjsr(Ca_SR, dCa_SR)) - bjsr(Ca_SR, dCa_SR)) / 2.0;
        Ca_SS = (sqrt(pow(bcss(Ca_SS, dCa_SS), 2.0) + 4.0 * ccss(Ca_SS, dCa_SS)) - bcss(Ca_SS, dCa_SS)) / 2.0;
        Ca_i = (sqrt(pow(bc(Ca_i, dCa_i), 2.0) + 4.0 * cc(Ca_i, dCa_i)) - bc(Ca_i, dCa_i)) / 2.0;
        Na_i = Na_i + dNa_i_dt * dt;
        K_i = K_i + dK_i_dt * dt;

        // Update gating variables - Rush Larsen
        X_r1 = x_r1_inf(V) - (x_r1_inf(V) - X_r1) * exp(-dt / tau_x_r1(V));
        X_r2 = x_r2_inf(V) - (x_r2_inf(V) - X_r2) * exp(-dt / tau_x_r2(V));
        X_s = x_s_inf(V) - (x_s_inf(V) - X_s) * exp(-dt / tau_x_s(V));
        r = r_inf(V) - (r_inf(V) - r) * exp(-dt / tau_r(V));
        s = s_inf(V) - (s_inf(V) - s) * exp(-dt / tau_s(V));
        m = m_inf(V) - (m_inf(V) - m) * exp(-dt / tau_m(V));
        h = h_inf(V) - (h_inf(V) - h) * exp(-dt / tau_h(V));
        j = j_inf(V) - (j_inf(V) - j) * exp(-dt / tau_j(V));
        d = d_inf(V) - (d_inf(V) - d) * exp(-dt / tau_d(V));
        f = f_inf(V) - (f_inf(V) - f) * exp(-dt / tau_f(V));
        f2 = f2_inf(V) - (f2_inf(V) - f2) * exp(-dt / tau_f2(V));
        fCass = fCass_inf(V) - (fCass_inf(V) - fCass) * exp(-dt / tau_fCass(V));

        // Update voltage
        V = V + (-I_total) * dt;

        time += dt;

        if((stim_count2 == stim_count_max2 || stim_count2 == stim_count_max2 - 1) && last_time <= 1000)
        {
            fprintf(fp_last_ap, "%lf\n", V);
            fprintf(fp_last_ap_time, "%lf\n", last_time);
            last_time += dt;
        }
    }

    fclose(fp_last_ap);
    fclose(fp_last_ap_time);

    fp_last = fopen("last-teste-3.txt", "w");
    fprintf(fp_last, "%lf\n", V);
    fprintf(fp_last, "%lf\n", X_r1);
    fprintf(fp_last, "%lf\n", X_r2);
    fprintf(fp_last, "%lf\n", X_s);
    fprintf(fp_last, "%lf\n", m);
    fprintf(fp_last, "%lf\n", h);
    fprintf(fp_last, "%lf\n", j);
    fprintf(fp_last, "%lf\n", d);
    fprintf(fp_last, "%lf\n", f);
    fprintf(fp_last, "%lf\n", f2);
    fprintf(fp_last, "%lf\n", fCass);
    fprintf(fp_last, "%lf\n", s);
    fprintf(fp_last, "%lf\n", r);
    fprintf(fp_last, "%lf\n", Ca_i);
    fprintf(fp_last, "%lf\n", Ca_SR);
    fprintf(fp_last, "%lf\n", Ca_SS);
    fprintf(fp_last, "%lf\n", R_prime);
    fprintf(fp_last, "%lf\n", Na_i);
    fprintf(fp_last, "%lf\n", K_i);
    fclose(fp_last);


    return 0;
}
