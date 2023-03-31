/*-----------------------------------------------------
Bueno_Orovio (Minimal VEntricular) model
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

// Parameters to reproduce M cells
double u_o = 0;
double u_u = 1.61;
double theta_v = 0.3;
double theta_w = 0.13;
double theta_vminus = 0.1;
double theta_o = 0.005;
double tau_v1minus = 80;
double tau_v2minus = 1.4506;
double tau_vplus = 1.4506;
double tau_w1minus = 70;
double tau_w2minus = 8;
double k_wminus = 200;
double u_wminus = 0.016;
double tau_wplus = 280;
double tau_fi = 0.078;
double tau_o1 = 410;
double tau_o2 = 7;
double tau_so1 = 91;
double tau_so2 = 0.8;
double k_so = 2.1;
double u_so = 0.6;
double tau_s1 = 2.7342;
double tau_s2 = 4;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 3.3849;
double tau_winf = 0.01;
double w_infstar = 0.5;

double D = 1.171; // +- 0.0221 cm^2/s  human ventricular diffusion coefficient

// Standar Heaviside function
double H(double x, double y)
{
    if (x > y)
    {
        return 1;
    }
    else if (x < y)
    {
        return 0;
    }
    else
        return 0.5;
}

double v_inf_function(double x, double y)
{
    if (x < y)
    {
        return 1;
    }
    else
        return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 9)
    {
        fprintf(stderr, "Usage: %s <simulation time (ms)> <S1 period (ms)> <S2 period (ms)> <S1 total stims> <S2 total stims> <time between S1 and S2 (ms)> <alpha Ca> <beta K>\n", argv[0]);
        exit(1);
    }

    // Pacing parameters
    double simulation_time = atof(argv[1]);
    double period = atof(argv[2]);
    double period2 = atof(argv[3]);
    int stim_count_max = atof(argv[4]);
    int stim_count_max2 = atof(argv[5]);
    double time_until_next_stim = atof(argv[6]);

    // Discretization
    double dt = 0.02;

    // Parameters
    double U;
    double V;
    double W;
    double S;

    double J_fi, J_so, J_si, J = 0, Istim;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt = 0, dv_dt = 0, ds_dt = 0, dw_dt = 0;
    
    double time = 0.0;
    double last_time = 0.0;
    int step = 0;

    double a = atof(argv[7]);   // Ca++
    double b = atof(argv[8]);   // K+

    // Stimulation parameters
    double stimduration = 2.0;
    double stimstrength = 2.0;
    double tbegin = 0.0;
    double tend = tbegin + stimduration;
    int stim_count = 1;
    int stim_count2 = 0;

    // Files to save V and time
    FILE *fp_v, *fp_t, *fp_last, *fp_last_ap = fopen("last-ap.txt", "w"), *fp_last_ap_time = fopen("last-ap-time.txt", "w");
    fp_v = fopen("V.txt", "w");
    fp_t = fopen("t.txt", "w");

    // Initial conditions
    U = 0;
    V = 1;
    W = 1;
    S = 0;

    // Forward Euler
    for (step = 0; time < simulation_time; step++)
    {
        // Save V and time
        fprintf(fp_v, "%lf\n", 85.7*U-84);
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

        tau_vminus = (1 - H(U, theta_vminus)) * tau_v1minus + H(U, theta_vminus) * tau_v2minus;
        tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1 + tanh(k_wminus * (U - u_wminus))) / 2;
        tau_so = (tau_so1 + (tau_so2 - tau_so1) * (1 + tanh(k_so * (U - u_so))) / 2) * b;
        tau_s = (1 - H(U, theta_w)) * tau_s1 + H(U, theta_w) * tau_s2;
        tau_o = (1 - H(U, theta_o)) * tau_o1 + H(U, theta_o) * tau_o2;

        J_fi = -V * H(U, theta_v) * (U - theta_v) * (u_u - U) / tau_fi;                 // fast inward: Na+
        J_so = (U - u_o) * ((1 - H(U, theta_w)) / tau_o) + (H(U, theta_w) / tau_so);    // slow outward: K+
        J_si = -H(U, theta_w) * W * S / (tau_si*a);                                         // slow inward: Ca++
        
        J = J_fi + J_so + J_si;

        v_inf = v_inf_function(U, theta_vminus);
        w_inf = (1 - H(U, theta_o)) * (1 - U / tau_winf) + H(U, theta_o) * w_infstar;

        du_dt = -J + Istim;
        dv_dt = (1 - H(U, theta_v)) * (v_inf - V) / tau_vminus - H(U, theta_v) * V / tau_vplus;
        dw_dt = (1 - H(U, theta_w)) * (w_inf - W) / tau_wminus - H(U, theta_w) * W / tau_wplus;
        ds_dt = ((1 + tanh(k_s * (U - u_s))) / 2 - S) / tau_s;

        // Update variables
        U = U + du_dt * dt;
        V = V + dv_dt * dt;
        W = W + dw_dt * dt;
        S = S + ds_dt * dt;

        time += dt;

        if((stim_count2 == stim_count_max2 || stim_count2 == stim_count_max2 - 1) && last_time <= 1000)
        {
            fprintf(fp_last_ap, "%lf\n", 85.7*U-84);
            fprintf(fp_last_ap_time, "%lf\n", last_time);
            last_time += dt;
        }
    }

    fclose(fp_last_ap);
    fclose(fp_last_ap_time);

    return 0;
}