/**
 * @file Matter_Oscillations.cpp
 * @author James Page (jp643@sussex.ac.uk)
 * @brief Computing (anti)neutrino oscillation probabilities with matter effects,
 * assuming constant electron density.
 * @version 1
 * @date 2021-05-29
 * 
 * 
 */

//include
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <globes/globes.h>


// General flavour functions
double Oscillation_Prob(const std::vector<double>& consts, double L, double E, double rho,
                        int init_flavour, int final_flavour, int anti);
double Oscillation_Prob_Vac(const std::vector<std::vector<std::vector<double> > >& U_PMNS, double L, double E,
                            int init_flavour, int final_flavour, int anti);
std::vector<std::vector<std::vector<double> > > calculate_PMNS();
std::vector<double> compute_constants(const std::vector<std::vector<std::vector<double> > >& U_PMNS, int init_flavour,int final_flavour);

// Specific flavour functions
std::vector<double> e_e_Survival_Prob_Constants();
std::vector<double> mu_mu_Survival_Prob_Constants();
std::vector<double> mu_e_Transition_Prob_Constants();
double anti_e_e_Survival_Prob(const std::vector<double>& consts, double rho, double E, double L);
double mu_mu_Survival_Prob(const std::vector<double>& consts, double rho, double E, double L);
double mu_e_Transition_Prob(const std::vector<double>& consts, double rho, double E, double L);

                /* ---------------------------- */
                /* ----- Global variables ----- */
                /* ---------------------------- */

// Conversion from km to eV^-1 (from hbar = 6.58211957e-16 eV.s, c = 299792458 m/s)
// double L_FACTOR = 5.06773e+09;

// Convert matter density (g/cm^3) to matter potential (eV)
// by multiplying by these factors (from GLOBES-3.0.11/src/glb_probability.h):
const double GLB_V_FACTOR_ = 7.5e-14;   /* Conversion factor for matter potentials */
const double GLB_Ne_MANTLE_ = 0.5;     /* Effective electron numbers for calculation */
// eV to km conversion (also from GLoBES: GLOBES-3.0.11/globes/globes.h):
const double GLB_EV_TO_KM_FACTOR_ = 1.9747235e-10;

// Oscillation constants
const double theta12 = std::asin(sqrt(0.307));
const double theta13 = std::asin(sqrt(2.18e-2));
const double theta23 = std::asin(sqrt(0.545));      // Normal Hierarchy
// const double theta23 = std::asin(sqrt(0.547));      // Inverted Hierarchy
const double delta = 1.36 * M_PI;                   // Not sure about the hierarchy (NH I think here)
const double m21 = 7.53e-5;                         // (Delta m_21^2, in eV^2)
const double m31 = 0.0025283;                       // (Delta m_31^2, in eV^2)


                /* ---------------------------- */
                /* ------ Main Functions ------ */
                /* ---------------------------- */

int main(int argc, char *argv[]) {
    // Read in arguments
    std::vector<std::vector<double> > input_data = read_args(argc, argv, 16);
    int init_flavour = atoi(argv[13]); // 0=e, 1=mu, 2=tau
    int final_flavour = atoi(argv[14]); // 0=e, 1=mu, 2=tau
    int anti = atoi(argv[15]); // -1 = antineutrino, 1 = neutrino
    std::string output_file_address = argv[16];

    /* ---------------------------------------------------------- */

    // Compute PMNS matrix (kinda part of my initialisation, unless I use the mixing angles directly,
    // which I would definitely for one harcoded flavour transition)
    std::vector<std::vector<std::vector<double> > > U_PMNS = calculate_PMNS();

    // Initialise my algorithm with pre-loop constants
    std::vector<double> consts = compute_constants(U_PMNS, init_flavour,final_flavour);

    // Initialise flavour specific constants
    std::vector<double> consts_e_e = e_e_Survival_Prob_Constants();
    std::vector<double> consts_m_m = mu_mu_Survival_Prob_Constants();
    std::vector<double> consts_m_e = mu_e_Transition_Prob_Constants();

    // Initialise GLoBES
    glbInit(argv[0]); /* Initialize GLoBES library */
    glb_params true_values = glbAllocParams();
    glbDefineParams(true_values, theta12, theta13, theta23, delta, m21, m31);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* ---------------------------------------------------------- */

    // Set up data saving
    std::vector<double> P, P_vac, P_globes, P_vac_globes, P_specific;
    std::vector<long double> time_P, time_P_vac, time_P_globes, time_P_vac_globes, time_P_specific;
    clock_t c_start, c_end;
    unsigned int N = input_data.size();

    double P_temp, L, E, rho;
    for (unsigned int i = 0; i < N+1; ++i) {
        // Unpack input data
        L = input_data.at(i).at(0);
        E = input_data.at(i).at(1);
        rho = input_data.at(i).at(2);

        // Run oscillation function
        c_start = clock();
        P_temp = Oscillation_Prob(consts, L, E, rho, init_flavour, final_flavour, anti);
        c_end = clock();
        // Record info
        P.push_back(P_temp);
        time_P.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
        
        // Run oscillation function
        c_start = clock();
        P_temp = Oscillation_Prob_Vac(U_PMNS, L, E, init_flavour, final_flavour, anti);
        c_end = clock();
        // Record info
        P_vac.push_back(P_temp);
        time_P_vac.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
        
        // Run oscillation function
        c_start = clock();
        P_temp = glbConstantDensityProbability(init_flavour + 1, final_flavour + 1, anti, E, L, rho);  // flavours + 1 to mine, and energy in GeV
        c_end = clock();
        // Record info
        P_globes.push_back(P_temp);
        time_P_globes.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
        
        // Run oscillation function
        c_start = clock();
        P_temp = glbVacuumProbability(init_flavour + 1, final_flavour + 1, anti, E, L);  // flavours + 1 to mine, and energy in GeV
        c_end = clock();
        // Record info
        P_vac_globes.push_back(P_temp);
        time_P_vac_globes.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
    }

        /* ---------------------------------------------------------- */


    if (init_flavour == 0 && final_flavour == 0 && anti == -1) {
        for (unsigned int i = 0; i < N+1; ++i) {
            // Run oscillation function
            c_start = clock();
            P_temp = anti_e_e_Survival_Prob(consts_e_e, rho, E, L);
            c_end = clock();
            // Record info
            P_specific.push_back(P_temp);
            time_P_specific.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
        }
    } else if (init_flavour == 1 && anti == 1) {
        if (final_flavour == 1) {
            for (unsigned int i = 0; i < N+1; ++i) {
                // Run oscillation function
                c_start = clock();
                P_temp = mu_mu_Survival_Prob(consts_m_m, rho, E, L);
                c_end = clock();
                // Record info
                P_specific.push_back(P_temp);
                time_P_specific.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
            }
        } else if (final_flavour == 0) {
            for (unsigned int i = 0; i < N+1; ++i) {
                // Run oscillation function
                c_start = clock();
                P_temp = mu_e_Transition_Prob(consts_m_e, rho, E, L);
                c_end = clock();
                // Record info
                P_specific.push_back(P_temp);
                time_P_specific.push_back(1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);
            }
        }
    }

    /* ---------------------------------------------------------- */

    // Open file to print results to
    std::ofstream datafile;
    datafile.open(output_file_address, std::ios::app);

    // Print results to file
    for (unsigned int i = 0; i < N+1; ++i) {
        // Unpack input data
        L = input_data.at(i).at(0);
        E = input_data.at(i).at(1);
        rho = input_data.at(i).at(2);

        datafile << anti << init_flavour << final_flavour << " " << E << " " << rho << " " << L << " "
        << P.at(i) << " " << time_P.at(i) << " "
        << P_vac.at(i) << " " << time_P_vac.at(i) << " "
        << P_globes.at(i) << " " << time_P_globes.at(i) << " "
        << P_vac_globes.at(i) << " " << time_P_vac_globes.at(i);
                 
        if (P_specific.size() > 0) {
            datafile << " " << P_specific.at(i) << " " << time_P_specific.at(i);
        }
        datafile << std::endl;
    }

    return 0;
}


/* -------------------- MISCELLANEOUS --------------------- */

/**
 * @brief Create list of input data to use in oscillation functions
 * 
 * @param argc 
 * @param argv 
 * @return std::vector<std::vector<double> > = {{L_0, E_0, rho_0}, {L_1, E_1, rho_1}, ..., {L_n, E_n, rho_n}}
 */
std::vector<std::vector<double> > read_args(int argc, char *argv[], unsigned int n_args_wanted) {
    // Read in arguments
    if (argc != n_args_wanted) {
        std::cout << "Need " << n_args_wanted << " arguments, not " << argc << std::endl;
        exit(1);
    }
    
    double L_min = atof(argv[1]); // km
    double L_max = atof(argv[2]); // km
    int L_N = atoi(argv[3]); // number of data points between L_min and L_max
    bool L_log = atoi(argv[4]);
    double E_min = atof(argv[5]) * 1e-3; // convert MeV to GeV (for consistency with GLoBES functions)
    double E_max = atof(argv[6]) * 1e-3; // convert MeV to GeV (for consistency with GLoBES functions)
    int E_N = atoi(argv[7]); // number of data points between E_min and E_max
    bool E_log = atoi(argv[8]);
    double rho_min = atof(argv[9]); // g/cm^3
    double rho_max = atof(argv[10]); // g/cm^3
    int rho_N = atoi(argv[11]); // number of data points between rho_min and rho_max
    bool rho_log = atoi(argv[12]);

    // Set up input data
    if (L_N < 0 || E_N < 0 || rho_N < 0) {
        std::cout << "Number of datapoints E_N, rho_N and L_N must be at least 0, not " << E_N << ", " << rho_N << ", " << L_N << std::endl;
        exit(1);
    }

    std::vector<double> L_steps, E_steps, rho_steps;
    double L_n_max = log(L_max - L_min + 1.0);
    double E_n_max = log(E_max - E_min + 1.0);
    double rho_n_max = log(rho_max - rho_min + 1.0);
    double n;
    for (unsigned int i = 0; i < L_N; ++i) {
        if (L_log) {
            n = i * (L_n_max / L_N);
            L_steps.push_back(L_min + pow(10.0, n) - 1.0);
        } else {
            L_steps.push_back(L_min + i * ((L_max - L_min) / L_N));
        }
    }
    for (unsigned int i = 0; i < E_N; ++i) {
        if (E_log) {
            E_steps.push_back(E_min + 0);
        } else {
            L_steps.push_back(E_min + i * ((E_max - E_min) / E_N));
        }
    }
    for (unsigned int i = 0; i < rho_N; ++i) {
        if (rho_log) {
            rho_steps.push_back(rho_min + 0);
        } else {
            rho_steps.push_back(rho_min + i * ((rho_max - rho_min) / rho_N));
        }
    }

    std::vector<std::vector<double> > input_data;
    double L = L_min;
    double E = E_min;
    double rho = rho_min;
    for (unsigned int i = 0; i < L_N; ++i) {
        L = L_steps.at(i);
        for (unsigned int j = 0; j < E_N; ++j) {
            E = E_steps.at(j);
            for (unsigned int k = 0; k < rho_N; ++k) {
                rho = rho_steps.at(k);
                std::vector<double> line_data = {L, E, rho};
                input_data.push_back(line_data);
            }
        }
    }

    return input_data;
}


/* -------------------- PRE-LOOP CONSTANTS ------------------- */

/**
 * @brief Computes PMNS matrix elements from mixing angles and CP phase (delta).
 * These are all global variables. Returns vector of two 3x3 vectors of vectors
 * representing the real part of the PMNS matrix, and the imaginary part respectovely.
 * 
 * @return std::vector<std::vector<std::vector<double>>> 
 */
std::vector<std::vector<std::vector<double> > > calculate_PMNS() {
    // Compute useful constants
    double s12 = std::sin(theta12);
    double s13 = std::sin(theta13);
    double s23 = std::sin(theta23);
    double c12 = std::cos(theta12);
    double c13 = std::cos(theta13);
    double c23 = std::cos(theta23);
    double c_delta = std::cos(delta);
    double s_delta = std::sin(delta);

    // Declare vector of 2 3x3 matrices filled with zeros
    std::vector<std::vector<std::vector<double> > > U (2, std::vector<std::vector<double> > (3, std::vector<double> (3, 0.0)));

    // Real part of PMNS matrix
    U.at(0).at(0).at(0) = c12 * c13;                               U.at(0).at(0).at(1) = s12 * c13;                               U.at(0).at(0).at(2) = s13 * c_delta;
    U.at(0).at(1).at(0) = - s12 * c23 - s13 * s23 * c12 * c_delta; U.at(0).at(1).at(1) = c12 * c23 - s12 * s13 * s23 * c_delta;   U.at(0).at(1).at(2) = s23 * c13;
    U.at(0).at(2).at(0) = s12 * s23 - s13 * c12 * c23 * c_delta;   U.at(0).at(2).at(1) = - s23 * c12 - s12 * s13 * c23 * c_delta; U.at(0).at(2).at(2) = c13 * c23;

    // Imaginary part of PMNS matrix
    U.at(1).at(0).at(0) = 0.0;                                     U.at(1).at(0).at(1) = 0.0;                                     U.at(1).at(0).at(2) = - s13 * s_delta;
    U.at(1).at(1).at(0) = - s13 * s23 * c12 * s_delta;             U.at(1).at(1).at(1) = - s12 * s13 * s23 * s_delta;             U.at(1).at(1).at(2) = 0.0;
    U.at(1).at(2).at(0) = - s13 * c12 * c23 * s_delta;             U.at(1).at(2).at(1) = - s12 * s13 * c23 * s_delta;             U.at(1).at(2).at(2) = 0.0;

    return U;
}

/**
 * @brief Returns vector of relevent constants for oscillation calculation (D ant T are divided by 3 w.r.t my paper).
 * 
 * @param m21 
 * @param m31 
 * @param PMNS_values 
 * @param init_flavour 
 * @param final_flavour 
 * @return std::vector<double> = {a0, a1, H_ee, Y_ee, H_r, H_i, Y_r, Y_i, D, T_r, T_i}
 */
std::vector<double> compute_constants(const std::vector<std::vector<std::vector<double> > >& U_PMNS, int init_flavour,int final_flavour) {
    // initialise vector
    std::vector<double> vals;

    // scalar values
    vals.push_back((1.0/27.0) * (m21*m21*m21 + m31*m31*m31) - (1.0/18.0) * (m21*m21 * m31 + m21 * m31*m31)); // a0
    vals.push_back((1.0/9.0) * (m21*m21 + m31*m31 - m21 * m31));  // a1

    /* ------------------------ */

    // Unpack data in PMNS matrix (real and imginary parts)
    std::vector<std::vector<double> > U_r = U_PMNS.at(0);
    std::vector<std::vector<double> > U_i = U_PMNS.at(1);

    // Compute relevent matrix components
    double H_ee = 0.0;
    double Y_ee = 0.0;
    double H_r = 0.0;
    double H_i = 0.0;
    double Y_r = 0.0;
    double Y_i = 0.0;
    double mf1[3] = {0.0, m21, m31};
    double mf1_2[3] = {2.0*m21*m31, m21*m21, m31*m31};
    for (unsigned int f = 0; f < 3; ++f) {
        H_ee +=                 mf1[f] * (U_r.at(0).at(f)*U_r.at(0).at(f) + U_i.at(0).at(f)*U_i.at(0).at(f) - (1.0/3.0));
        Y_ee += (1.0/3.0) *     mf1_2[f] * (U_r.at(0).at(f)*U_r.at(0).at(f) + U_i.at(0).at(f)*U_i.at(0).at(f) - (1.0/3.0));
        H_r +=                  mf1[f] * (U_r.at(final_flavour).at(f)*U_r.at(init_flavour).at(f) + U_i.at(final_flavour).at(f)*U_i.at(init_flavour).at(f));
        Y_r += (1.0/3.0) *      mf1_2[f] * (U_r.at(final_flavour).at(f)*U_r.at(init_flavour).at(f) + U_i.at(final_flavour).at(f)*U_i.at(init_flavour).at(f));
        if (init_flavour == final_flavour) {
            H_r -= (1.0/3.0) *  mf1[f];
            Y_r -= (1.0/9.0) *  mf1_2[f];
        } else {
            H_i +=              mf1[f] * (U_i.at(final_flavour).at(f)*U_r.at(init_flavour).at(f) - U_r.at(final_flavour).at(f)*U_i.at(init_flavour).at(f));
            Y_i += (1.0/3.0) *  mf1_2[f] * (U_i.at(final_flavour).at(f)*U_r.at(init_flavour).at(f) - U_r.at(final_flavour).at(f)*U_i.at(init_flavour).at(f));
        }
    }
    vals.push_back(H_ee);   vals.push_back(Y_ee);   vals.push_back(H_r);
    vals.push_back(H_i);    vals.push_back(Y_r);    vals.push_back(Y_i);

    /* ---------------- */

    // Relevent components of matrix D (diagonal matrix)
    double D[3] = {2.0/3.0, -1.0/3.0, -1.0/3.0};
    if (init_flavour != final_flavour) {
        vals.push_back(0.0);
    } else {
        vals.push_back(D[init_flavour]);
    }

    // Compute relevent T matrix component
    double T_r = 0.0;
    double T_i = 0.0;
    if (init_flavour == final_flavour) {
        if (init_flavour == 0) {
            T_r = (2.0/3.0) * H_ee;
        } else {
            T_r = -(2.0/3.0) * (H_ee + H_r);
        }
    } else {
        T_r = (1.0/3.0) * H_r;
        T_i = (1.0/3.0) * H_i;
        if ((init_flavour * final_flavour) != 0) {
            T_r *= -2.0;
            T_i *= -2.0;
        }
    }
    vals.push_back(T_r);    vals.push_back(T_i);

    return vals;
}


/* -------------------- OSCILLATION PROBS ----------------------- */


/**
 * @brief Compute oscillation probability with my algorithm.
 * 
 * @param consts Pre-computed constants (only depend on oscillation parameters and particular flavour transition).
 * @param L Baseline (km)
 * @param E (anti)neutrino energy (GeV)
 * @param rho Matter density (g/cm^3). Use same conversion to matter potential and GLoBES.
 * @param init_flavour Initial neutrino flavour (0=e, 1=mu, 2=tau).
 * @param final_flavour Final neutrino flavour (0=e, 1=mu, 2=tau).
 * @param anti 1=neutrino, -1=antineutrino.
 * @return double 
 */
double Oscillation_Prob(const std::vector<double>& consts, double L, double E, double rho,
                        int init_flavour, int final_flavour, int anti) {

    // Check anti value
    if (anti != 1 && anti != -1) {
        std::cout << "anti should be 1 or -1, not " << anti << std::endl;
        exit(1);
    }

    // convert units to eV
    E *= 1e9; //(GeV to eV)
    L /= GLB_EV_TO_KM_FACTOR_; //(km to eV^-1)

    // Unpack constants {-a0, -a1, H_ee, Y_ee, H_r, H_i, Y_r, Y_i, D, T_r, T_i}
    double a0 = consts.at(0);
    double a1 = consts.at(1);
    double H_ee = consts.at(2);
    double Y_ee = consts.at(3);
    double H_r = consts.at(4);
    double H_i = consts.at(5);
    double Y_r = consts.at(6);
    double Y_i = consts.at(7);
    double D = consts.at(8);
    double T_r = consts.at(9);
    double T_i = consts.at(10);

    if (rho != 0.0) {
        // Convert matter density (g/cm^3) to matter potential (eV): V = GLB_V_FACTOR_ * GLB_Ne_MANTLE_ * rho
        // And: A_CC = +/- 2 * E * V
        double A_CC = anti * 2.0 * E * GLB_V_FACTOR_ * GLB_Ne_MANTLE_ * rho; // (eV^2)

        // Correct constants for matter effects (D = 0 for transition prob)
        a0 += 0.5 * Y_ee * A_CC + (1.0/6.0) * H_ee * A_CC*A_CC + (1.0/27.0) * A_CC*A_CC*A_CC; // a0
        a1 += (1.0/3.0) * H_ee * A_CC + (1.0/9.0) * A_CC*A_CC; // a1
        Y_r += A_CC * T_r;

        if (init_flavour == final_flavour) {
            Y_r += (1.0/3.0) * A_CC*A_CC * D;
            H_r += A_CC * D;
        } else {
            Y_i += A_CC * T_i;
        }
    }

    // Different cases for survival and transition probabilities
    double P;
    double arcCos = (1.0/3.0) * acos(a0 / (sqrt(a1) * a1));
    double preFact = 2.0 * sqrt(a1);
    double PHASE = (2.0/3.0) * M_PI;

    if (init_flavour == final_flavour) {
        // Get eigenvalues of H, and constants X
        double eigen[3];
        double X[3];
        for(int i=0; i<3; ++i){
            eigen[i] = preFact * cos(arcCos - PHASE * i);
            X[i] = 1.0 + (eigen[i] * H_r + Y_r) / (eigen[i]*eigen[i] - a1);
        }

        double s2_10 = sin(((eigen[1] - eigen[0]) * L) / (4.0 * E));
        double s2_20 = sin(((eigen[2] - eigen[0]) * L) / (4.0 * E));
        double s2_21 = sin(((eigen[2] - eigen[1]) * L) / (4.0 * E));

        // Compute probability
        P = 1.0 - (4.0/9.0) * (X[1]*X[0]*s2_10*s2_10 + X[2]*X[0]*s2_20*s2_20 + X[2]*X[1]*s2_21*s2_21);

    } else {
        // Get eigenvalues of H, and constants X
        double eigen[3];
        double R_X[3];
        double I_X[3];
        for(int i=0; i<3; ++i){
            eigen[i] = preFact * cos(arcCos - PHASE * i);
            R_X[i] = (eigen[i] * H_r + Y_r) / (eigen[i]*eigen[i] - a1);
            I_X[i] = (eigen[i] * H_i + Y_i) / (eigen[i]*eigen[i] - a1);
        }

        double Theta_10 = ((eigen[1] - eigen[0]) * L) / (2.0 * E);
        double Theta_20 = ((eigen[2] - eigen[0]) * L) / (2.0 * E);
        double Theta_21 = ((eigen[2] - eigen[1]) * L) / (2.0 * E);

        double s2_10 = sin(Theta_10 / 2.0);
        double s2_20 = sin(Theta_20 / 2.0);
        double s2_21 = sin(Theta_21 / 2.0);
        double s_10 = sin(Theta_10);
        double s_20 = sin(Theta_20);
        double s_21 = sin(Theta_21);

        // Compute probability
        P  = (2.0/9.0) * (- 2.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * s2_10*s2_10
                                 + (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * s2_20*s2_20
                                 + (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * s2_21*s2_21)
                         + anti * ((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * s_10
                                 + (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * s_20
                                 + (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * s_21));

    }

    return P;
}

/**
 * @brief Vacuum oscillation calculation (standard, using PMNS matric elements).
 * 
 * @param m21 Mass difference (eV^2).
 * @param m31 Mass difference (eV^2).
 * @param PMNS_values Vector containing PMNS matrix elements.
 * @param L Baseline (km)
 * @param E (anti)neutrino energy (GeV)
 * @param init_flavour Initial neutrino flavour (0=e, 1=mu, 2=tau).
 * @param final_flavour Final neutrino flavour (0=e, 1=mu, 2=tau).
 * @param anti 1=neutrino, -1=antineutrino.
 * @return double 
 */
double Oscillation_Prob_Vac(const std::vector<std::vector<std::vector<double> > >& U_PMNS, double L, double E,
                            int init_flavour, int final_flavour, int anti) {
    // convert all units to eV
    E *= 1e9; //(GeV to eV)
    L /= GLB_EV_TO_KM_FACTOR_; //(km to eV^-1)


    // Unpack data in PMNS matrix (real and imginary parts)
    std::vector<std::vector<double> > U_r = U_PMNS.at(0);
    std::vector<std::vector<double> > U_i = U_PMNS.at(1);

    // Compute oscillation probability: different case for survival and transition
    double P;
    double m32 = m31 - m21;

    if (init_flavour == final_flavour) {
        double s2_21 = sin((m21 * L) / (4.0 * E));
        double s2_31 = sin((m31 * L) / (4.0 * E));
        double s2_32 = sin((m32 * L) / (4.0 * E));

        // std::cout << "sin((m21 * L) / (4 * E)) = " << s2_21 << std::endl;
        // std::cout << "sin((m31 * L) / (4 * E)) = " << s2_31 << std::endl;
        // std::cout << "sin((m32 * L) / (4 * E)) = " << s2_32 << std::endl;

        double U_mag2[3];
        for(int i=0; i<3; ++i){
            U_mag2[i] = U_r.at(init_flavour).at(i)*U_r.at(init_flavour).at(i) + U_i.at(init_flavour).at(i)*U_i.at(init_flavour).at(i);

            // std::cout << "|U_" << init_flavour << i << "| = " <<  U_mag2[i] << std::endl;
        }

        // Compute probability
        P = 1.0 - 4.0 * (U_mag2[1] * U_mag2[0] * s2_21*s2_21
                       + U_mag2[2] * U_mag2[0] * s2_31*s2_31
                       + U_mag2[2] * U_mag2[1] * s2_32*s2_32);

    } else {

        double Theta_21 = (m21 * L) / (2.0 * E);
        double Theta_31 = (m31 * L) / (2.0 * E);
        double Theta_32 = (m32 * L) / (2.0 * E);

        double s2_21 = sin(Theta_21 / 2.0);
        double s2_31 = sin(Theta_31 / 2.0);
        double s2_32 = sin(Theta_32 / 2.0);
        double s_21 = sin(Theta_21);
        double s_31 = sin(Theta_31);
        double s_32 = sin(Theta_32);

        // Get relevent PMNS matrix element products
        double U_Re[3];
        double U_Im[3];
        for (int k=0; k<3; ++k) {
            U_Re[k] =  U_r.at(final_flavour).at(k) * U_r.at(init_flavour).at(k) + U_i.at(final_flavour).at(k) * U_i.at(init_flavour).at(k);
            U_Im[k] =  U_i.at(final_flavour).at(k) * U_r.at(init_flavour).at(k) - U_r.at(final_flavour).at(k) * U_i.at(init_flavour).at(k);
        }

        // Compute probability
        P  =       - 4.0 * ((U_Re[1]*U_Re[0] + U_Im[1]*U_Im[0]) * s2_21*s2_21
                          + (U_Re[2]*U_Re[0] + U_Im[2]*U_Im[0]) * s2_31*s2_31
                          + (U_Re[2]*U_Re[1] + U_Im[2]*U_Im[1]) * s2_32*s2_32)
            + anti * 2.0 * ((U_Im[1]*U_Re[0] - U_Re[1]*U_Im[0]) * s_21
                          + (U_Im[2]*U_Re[0] - U_Re[2]*U_Im[0]) * s_31
                          + (U_Im[2]*U_Re[1] - U_Re[2]*U_Im[1]) * s_32);

    }

    return P;
}




/* ---------------------------------------------------------------------------------- */
/* ----------------------- SPECIFIC FLAVOUR TRANSITIONS ----------------------------- */
/* ---------------------------------------------------------------------------------- */


/* ----------------- Constant functions -------------------- */

/**
 * @brief Computes constants needed for anti-electron neutrino survival probability function.
 * 
 * @return std::vector<double> = {a0, a1, H_ee, Y_ee}
 */
std::vector<double> e_e_Survival_Prob_Constants() {
    // Useful constants
    double s12 = std::sin(theta12);
    double s13 = std::sin(theta13);
    double s23 = std::sin(theta23);
    double c12 = std::cos(theta12);
    double c13 = std::cos(theta13);
    double c23 = std::cos(theta23);

    // Compute vacuum constants
    std::vector<double> consts;
    consts.push_back((1.0/27.0) * (m21*m21*m21 + m31*m31*m31) - (1.0/18.0) * (m21*m21 * m31 + m21 * m31*m31)); // a0
    consts.push_back((1.0/9.0) * (m21*m21 + m31*m31 - m21 * m31)); // a1
    consts.push_back(m21 * (s12*s12 * c13*c13 - (1.0/3.0)) + m31 * (s13*s13 - (1.0/3.0))); // H_ee
    consts.push_back((1.0/3.0) * (m21*m21 * (s12*s12 * c13*c13 - (1.0/3.0))
                                    + m31*m31 * (s13*s13 - (1.0/3.0))
                                    + 2.0 * m21 * m31 * (c12*c12 * c13*c13 - (1.0/3.0)))); // Y_ee

    return consts;
}

/**
 * @brief Computes constants needed for muon neutrino survival probability function.
 * 
 * @return std::vector<double> = {a0, a1, H_ee, Y_ee, H_mm, Y_mm}
 */
std::vector<double> mu_mu_Survival_Prob_Constants() {
    // Useful constants
    double s12 = std::sin(theta12);
    double s13 = std::sin(theta13);
    double s23 = std::sin(theta23);
    double c12 = std::cos(theta12);
    double c13 = std::cos(theta13);
    double c23 = std::cos(theta23);

    // Compute vacuum constants
    std::vector<double> consts;
    consts.push_back((1.0/27.0) * (m21*m21*m21 + m31*m31*m31) - (1.0/18.0) * (m21*m21 * m31 + m21 * m31*m31)); // a0
    consts.push_back((1.0/9.0) * (m21*m21 + m31*m31 - m21 * m31)); // a1
    consts.push_back(m21 * (s12*s12 * c13*c13 - (1.0/3.0)) + m31 * (s13*s13 - (1.0/3.0))); // H_ee
    consts.push_back((1.0/3.0) * (m21*m21 * (s12*s12 * c13*c13 - (1.0/3.0))
                                + m31*m31 * (s13*s13 - (1.0/3.0))
                                + 2.0 * m21 * m31 * (c12*c12 * c13*c13 - (1.0/3.0)))); // Y_ee

    /*~~~~~~~~ mu -> mu Extra stuff ~~~~~~~~~*/

    double c_delta = cos(delta);

    consts.push_back(m21 * (c12*c12 * c23*c23 + s12*s12 * s13*s13 * s23*s23 - 2.0 * s12 * s13 * s23 * c12 * c23 * c_delta - (1.0/3.0))
                   + m31 * (s23*s23 * c13*c13 - (1.0/3.0))); // H_mm
    consts.push_back((1.0/3.0) * (m21*m21 * (c12*c12 * c23*c23 + s12*s12 * s13*s13 * s23*s23 - 2.0 * s12 * s13 * s23 * c12 * c23 * c_delta - (1.0/3.0))
                                + m31*m31 * (s23*s23 * c13*c13 - (1.0/3.0))
                                + 2.0 * m21 * m31 * (s12*s12 * c23*c23 + s13*s13 * s23*s23 * c12*c12 + 2.0 * s12 * s13 * s23 * c12 * c23 * c_delta - (1.0/3.0)))); // Y_mm

    return consts;
}

/**
 * @brief Computes constants needed for muon neutrino to electron neutrino transition probability function.
 * 
 * @return std::vector<double> = {a0, a1, H_ee, Y_ee, R[H_em], I[H_em], R[Y_em], I[Y_em]}
 */
std::vector<double> mu_e_Transition_Prob_Constants() {
    // Useful constants
    double s12 = std::sin(theta12);
    double s13 = std::sin(theta13);
    double s23 = std::sin(theta23);
    double c12 = std::cos(theta12);
    double c13 = std::cos(theta13);
    double c23 = std::cos(theta23);

    // Compute vacuum constants
    std::vector<double> consts;
    consts.push_back((1.0/27.0) * (m21*m21*m21 + m31*m31*m31) - (1.0/18.0) * (m21*m21 * m31 + m21 * m31*m31)); // a0
    consts.push_back((1.0/9.0) * (m21*m21 + m31*m31 - m21 * m31)); // a1
    consts.push_back(m21 * (s12*s12 * c13*c13 - (1.0/3.0)) + m31 * (s13*s13 - (1.0/3.0))); // H_ee
    consts.push_back((1.0/3.0) * (m21*m21 * (s12*s12 * c13*c13 - (1.0/3.0))
                                    + m31*m31 * (s13*s13 - (1.0/3.0))
                                    + 2.0 * m21 * m31 * (c12*c12 * c13*c13 - (1.0/3.0)))); // Y_ee

    /*~~~~~~~~ mu -> e Extra stuff ~~~~~~~~~*/

    double c_delta = cos(delta);
    double s_delta = sin(delta);

    consts.push_back(m21 * s12 * c13 * (c12 * c23 - s12 * s23 * s13 * c_delta)
                    + m31 * s13 * s23 * c13 * c_delta); // R[H_em]
    consts.push_back(m21 * s12*s12 * s13 * s23 * c13 * s_delta - m31 * s13 * s23 * c13 * s_delta); // I[H_em]

    double R_X2 = (1.0/3.0) * s12 * c13 * (c12 * c23 - s12 * s13 * s23 * c_delta);
    double R_X3 = (1.0 / 3.0) * s13 * s23 * c13 * c_delta;
    double R_X23 = - (2.0/3.0) * c12 * c13 * (s12 * c23 + s13 * s23 * c12 * c_delta);
    consts.push_back(m21*m21 * R_X2 + m31*m31 * R_X3 + m21 * m31 * R_X23); // R[Y_em]

    double I_X2 = (1.0/3.0) * s12*s12 * s13 * s23 * c13 * s_delta;
    double I_X3 = - (1.0 / 3.0) * s13 * s23 * c13 * s_delta;
    double I_X23 = (2.0/3.0) * s13 * s23 * c12*c12 * c13 * s_delta; 
    consts.push_back(m21*m21 * I_X2 + m31*m31 * I_X3 + m21 * m31 * I_X23); // I[Y_em]

    return consts;
}


/* ------------------ Oscillations Functions ----------------- */

/**
 * @brief Compute survival probability of anti-electron neutrinos, with matter effects.
 * 
 * @param constants Precomputed, independent of neutrino and matter effects, packaged as:
 *   constants = = {a0, a1, H_ee, Y_ee}
 * @param rho Matter density (g/cm^3). Use same conversion to matter potential and GLoBES.
 * @param E Antineutrino energy (GeV).
 * @param L Baseline (km).
 * @param anti true=antineutrino, false=neutrino.
 * @return double 
 */
double anti_e_e_Survival_Prob(const std::vector<double>& consts, double rho, double E, double L) {
    // convert units to eV
    E *= 1e9; //(GeV to eV)
    L /= GLB_EV_TO_KM_FACTOR_; //(km to eV^-1)

    // Initialise constants with vacuum values
    double a0 = consts.at(0);
    double a1 = consts.at(1);
    double H_ee = consts.at(2);
    double Y_ee = consts.at(3);

    // If not vacuum, make corrections
    if(rho != 0.0){
        double A_CC = -2.0 * E * GLB_V_FACTOR_ * GLB_Ne_MANTLE_ * rho; // (eV^2)
        // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
        double alpha_1 = (1.0/3.0) * H_ee * A_CC + (1.0/9.0) * A_CC*A_CC;
        a0 += 0.5 * Y_ee * A_CC + (1.0/6.0) * H_ee * A_CC*A_CC + (1.0/27.0) * A_CC*A_CC*A_CC; // a0
        a1 += alpha_1; // a1
        Y_ee += 2.0 * alpha_1;
        H_ee += (2.0/3.0) * A_CC;
    }


    // Get eigenvalues of H, and constants X and theta
    double eigen[3];
    double X[3];

    double arcCos = (1.0/3.0) * acos(a0 / (sqrt(a1) * a1));
    double preFact = 2.0 * sqrt(a1);
    double PHASE = (2.0/3.0) * M_PI;

    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - PHASE * i);
        X[i] = 1.0 + (eigen[i] * H_ee + Y_ee) / (eigen[i]*eigen[i] - a1);
    }

    double s_10 = sin(((eigen[1] - eigen[0]) * L) / (4.0 * E));
    double s_20 = sin(((eigen[2] - eigen[0]) * L) / (4.0 * E));
    double s_21 = sin(((eigen[2] - eigen[1]) * L) / (4.0 * E));


    // Compute probability
    return 1.0 - (4.0/9.0) * (X[1]*X[0]*s_10*s_10 + X[2]*X[0]*s_20*s_20 + X[2]*X[1]*s_21*s_21);
}

/**
 * @brief Compute survival probability of muon neutrinos, with matter effects.
 * 
 * @param constants Precomputed, independent of neutrino and matter effects, packaged as:
 *   constants = {a0, a1, H_ee, Y_ee, H_mm, Y_mm}
 * @param rho Matter density (g/cm^3). Use same conversion to matter potential and GLoBES.
 * @param E Antineutrino energy (GeV).
 * @param L Baseline (km).
 * @param anti true=antineutrino, false=neutrino.
 * @return double 
 */
double mu_mu_Survival_Prob(const std::vector<double>& consts, double rho, double E, double L) {
    // convert units to eV
    E *= 1e9; //(GeV to eV)
    L /= GLB_EV_TO_KM_FACTOR_; //(km to eV^-1)

    // Initialise constants with vacuum values
    double a0 = consts.at(0);
    double a1 = consts.at(1);
    double H_mm = consts.at(4);
    double Y_mm = consts.at(5);

    // If not vacuum, make corrections
    if(rho != 0.0){
        double A_CC = 2.0 * E * GLB_V_FACTOR_ * GLB_Ne_MANTLE_ * rho; // (eV^2)
        // Compute new values for H_mm, Y_mm, a0 and a1
        a0 += 0.5 * consts.at(3) * A_CC + (1.0/6.0) * consts.at(2) * A_CC*A_CC + (1.0/27.0) * A_CC*A_CC*A_CC; // a0
        a1 += (1.0/3.0) * consts.at(2) * A_CC + (1.0/9.0) * A_CC*A_CC; // a1
        Y_mm -= (2.0/3.0) * (consts.at(2) + H_mm) * A_CC + (1.0/9.0) * A_CC*A_CC;
        H_mm -= (1.0/3.0) * A_CC;
    }


    // Get eigenvalues of H, and constants X and theta
    double eigen[3];
    double X[3];

    double arcCos = (1.0/3.0) * acos(a0 / (sqrt(a1) * a1));
    double preFact = 2.0 * sqrt(a1);
    double PHASE = (2.0/3.0) * M_PI;

    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - PHASE * i);
        X[i] = 1.0 + (eigen[i] * H_mm + Y_mm) / (eigen[i]*eigen[i] - a1);
    }

    double s_10 = sin(((eigen[1] - eigen[0]) * L) / (4.0 * E));
    double s_20 = sin(((eigen[2] - eigen[0]) * L) / (4.0 * E));
    double s_21 = sin(((eigen[2] - eigen[1]) * L) / (4.0 * E));


    // Compute probability
    return 1.0 - (4.0/9.0) * (X[1]*X[0]*s_10*s_10 + X[2]*X[0]*s_20*s_20 + X[2]*X[1]*s_21*s_21);
}

/**
 * @brief Compute transition probability of muon neutrinos to electron neutrinos,
 * with matter effects.
 * 
 * @param constants Precomputed, independent of neutrino and matter effects, packaged as:
 *   constants = {a0, a1, H_ee, Y_ee, R[H_em], I[H_em], R[Y_em], I[Y_em]}
 * @param rho Matter density (g/cm^3). Use same conversion to matter potential and GLoBES.
 * @param E Neutrino energy (GeV).
 * @param L Baseline (km).
 * @return double 
 */
double mu_e_Transition_Prob(const std::vector<double>& consts, double rho, double E, double L) {
    // convert units to eV
    E *= 1e9; //(GeV to eV)
    L /= GLB_EV_TO_KM_FACTOR_; //(km to eV^-1)

    // Initialise constants with vacuum values
    double a0 = consts.at(0);
    double a1 = consts.at(1);
    double R_H_em = consts.at(4);
    double I_H_em = consts.at(5);
    double R_Y_em = consts.at(6);
    double I_Y_em = consts.at(7);

    // If not vacuum, make corrections
    double A_CC = 0.0;
    if(rho != 0.0){
        double A_CC = 2.0 * E * GLB_V_FACTOR_ * GLB_Ne_MANTLE_ * rho; // (eV^2)
        // Compute new values for a0 and a1
        a0 += 0.5 * consts.at(3) * A_CC + (1.0/6.0) * consts.at(2) * A_CC*A_CC + (1.0/27.0) * A_CC*A_CC*A_CC; // a0
        a1 += (1.0/3.0) * consts.at(2) * A_CC + (1.0/9.0) * A_CC*A_CC; // a1
    }

    // Get eigenvalues of H, and constants X and theta
    double eigen[3];
    double R_X[3];
    double I_X[3];

    double arcCos = (1.0/3.0) * acos(a0 / (sqrt(a1) * a1));
    double preFact = 2.0 * sqrt(a1);
    double PHASE = (2.0/3.0) * M_PI;

    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - PHASE * i);
        R_X[i] = ((eigen[i] + (A_CC / 3.0)) * R_H_em + R_Y_em) / (eigen[i]*eigen[i] - a1);
        I_X[i] = ((eigen[i] + (A_CC / 3.0)) * I_H_em + I_Y_em) / (eigen[i]*eigen[i] - a1);
    }

    double Theta_10 = ((eigen[1] - eigen[0]) * L) / (2.0 * E);
    double Theta_20 = ((eigen[2] - eigen[0]) * L) / (2.0 * E);
    double Theta_21 = ((eigen[2] - eigen[1]) * L) / (2.0 * E);

    double s2_10 = sin(Theta_10 / 2.0);
    double s2_20 = sin(Theta_20 / 2.0);
    double s2_21 = sin(Theta_21 / 2.0);
    double s_10 = sin(Theta_10);
    double s_20 = sin(Theta_20);
    double s_21 = sin(Theta_21);

    // Compute probability
    return (2.0/9.0) * (- 2.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * s2_10*s2_10
                               + (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * s2_20*s2_20
                               + (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * s2_21*s2_21)
                              + ((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * s_10
                               + (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * s_20
                               + (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * s_21));
}
