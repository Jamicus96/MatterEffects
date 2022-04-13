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
#include <globes/globes.h>


double Oscillation_Prob(std::vector<double> consts, double L, double E, double rho,
                        int init_flavour, int final_flavour, int anti);
double Oscillation_Prob_Vac(double m21, double m31, double PMNS_values[18], double L, double E,
                            int init_flavour, int final_flavour, int anti);
std::vector<double> compute_constants(double m21, double m31, double PMNS_values[18], int init_flavour,int final_flavour);
double* read_data(std::string filename, int num);

/* ---------------------------- */
/* ----- Global variables ----- */
/* ---------------------------- */

// Universal constants
double hbar = 6.58211957e-16;   // (eV.s)
double c = 299792458;           // (m/s)

// Convert matter density (g/cm^3) to matter potential (eV)
// by multiplying by these factors (from GLOBES-3.0.11/src/glb_probability.h):
double GLB_V_FACTOR = 7.5e-14;   /* Conversion factor for matter potentials */
double GLB_Ne_MANTLE = 0.5;     /* Effective electron numbers for calculation */


/* ---------------------------- */
/* ------ Main Functions ------ */
/* ---------------------------- */

int main(int argc, char *argv[]) {
    // Read in arguments
    double E = atof(argv[1]); // MeV
    double L_min = atof(argv[2]); // km
    double L_max = atof(argv[3]); // km
    int N = atof(argv[4]); // number of data points between L_min and L_max
    double rho = atof(argv[5]); // g/cm^3
    int init_flavour = atoi(argv[6]); // 0=e, 1=mu, 2=tau
    int final_flavour = atoi(argv[7]); // 0=e, 1=mu, 2=tau
    int anti = atoi(argv[8]); // -1 = antineutrino, 1 = neutrino

    // Read in calculated constants from files
    double* osc_consts = read_data("oscillation_constants.txt", 6);
    double* PMNS_values = read_data("PMNS.txt", 18);

    // Unpack osc_consts
    double theta12 = osc_consts[0];
    double theta13 = osc_consts[1];
    double theta23 = osc_consts[2];
    double delta = osc_consts[3];
    double m21 = osc_consts[4];
    double m31 = osc_consts[5];

    // Initialse my algorithm with pre-loop constants
    std::vector<double> consts = compute_constants(m21, m31, PMNS_values, init_flavour,final_flavour);

    // Initialise GLoBES
    glbInit(argv[0]); /* Initialize GLoBES library */
    glb_params true_values = glbAllocParams();
    glbDefineParams(true_values, theta12, theta13, theta23, delta, m21, m31);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // Open file to print results to
    std::ofstream datafile;
    datafile.open("results.txt", std::ios::app);

    // Compute oscillation probabilities via various methods, for various baselines
    if (N < 0) {
        std::cout << "Number of datapoints N must be at least 0, not " << N << std::cout;
    }
    double L_step = 0.0;
    if (N > 0) {L_step = (L_max - L_min) / N;}
    double L = L_min;
    for (unsigned int i = 0; i < N+1; ++i) {
        double P = Oscillation_Prob(consts, L, E, rho, init_flavour, final_flavour, anti);
        double P_vac = Oscillation_Prob_Vac(m21, m31, PMNS_values, L, E, init_flavour, final_flavour, anti);
        double P_globes = glbConstantDensityProbability(init_flavour + 1, final_flavour + 1, anti, E * 1e-3, L, rho);   // flavours + 1 to mine, and energy in GeV
        double P_vac_globes = glbVacuumProbability(init_flavour + 1, final_flavour + 1, anti, E * 1e-3, L);             // flavours + 1 to mine, and energy in GeV

        // std::cout << "P = " << P << ", P_vac = " << P_vac << ", P_globes = " << P_globes << std::endl;

        // Print results to file
        datafile << E << " " << rho << " " << L << " " << P << " " << P_vac << " " << P_globes << " " << P_vac_globes << std::endl;
        // Step baseline forward
        L += L_step;
    }

    return 0;
}



/**
 * @brief Read data from text file, with one number per line, into array of given length
 * 
 * @param filename 
 * @param num 
 * @return double*
 */
double* read_data(std::string filename, int num) {
    //Create an input file stream
    std::ifstream infile(filename.c_str());

    double* array = new double[num];
    double a;
    int i = 0;
    while (infile >> a) {
        array[i] = a;
        ++i;
    }
    
    return array;
}


/* -------------------- PRE-LOOP CONSTANTS ------------------- */

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
std::vector<double> compute_constants(double m21, double m31, double PMNS_values[18], int init_flavour,int final_flavour) {
    // initialise vector
    std::vector<double> vals;

    // scalar values
    vals.push_back(- (2.0/27.0) * (m21*m21*m21 + m31*m31*m31) + (1.0/9.0) * (m21*m21 * m31 + m21 * m31*m31)); // a0
    vals.push_back((1.0/3.0) * (m21 * m31 - m21*m21 - m31*m31));  // a1

    /* ------------------------ */

    // Unpack data in PMNS matrix (real and imginary parts)
    double U_r[3][3];
    double U_i[3][3];
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            U_r[i][j] = PMNS_values[3 * i + j];
            U_i[i][j] = PMNS_values[3 * i + j + 9];
        }
    }

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
        H_ee +=                 mf1[f] * (U_r[0][f]*U_r[0][f] + U_i[0][f]*U_i[0][f] - (1.0/3.0));
        Y_ee += (1.0/3.0) *     mf1_2[f] * (U_r[0][f]*U_r[0][f] + U_i[0][f]*U_i[0][f] - (1.0/3.0));
        H_r +=                  mf1[f] * (U_r[init_flavour][f]*U_r[final_flavour][f] + U_i[init_flavour][f]*U_i[final_flavour][f]);
        Y_r += (1.0/3.0) *      mf1_2[f] * (U_r[init_flavour][f]*U_r[final_flavour][f] + U_i[init_flavour][f]*U_i[final_flavour][f]);
        if (init_flavour == final_flavour) {
            H_r -= (1.0/3.0) *  mf1[f];
            Y_r -= (1.0/9.0) *  mf1_2[f];
        } else {
            H_i +=              mf1[f] * (U_i[init_flavour][f]*U_r[final_flavour][f] - U_r[init_flavour][f]*U_i[final_flavour][f]);
            Y_i += (1.0/3.0) *  mf1_2[f] * (U_i[init_flavour][f]*U_r[final_flavour][f] - U_r[init_flavour][f]*U_i[final_flavour][f]);
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
        T_r = (1.0/3,0) * H_r;
        T_i = (1.0/3,0) * H_i;
        if ((init_flavour * final_flavour) != 0) {
            T_r *= -(2.0/3.0);
            T_i *= -(2.0/3.0);
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
 * @param E (anti)neutrino energy (MeV)
 * @param rho Matter density (g/cm^3). Use same conversion to matter potential and GLoBES.
 * @param init_flavour Initial neutrino flavour (0=e, 1=mu, 2=tau).
 * @param final_flavour Final neutrino flavour (0=e, 1=mu, 2=tau).
 * @param anti 1=neutrino, -1=antineutrino.
 * @return double 
 */
double Oscillation_Prob(std::vector<double> consts, double L, double E, double rho,
                        int init_flavour, int final_flavour, int anti) {

    // Check anti value
    if (anti != 1 && anti != -1) {
        std::cout << "anti should be 1 or -1, not " << anti << std::endl;
        exit(1);
    }

    // convert units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)

    // Unpack constants {a0, a1, H_ee, Y_ee, H_r, H_i, Y_r, Y_i, D, T_r, T_i}
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
        // Convert matter density (g/cm^3) to matter potential (eV): V = GLB_V_FACTOR * GLB_Ne_MANTLE * rho
        // And: A_CC = +/- 2 * E * V
        double A_CC = anti * 2.0 * E * GLB_V_FACTOR * GLB_Ne_MANTLE * rho; // (eV^2)

        // Correct constants for matter effects (D = 0 for transition prob)
        a0 += -Y_ee * A_CC - (1.0/3.0) * H_ee * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
        a1 += -H_ee * A_CC - (1.0/3.0) * A_CC*A_CC;
        Y_r += A_CC * T_r + (1.0/3.0) * A_CC*A_CC * D;
        Y_i += A_CC * T_i;
        H_r += A_CC * D;
    }

    // Different cases for survival and transition probabilities
    double P;
    double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
    double preFact = 2.0 * sqrt(- a1 / 3.0);

    if (init_flavour == final_flavour) {
        // Get eigenvalues of H, and constants X
        double eigen[3];
        double X[3];
        for(int i=0; i<3; ++i){
            eigen[i] = preFact * cos(arcCos - (2.0 * M_PI * i) / 3.0);
            X[i] = (1.0/3.0) + (eigen[i] * H_r + Y_r) / (3.0 * eigen[i]*eigen[i] + a1);
            // std::cout << "(X_" << i << ")_" << init_flavour << final_flavour << " = " <<  X[i] << std::endl;
        }

        double s2_10 = sin(((eigen[1] - eigen[0]) * L) / (4.0 * E));
        double s2_20 = sin(((eigen[2] - eigen[0]) * L) / (4.0 * E));
        double s2_21 = sin(((eigen[2] - eigen[1]) * L) / (4.0 * E));

        // std::cout << "sin(((eigen[1] - eigen[0]) * L) / (4 * E)) = " << s2_10 << std::endl;
        // std::cout << "sin(((eigen[2] - eigen[0]) * L) / (4 * E)) = " << s2_20 << std::endl;
        // std::cout << "sin(((eigen[2] - eigen[1]) * L) / (4 * E)) = " << s2_21 << std::endl;

        // Compute probability
        P = 1.0 - 4.0 * (X[1]*X[0]*s2_10*s2_10 + X[2]*X[0]*s2_20*s2_20 + X[2]*X[1]*s2_21*s2_21);

    } else {
        // Get eigenvalues of H, and constants X
        double eigen[3];
        double R_X[3];
        double I_X[3];
        for(int i=0; i<3; ++i){
            eigen[i] = preFact * cos(arcCos - (2.0 * M_PI * i) / 3.0);
            R_X[i] = (eigen[i] * H_r + Y_r) / (3.0 * eigen[i]*eigen[i] + a1);
            I_X[i] = (eigen[i] * H_i + Y_i) / (3.0 * eigen[i]*eigen[i] + a1);
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
        P  =        - 4.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * s2_10*s2_10
                           + (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * s2_20*s2_20
                           + (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * s2_21*s2_21)
             + anti * 2.0 * ((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * s_10
                           + (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * s_20
                           + (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * s_21);

    }

    return P;
}


/**
 * @brief Vacuum oscillation calculation (standard, using PMNS matric elements).
 * 
 * @param m21 Mass difference (eV^2).
 * @param m31 Mass difference (eV^2).
 * @param PMNS_values Vector containing PMNS matrix elements.
 * @param L 
 * @param E 
 * @param init_flavour Initial neutrino flavour (0=e, 1=mu, 2=tau).
 * @param final_flavour Final neutrino flavour (0=e, 1=mu, 2=tau).
 * @param anti 1=neutrino, -1=antineutrino.
 * @return double 
 */
double Oscillation_Prob_Vac(double m21, double m31, double PMNS_values[18], double L, double E,
                            int init_flavour, int final_flavour, int anti) {
    // convert all units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)


    // Unpack data in PMNS matrix (real and imginary parts)
    double U_r[3][3];
    double U_i[3][3];
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            U_r[i][j] = PMNS_values[3 * i + j];
            U_i[i][j] = PMNS_values[3 * i + j + 9];
        }
    }

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
            U_mag2[i] = U_r[init_flavour][i]*U_r[init_flavour][i] + U_i[init_flavour][i]*U_i[init_flavour][i];

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
            U_Re[k] =  U_r[init_flavour][k] * U_r[final_flavour][k] + U_i[init_flavour][k] * U_i[final_flavour][k];
            U_Im[k] =  U_i[init_flavour][k] * U_r[final_flavour][k] - U_r[init_flavour][k] * U_i[final_flavour][k];
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