/**
 * @file Matter_Oscillations.cpp
 * @author James Page (jp643@sussex.ac.uk)
 * @brief Computing anti-electron neutrino survival probability with matter effects,
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


double Survival_Prob(double constants[4], double Ne, double E, double L, bool anti=true);
double Survival_Prob_Vac(std::string filename, double E, double L);
double Survival_Prob_Globes(double Ne, double E, double L, bool anti=true);
double Survival_Prob_Vac_Globes(double E, double L, bool anti=true);
double Transition_Prob(double constants[8], double Ne, double E, double L, bool anti=true);
double Transition_Prob_Vac(std::string filename, double E, double L, bool anti=true);
double Transition_Prob_Globes(double Ne, double E, double L, bool anti=true);
double Transition_Prob_Vac_Globes(double E, double L, bool anti=true);
double* read_data(std::string filename, int num);

std::vector<std::vector<double>> compute_constants(double m21, double m31);
std::vector<std::vector<std::vector<double>>> compute_vac_matrices(double PMNS_values[18], double m21, double m31);
std::vector<std::vector<std::vector<double>>> compute_T_matrices(std::vector<std::vector<double>> H_r, std::vector<std::vector<double>> H_i);

// Global variables
double GF = 1.663788e-14;       // (eV^-2)
double hbar = 6.58211957e-16;   // (eV.s)
double c = 299792458;           // (m/s)



int main(int argc, char *argv[]) {
    // Read in calculated constants
    double* osc_consts = read_data("oscillation_constants.txt", 6);
    double* PMNS_values = read_data("PMNS.txt", 18);

    // Read in arguments
    double E = atof(argv[1]); // MeV
    double L = atof(argv[2]); // km
    double rho = atof(argv[3]); // g/cm^3
    int init_flavour = atoi(argv[4]); // 0=e, 1=mu, 2=tau
    int final_flavour = atoi(argv[5]); // 0=e, 1=mu, 2=tau
    int anti = atoi(argv[6]); // -1 = antineutrino, 1 = neutrino

    // Convert matter density of matter potential (from GLOBES-3.0.11/src/glb_probability.h):
    double GLB_V_FACTOR = 7.5e-14;   /* Conversion factor for matter potentials */
    double GLB_Ne_MANTLE = 0.5;     /* Effective electron numbers for calculation */
    double V = GLB_V_FACTOR * GLB_Ne_MANTLE * rho;  // Matter potential (eV)

    // Initialse my algorithm with pre-loop constants
    double m21 = osc_consts[4];
    double m31 = osc_consts[5];
    std::vector<std::vector<double>> consts = compute_constants(m21, m31);
    std::vector<double> a = consts.at(0);
    std::vector<double> D = consts.at(1);
    std::vector<std::vector<std::vector<double>>> vac_matrices = compute_vac_matrices(PMNS_values[18], m21, m31);
    std::vector<std::vector<double>> H_r = vac_matrices.at(0);
    std::vector<std::vector<double>> H_i = vac_matrices.at(1);
    std::vector<std::vector<double>> Y_r = vac_matrices.at(2);
    std::vector<std::vector<double>> Y_i = vac_matrices.at(3);
    std::vector<std::vector<std::vector<double>>> T = compute_T_matrices(vac_matrices.at(0), vac_matrices.at(1));

    double P;
    double P_vac;
    double P_globes;
    double P_vac_globes;
    glbInit(argv[0]); /* Initialize GLoBES library */

    if (type == 0) {
        P = Survival_Prob(survival_constants, Ne, E, L, true);
        P_vac = Survival_Prob_Vac("PMNS_m2_data.txt", E, L);
        P_globes = Survival_Prob_Globes(Ne, E, L, true);
        P_vac_globes = Survival_Prob_Vac_Globes(E, L, true);
    } else if (type == 1) {
        P = Transition_Prob(transition_constants, Ne, E, L, true);
        P_vac = Transition_Prob_Vac("PMNS_m2_data.txt", E, L, true);
        P_globes = Transition_Prob_Globes(Ne, E, L, true);
        P_vac_globes = Transition_Prob_Vac_Globes(E, L, true);
    }
    
    // std::cout << "P = " << P << ", P_vac = " << P_vac << ", P_globes = " << P_globes << std::endl;

    // Print results to file
    std::ofstream datafile;
    datafile.open("results.txt", std::ios::app);
    datafile << P << " " << P_vac << " " << P_globes << " " << P_vac_globes << std::endl;

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
    for (unsigned int f = 0; f < 3; ++f) {
        H_ee +=                 mf1.at(f) * (U_r[0][f]*U_r[0][f] + U_i[0][f]*U_i[0][f] - (1.0/3.0));
        Y_ee += (1.0/3.0) *     mf1_2.at(f) * (U_r[0][f]*U_r[0][f] + U_i[0][f]*U_i[0][f] - (1.0/3.0));
        H_r +=                  mf1.at(f) * (U_r[init_flavour][f]*U_r[final_flavour][f] + U_i[init_flavour][f]*U_i[final_flavour][f]);
        Y_r += (1.0/3.0) *      mf1_2.at(f) * (U_r[init_flavour][f]*U_r[final_flavour][f] + U_i[init_flavour][f]*U_i[final_flavour][f]);
        if (init_flavour == final_flavour) {
            H_r -= (1.0/3.0) *  mf1.at(f);
            Y_r -= (1.0/9.0) *  mf1_2.at(f);
        } else {
            H_i +=              mf1.at(f) * (U_i[init_flavour][f]*U_r[final_flavour][f] - U_r[init_flavour][f]*U_i[final_flavour][f]);
            Y_i += (1.0/3.0) *  mf1_3.at(f) * (U_i[init_flavour][f]*U_r[final_flavour][f] - U_r[init_flavour][f]*U_i[final_flavour][f]);
        }
    }
    vals.push_back(H_ee);   vals.push_back(Y_ee);   vals.push_back(H_r);
    vals.push_back(H_i);    vals.push_back(Y_r);    vals.push_back(Y_i);

    /* ---------------- */

    // Relevent components of matrix D (diagonal matrix)
    std::vector<double> D = {2.0/3.0, -1.0/3.0, -1.0/3.0};
    if (init_flavour != final_flavour) {
        vals.push_back(0.0);
    } else {
        vals.push_back(D.at(init_flavour));
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


/* -------------------- SURVIVAL PROBS ----------------------- */

double Oscillation_Prob(std::vector<std::vector<double>> consts, std::vector<std::vector<std::vector<double>>> vac_matrices,
                        std::vector<std::vector<std::vector<double>>> T, double L, double E, double V,
                        int init_flavour,int final_flavour, int anti) {

    // Unpack constants
    std::vector<double> a = consts.at(0);
    std::vector<double> D = consts.at(1);
    std::vector<std::vector<double>> H_r = vac_matrices.at(0);
    std::vector<std::vector<double>> H_i = vac_matrices.at(1);
    std::vector<std::vector<double>> Y_r = vac_matrices.at(2);
    std::vector<std::vector<double>> Y_i = vac_matrices.at(3);
    std::vector<std::vector<double>> T_r = T.at(0);
    std::vector<std::vector<double>> T_i = T.at(1);

    if (V != 0.0) {
        Y_r.at(final_flavour).at(init_flavour)
    }
}


/**
 * @brief Compute survival probability of (anti)electron neutrinos, with matter effects.
 * 
 * @param constants Precomputed, independent of neutrino and matter effects, packaged as:
 *   constants[0] << H_ee_vac,
 *   constants[1] << Y_vac,
 *   constants[2] << a0_vac,
 *   constants[3] << a1_vac,
 * @param Ne Matter electron density (m^-3).
 * @param E Neutrino energy (MeV).
 * @param L Propagation length (km).
 * @param anti true=antineutrino, false=neutrino.
 * @return double 
 */
double Survival_Prob(double constants[4], double Ne, double E, double L, bool anti) {
    // convert units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)

    // Initialise constants with vacuum values
    double H = constants[0];
    double Y = constants[1];
    double a0 = constants[2];
    double a1 = constants[3];

    // std::cout << "H_ee = " << H << std::endl;
    // std::cout << "Y_ee = " << Y << std::endl;
    // std::cout << "a0 = " << a0 << std::endl;
    // std::cout << "a1 = " << a1 << std::endl;

    // If not vacuum, make corrections
    if(Ne != 0.0){
        // convert units to eV
        //std::cout << "Ne = " << Ne << std::endl;
        Ne *= c*c*c * hbar*hbar*hbar; //(m^-3 to ev^3)
        //std::cout << "Ne = " << Ne << std::endl;

        // Work out ACC in (eV)^2
        double A_CC = 2 * sqrt(2) * E * GF * Ne;
        if(anti){
            A_CC *= -1;
        }
        //std::cout << "A_CC = " << A_CC << std::endl;

        // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
        double alpha_1 = -H * A_CC - (1.0/3.0) * A_CC*A_CC;

        a0 += -Y * A_CC - (1.0/3.0) * H * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
        a1 += -H * A_CC - (1.0/3.0) * A_CC*A_CC;
        Y -= (2.0/3.0) * alpha_1;
        H += (2.0/3.0) * A_CC;
    }


    // Get eigenvalues of H, and constants X and theta
    double eigen[3];
    double X[3];

    double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
    double preFact = 2.0 * sqrt(- a1 / 3.0);

    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - (2.0 * M_PI * i) / 3.0);
        //std::cout << "E_" << i << " = " << eigen[i] << std::endl;

        X[i] = (1.0/3.0) + (eigen[i] * H + Y) / (3.0 * eigen[i]*eigen[i] + a1);
        //std::cout << "X_" << i << " = " << X[i] << std::endl;
    }

    // std::cout << "E_10 = " << eigen[1] - eigen[0] << std::endl;
    // std::cout << "E_20 = " << eigen[2] - eigen[0] << std::endl;
    // std::cout << "E_21 = " << eigen[2] - eigen[1] << std::endl;

    double s_10 = sin(((eigen[1] - eigen[0]) * L) / (4.0 * E));
    double s_20 = sin(((eigen[2] - eigen[0]) * L) / (4.0 * E));
    double s_21 = sin(((eigen[2] - eigen[1]) * L) / (4.0 * E));


    // Compute probability
    double P = 1.0 - 4.0 * (X[1]*X[0]*s_10*s_10 + X[2]*X[0]*s_20*s_20 + X[2]*X[1]*s_21*s_21);

    return P;
}


/**
 * @brief Survival Probability in a Vacuum (standard calc).
 * 
 * @param filename for PMNS matrix elements and mass differences.
 * @param E Neutrino energy (MeV).
 * @param L Propagation length (km).
 * @return double 
 */
double Survival_Prob_Vac(std::string filename, double E, double L) {
    // convert all units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)

    // Read in elements of the top line of the PMNS matrix and mass differences
    double* dataPoint = read_data(filename, 21);
    double U_Re;
    double U_Im;
    double U_mag2[3];
    for(int i=0; i<3; ++i){
        U_Re = *(dataPoint + i);
        U_Im = *(dataPoint + 9 + i);
        U_mag2[i] = U_Re*U_Re + U_Im*U_Im;
    }

    double m21 = *(dataPoint + 18);
    double m31 = *(dataPoint + 19);
    double m32 = *(dataPoint + 20);

    // Compute survival probability
    double sin21 = sin((m21 * L) / (4.0 * E));
    double sin31 = sin((m31 * L) / (4.0 * E));
    double sin32 = sin((m32 * L) / (4.0 * E));

    double P = 1.0 - 4.0 * (U_mag2[1] * U_mag2[0] * sin21*sin21
                                + U_mag2[2] * U_mag2[0] * sin31*sin31
                                + U_mag2[2] * U_mag2[1] * sin32*sin32);

    return P;
}


/**
 * @brief Survival probability compted by Globes package
 * 
 * @param Ne Electron density (is converted to matter density before being passed to globes function)
 * @param E Neutrino energy
 * @param L Baseline
 * @param anti Is it an antineutrino?
 * @return double 
 */
double Survival_Prob_Globes(double Ne, double E, double L, bool anti) {
    // convert units to eV
    E /= 1e3; //(MeV to GeV)

    // Input parameters
    double delta13 = 1.36 * M_PI;

    double s12 = sqrt(0.297);
    double s13 = sqrt(0.0215);
    double s23 = sqrt(0.545);

    double theta12 = asin(s12);
    double theta13 = asin(s13);
    double theta23 = asin(s23);

    // Define mass differences
    double m21 = 7.53e-5; //(Delta m_21^2, in eV^2)
    double m31 = 0.0025283; //(Delta m_31^2, in eV^2)
    double m32 = m31 - m21; //(Delta m_32^2, in eV^2)

    // Initialise Globes stuff
    // char* exp_file = (char *)"/home/jamicus/Studies/PhD/Antinu/Matter_Effects/Globes/globes-3.0.11/examples/NFstandard.glb";
    /* Initialize experiment NFstandard.glb */
    // glbInitExperiment(exp_file, &glb_experiment_list[0], &glb_num_of_exps);
    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    /* Assign: theta12,theta13,theta23,deltacp,dm2solar,dm2atm */
    glbDefineParams(true_values, theta12, theta13, theta23, delta13, m21, m31);
    //glbSetDensityParams(true_values,1.0,GLB_ALL); // Matter scaling
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // Run Globes function
    int anitInt = 1;
    if(anti){anitInt = -1;}
    Ne *= 1e-6; //(m^-3 to cm^3)
    double rho = Ne * 4.821e-15;  //matter density g/cm^3
    //std::cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << std::endl;
    double P_globes = glbConstantDensityProbability(1, 1, anitInt, E, L, rho);

    return P_globes;
}

/**
 * @brief Survival probability in vacuum compted by Globes package
 * 
 * @param E Neutrino energy
 * @param L Baseline
 * @param anti Is it an antineutrino?
 * @return double 
 */
double Survival_Prob_Vac_Globes(double E, double L, bool anti) {
    // convert units to eV
    E /= 1e3; //(MeV to GeV)

    // Input parameters
    double delta13 = 1.36 * M_PI;

    double s12 = sqrt(0.297);
    double s13 = sqrt(0.0215);
    double s23 = sqrt(0.545);

    double theta12 = asin(s12);
    double theta13 = asin(s13);
    double theta23 = asin(s23);

    // Define mass differences
    double m21 = 7.53e-5; //(Delta m_21^2, in eV^2)
    double m31 = 0.0025283; //(Delta m_31^2, in eV^2)
    double m32 = m31 - m21; //(Delta m_32^2, in eV^2)

    // Initialise Globes stuff
    // char* exp_file = (char *)"/home/jamicus/Studies/PhD/Antinu/Matter_Effects/Globes/globes-3.0.11/examples/NFstandard.glb";
    /* Initialize experiment NFstandard.glb */
    // glbInitExperiment(exp_file, &glb_experiment_list[0], &glb_num_of_exps);
    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    /* Assign: theta12,theta13,theta23,deltacp,dm2solar,dm2atm */
    glbDefineParams(true_values, theta12, theta13, theta23, delta13, m21, m31);
    //glbSetDensityParams(true_values,1.0,GLB_ALL); // Matter scaling
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // Run Globes function
    int anitInt = 1;
    if(anti){anitInt = -1;}
    //std::cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << std::endl;
    double P_globes = glbVacuumProbability(1, 1, anitInt, E, L);

    return P_globes;
}


/* -------------------- TRANSITION PROBS ----------------------- */


/**
 * @brief Compute transition probability of (anti)muon neutrinos to (anti)electron neutrinos,
 * with matter effects.
 * 
 * @param constants Precomputed, independent of neutrino and matter effects, packaged as:
 *   constants[0] << H_ee_vac,
 *   constants[1] << Y_vac,
 *   constants[2] << a0_vac,
 *   constants[3] << a1_vac,
 *   constants[4] << R_H_em,
 *   constants[5] << I_H_em,
 *   constants[6] << R_Y_em,
 *   constants[7] << I_Y_em,
 * @param Ne Matter electron density (m^-3).
 * @param E Neutrino energy (MeV).
 * @param L Propagation length (km).
 * @param anti true=antineutrino, false=neutrino.
 * @return double 
 */
double Transition_Prob(double constants[8], double Ne, double E, double L, bool anti) {
    // convert units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)

    // Initialise constants with vacuum values
    double H_ee = constants[0];
    double Y_ee = constants[1];
    double a0 = constants[2];
    double a1 = constants[3];
    double R_H_em = constants[4];
    double I_H_em = constants[5];
    double R_Y_em = constants[6];
    double I_Y_em = constants[7];

    // std::cout << "H_ee = " << H_ee << std::endl;
    // std::cout << "Y_ee = " << Y_ee << std::endl;
    // std::cout << "a0 = " << a0 << std::endl;
    // std::cout << "a1 = " << a1 << std::endl;
    // std::cout << "R_H_em = " << R_H_em << std::endl;
    // std::cout << "I_H_em = " << I_H_em << std::endl;
    // std::cout << "R_Y_em = " << R_Y_em << std::endl;
    // std::cout << "I_Y_em = " << I_Y_em << std::endl;

    // If not vacuum, make corrections
    double A_CC = 0.0;
    if(Ne != 0.0){
        // convert units to eV
        //std::cout << "Ne = " << Ne << std::endl;
        Ne *= c*c*c * hbar*hbar*hbar; //(m^-3 to ev^3)
        //std::cout << "Ne = " << Ne << std::endl;

        // Work out ACC in (eV)^2
        A_CC = 2 * sqrt(2) * E * GF * Ne;
        if(anti){
            A_CC *= -1;
        }
        std::cout << "A_CC = " << A_CC << std::endl;

        // Compute new values for a0 and a1
        a0 += -Y_ee * A_CC - (1.0/3.0) * H_ee * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
        a1 += -H_ee * A_CC - (1.0/3.0) * A_CC*A_CC;
    }


    // Get eigenvalues of H, and constants X and theta
    double eigen[3];
    double R_X[3];
    double I_X[3];

    double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
    double preFact = 2.0 * sqrt(- a1 / 3.0);

    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - (2.0 * M_PI * i) / 3.0);
        // std::cout << "E_" << i << " = " << eigen[i] << std::endl;

        R_X[i] = ((eigen[i] + (A_CC / 3.0)) * R_H_em + R_Y_em) / (3.0 * eigen[i]*eigen[i] + a1);
        I_X[i] = ((eigen[i] + (A_CC / 3.0)) * I_H_em + I_Y_em) / (3.0 * eigen[i]*eigen[i] + a1);
        // std::cout << "X_" << i << " = " << R_X[i] << "+ i" << I_X[i] << std::endl;
    }

    // std::cout << "E_10 = " << eigen[1] - eigen[0] << std::endl;
    // std::cout << "E_20 = " << eigen[2] - eigen[0] << std::endl;
    // std::cout << "E_21 = " << eigen[2] - eigen[1] << std::endl;

    double Theta_10 = ((eigen[1] - eigen[0]) * L) / (2.0 * E);
    double Theta_20 = ((eigen[2] - eigen[0]) * L) / (2.0 * E);
    double Theta_21 = ((eigen[2] - eigen[1]) * L) / (2.0 * E);

    // std::cout << "theta_10 = " << Theta_10 << std::endl;
    // std::cout << "theta_20 = " << Theta_20 << std::endl;
    // std::cout << "theta_21 = " << Theta_21 << std::endl;

    double s2_10 = sin(Theta_10 / 2.0);
    double s2_20 = sin(Theta_20 / 2.0);
    double s2_21 = sin(Theta_21 / 2.0);
    double s_10 = sin(Theta_10);
    double s_20 = sin(Theta_20);
    double s_21 = sin(Theta_21);

    // temp test
    // R_X[2] = -0.272838; I_X[2] = -0.0681112;
    // R_X[1] = 0.31843; I_X[1] = -0.0287753;
    // R_X[0] = -0.0455914; I_X[0] = 0.0968868;

    double anitInt = 1.0;
    if(anti){anitInt = -1.0;}

    // Compute probability
    double P  =  - 4.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * s2_10*s2_10
                        + (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * s2_20*s2_20
                        + (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * s2_21*s2_21)
        + anitInt * 2.0 * ((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * s_10
                        + (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * s_20
                        + (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * s_21);

    return P;
}


/**
 * @brief Transition Probability of (anti)muon neutrinos to (anti)electron neutrinos,
 * in a Vacuum (standard calc).
 * 
 * @param filename for PMNS matrix elements and mass differences.
 * @param E Neutrino energy (MeV).
 * @param L Propagation length (km).
 * @return double 
 */
double Transition_Prob_Vac(std::string filename, double E, double L, bool anti) {
    // convert all units to eV
    E *= 1e6; //(MeV to eV)
    L *= 1e3 / (c * hbar); //(km to eV^-1)

    // Read in elements of the top line of the PMNS matrix and mass differences
    double* dataPoint = read_data(filename, 21);
    double U_Re[2][3];
    double U_Im[2][3];
    for (int i=0; i<2; ++i) {
        for (int j=0; j<3; ++j) {
            U_Re[i][j] = dataPoint[3*i + j];
            U_Im[i][j] = dataPoint[9 + 3*i + j];

            // std::cout << "U_" << i << j << " = " << U_Re[i][j] << "+ i" << U_Im[i][j] << std::endl;
        }
    }
    double U_em_Re[3];
    double U_em_Im[3];
    for (int k=0; k<3; ++k) {
        U_em_Re[k] =  U_Re[0][k] * U_Re[1][k] + U_Im[0][k] * U_Im[1][k];
        U_em_Im[k] =  U_Im[0][k] * U_Re[1][k] - U_Re[0][k] * U_Im[1][k];
        // std::cout << "U_em_" << k+1 << " = " << U_em_Re[k] << "+ i" << U_em_Im[k] << std::endl;
    }

    double m21 = *(dataPoint + 18);
    double m31 = *(dataPoint + 19);
    double m32 = *(dataPoint + 20);

    // std::cout << "m21 = " << m21 << std::endl;
    // std::cout << "m31 = " << m31 << std::endl;
    // std::cout << "m32 = " << m32 << std::endl;

    // Compute survival probability
    double sin21_2 = sin((m21 * L) / (4.0 * E));
    double sin31_2 = sin((m31 * L) / (4.0 * E));
    double sin32_2 = sin((m32 * L) / (4.0 * E));
    double sin21 = sin((m21 * L) / (2.0 * E));
    double sin31 = sin((m31 * L) / (2.0 * E));
    double sin32 = sin((m32 * L) / (2.0 * E));

    double anitInt = 1.0;
    if(anti){anitInt = -1.0;}

    // Compute probability
    double P  =  - 4.0 * ((U_em_Re[1]*U_em_Re[0] + U_em_Im[1]*U_em_Im[0]) * sin21_2*sin21_2
                        + (U_em_Re[2]*U_em_Re[0] + U_em_Im[2]*U_em_Im[0]) * sin31_2*sin31_2
                        + (U_em_Re[2]*U_em_Re[1] + U_em_Im[2]*U_em_Im[1]) * sin32_2*sin32_2)
        + anitInt * 2.0 * ((U_em_Im[1]*U_em_Re[0]- U_em_Re[1]*U_em_Im[0]) * sin21
                        + (U_em_Im[2]*U_em_Re[0] - U_em_Re[2]*U_em_Im[0]) * sin31
                        + (U_em_Im[2]*U_em_Re[1] - U_em_Re[2]*U_em_Im[1]) * sin32);

    return P;
}


/**
 * @brief Transition probability of (anti)muon neutrinos to (anti)electron neutrinos,
 * compted by Globes package
 * 
 * @param Ne Electron density (is converted to matter density before being passed to globes function)
 * @param E Neutrino energy
 * @param L Baseline
 * @param anti Is it an antineutrino?
 * @return double 
 */
double Transition_Prob_Globes(double Ne, double E, double L, bool anti) {
    // convert units to eV
    E /= 1e3; //(MeV to GeV)

    // Input parameters
    double delta13 = 1.36 * M_PI;

    double s12 = sqrt(0.297);
    double s13 = sqrt(0.0215);
    double s23 = sqrt(0.545);

    double theta12 = asin(s12);
    double theta13 = asin(s13);
    double theta23 = asin(s23);

    // Define mass differences
    double m21 = 7.5288e-05; //(Delta m_21^2, in eV^2)
    double m31 = 0.0025283; //(Delta m_31^2, in eV^2)
    double m32 = 0.00245301; //(Delta m_32^2, in eV^2)

    // Initialise Globes stuff
    // char* exp_file = (char *)"/home/jamicus/Studies/PhD/Antinu/Matter_Effects/Globes/globes-3.0.11/examples/NFstandard.glb";
    /* Initialize experiment NFstandard.glb */
    // glbInitExperiment(exp_file, &glb_experiment_list[0], &glb_num_of_exps);
    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    /* Assign: theta12,theta13,theta23,deltacp,dm2solar,dm2atm */
    glbDefineParams(true_values, theta12, theta13, theta23, delta13, m21, m31);
    //glbSetDensityParams(true_values,1.0,GLB_ALL); // Matter scaling
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // Run Globes function
    int anitInt = 1;
    if(anti){anitInt = -1;}
    Ne *= 1e-6; //(m^-3 to cm^3)
    double rho = Ne * 4.821e-15;  //matter density g/cm^3
    //std::cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << std::endl;
    double P_globes = glbConstantDensityProbability(2, 1, anitInt, E, L, rho);

    return P_globes;
}

/**
 * @brief Transition probability of (anti)muon neutrinos to (anti)electron neutrinos in vacuum,
 * compted by Globes package
 * 
 * @param E Neutrino energy
 * @param L Baseline
 * @param anti Is it an antineutrino?
 * @return double 
 */
double Transition_Prob_Vac_Globes(double E, double L, bool anti) {
    // convert units to eV
    E /= 1e3; //(MeV to GeV)

    // Input parameters
    double delta13 = 1.36 * M_PI;

    double s12 = sqrt(0.297);
    double s13 = sqrt(0.0215);
    double s23 = sqrt(0.545);

    double theta12 = asin(s12);
    double theta13 = asin(s13);
    double theta23 = asin(s23);

    // Define mass differences
    double m21 = 7.5288e-05; //(Delta m_21^2, in eV^2)
    double m31 = 0.0025283; //(Delta m_31^2, in eV^2)
    double m32 = 0.00245301; //(Delta m_32^2, in eV^2)

    // Initialise Globes stuff
    // char* exp_file = (char *)"/home/jamicus/Studies/PhD/Antinu/Matter_Effects/Globes/globes-3.0.11/examples/NFstandard.glb";
    /* Initialize experiment NFstandard.glb */
    // glbInitExperiment(exp_file, &glb_experiment_list[0], &glb_num_of_exps);
    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    /* Assign: theta12,theta13,theta23,deltacp,dm2solar,dm2atm */
    glbDefineParams(true_values, theta12, theta13, theta23, delta13, m21, m31);
    //glbSetDensityParams(true_values,1.0,GLB_ALL); // Matter scaling
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    // Run Globes function
    int anitInt = 1;
    if(anti){anitInt = -1;}
    //std::cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << std::endl;
    double P_globes = glbVacuumProbability(2, 1, anitInt, E, L);

    return P_globes;
}