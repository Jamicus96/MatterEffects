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
using namespace std;
#include <globes/globes.h>


double Survival_Prob(double constants[4], double Ne, double E, double L, bool anti=true);
double Survival_Prob_Vac(string filename, double E, double L);
double Survival_Prob_Globes(double Ne, double E, double L, bool anti=true);
double Survival_Prob_Vac_Globes(double E, double L, bool anti=true);
double Transition_Prob(double constants[8], double Ne, double E, double L, bool anti=true);
double Transition_Prob_Vac(string filename, double E, double L, bool anti=true);
double Transition_Prob_Globes(double Ne, double E, double L, bool anti=true);
double Transition_Prob_Vac_Globes(double E, double L, bool anti=true);
double* read_data(string filename, int num);

// Global variables
double GF = 1.663788e-14; //(eV^-2)
double hbar = 6.58211957e-16; //(eV.s)
double c = 299792458; //(m/s)


int main(int argc, char *argv[]) {
    // Read in calculated constants
    double* constantsPoint = read_data("Constants.txt", 8);
    double survival_constants[4];
    double transition_constants[8];
    for(int i=0; i<8; ++i){
        transition_constants[i] = *(constantsPoint + i);
        // cout << "const_" << i << " = " << transition_constants[i] << endl;
        if (i < 4) {
            survival_constants[i] = *(constantsPoint + i);
        }
    }

    // Call function to compute and print needed constants to new file, and compare to vac
    double E = atof(argv[1]); // MeV
    double L = atof(argv[2]); // km
    double Ne = atof(argv[3]); // m^-3
    int type = atoi(argv[4]); // 0=survival, 1=transition

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
    
    // cout << "P = " << P << ", P_vac = " << P_vac << ", P_globes = " << P_globes << endl;

    // Print results to file
    ofstream datafile;
    datafile.open("results.txt", std::ios::app);
    datafile << P << " " << P_vac << " " << P_globes << " " << P_vac_globes << endl;

    return 0;
}



/**
 * @brief Read data from text file, with one number per line, into array of given length
 * 
 * @param filename 
 * @param num 
 * @return double*
 */
double* read_data(string filename, int num) {
    //Create an input file stream
	ifstream input(filename, ios::in);

	double* array = new double[num];
	for(int i=0; i<num; ++i){
        input >> array[i];
    }
    
    return array;
}


/* -------------------- SURVIVAL PROBS ----------------------- */


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

    // cout << "H_ee = " << H << endl;
    // cout << "Y_ee = " << Y << endl;
    // cout << "a0 = " << a0 << endl;
    // cout << "a1 = " << a1 << endl;

    // If not vacuum, make corrections
    if(Ne != 0.0){
        // convert units to eV
        //cout << "Ne = " << Ne << endl;
        Ne *= c*c*c * hbar*hbar*hbar; //(m^-3 to ev^3)
        //cout << "Ne = " << Ne << endl;

        // Work out ACC in (eV)^2
        double A_CC = 2 * sqrt(2) * E * GF * Ne;
        if(anti){
            A_CC *= -1;
        }
        //cout << "A_CC = " << A_CC << endl;

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
        //cout << "E_" << i << " = " << eigen[i] << endl;

        X[i] = (1.0/3.0) + (eigen[i] * H + Y) / (3.0 * eigen[i]*eigen[i] + a1);
        //cout << "X_" << i << " = " << X[i] << endl;
    }

    // cout << "E_10 = " << eigen[1] - eigen[0] << endl;
    // cout << "E_20 = " << eigen[2] - eigen[0] << endl;
    // cout << "E_21 = " << eigen[2] - eigen[1] << endl;

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
double Survival_Prob_Vac(string filename, double E, double L) {
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
    //cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << endl;
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
    //cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << endl;
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

    // cout << "H_ee = " << H_ee << endl;
    // cout << "Y_ee = " << Y_ee << endl;
    // cout << "a0 = " << a0 << endl;
    // cout << "a1 = " << a1 << endl;
    // cout << "R_H_em = " << R_H_em << endl;
    // cout << "I_H_em = " << I_H_em << endl;
    // cout << "R_Y_em = " << R_Y_em << endl;
    // cout << "I_Y_em = " << I_Y_em << endl;

    // If not vacuum, make corrections
    double A_CC = 0.0;
    if(Ne != 0.0){
        // convert units to eV
        //cout << "Ne = " << Ne << endl;
        Ne *= c*c*c * hbar*hbar*hbar; //(m^-3 to ev^3)
        //cout << "Ne = " << Ne << endl;

        // Work out ACC in (eV)^2
        A_CC = 2 * sqrt(2) * E * GF * Ne;
        if(anti){
            A_CC *= -1;
        }
        cout << "A_CC = " << A_CC << endl;

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
        // cout << "E_" << i << " = " << eigen[i] << endl;

        R_X[i] = ((eigen[i] + (A_CC / 3.0)) * R_H_em + R_Y_em) / (3.0 * eigen[i]*eigen[i] + a1);
        I_X[i] = ((eigen[i] + (A_CC / 3.0)) * I_H_em + I_Y_em) / (3.0 * eigen[i]*eigen[i] + a1);
        // cout << "X_" << i << " = " << R_X[i] << "+ i" << I_X[i] << endl;
    }

    // cout << "E_10 = " << eigen[1] - eigen[0] << endl;
    // cout << "E_20 = " << eigen[2] - eigen[0] << endl;
    // cout << "E_21 = " << eigen[2] - eigen[1] << endl;

    double Theta_10 = ((eigen[1] - eigen[0]) * L) / (2.0 * E);
    double Theta_20 = ((eigen[2] - eigen[0]) * L) / (2.0 * E);
    double Theta_21 = ((eigen[2] - eigen[1]) * L) / (2.0 * E);

    // cout << "theta_10 = " << Theta_10 << endl;
    // cout << "theta_20 = " << Theta_20 << endl;
    // cout << "theta_21 = " << Theta_21 << endl;

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
double Transition_Prob_Vac(string filename, double E, double L, bool anti) {
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

            // cout << "U_" << i << j << " = " << U_Re[i][j] << "+ i" << U_Im[i][j] << endl;
        }
    }
    double U_em_Re[3];
    double U_em_Im[3];
    for (int k=0; k<3; ++k) {
        U_em_Re[k] =  U_Re[0][k] * U_Re[1][k] + U_Im[0][k] * U_Im[1][k];
        U_em_Im[k] =  U_Im[0][k] * U_Re[1][k] - U_Re[0][k] * U_Im[1][k];
        // cout << "U_em_" << k+1 << " = " << U_em_Re[k] << "+ i" << U_em_Im[k] << endl;
    }

    double m21 = *(dataPoint + 18);
    double m31 = *(dataPoint + 19);
    double m32 = *(dataPoint + 20);

    // cout << "m21 = " << m21 << endl;
    // cout << "m31 = " << m31 << endl;
    // cout << "m32 = " << m32 << endl;

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
    //cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << endl;
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
    //cout << ((sqrt(2) * GF) / (7.5e-14 * 0.5)) * 1e6 * c*c*c * hbar*hbar*hbar << endl;
    double P_globes = glbVacuumProbability(2, 1, anitInt, E, L);

    return P_globes;
}