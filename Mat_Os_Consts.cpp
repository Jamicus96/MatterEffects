/**
 * @file Mat_Os_Consts.cpp
 * @author James Page (jp643@sussex.ac.uk)
 * @brief Return constants needed for computation.
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
//#include "cplx.hpp"


void Survival_Prob_Constants(std::string filename, double data[20]);
double* read_data(std::string filename, int num);
double Real_Yem(double m21, double m31, double s12, double s13, double s23, double delta);
double Im_Yem(double m21, double m31, double s12, double s13, double s23, double delta);


int main() {
    double Nwater = 3.3e29;

    // Read in data from text file
    double* datapoint = read_data("PMNS_m2_data.txt", 20);
    double data[20];
    for(int i=0; i<20; ++i){
        data[i] = *(datapoint + i);
    }

    // Call function to compute and print needed constants to new file
    Survival_Prob_Constants("Constants.txt", data);

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


/* ------------------------------------------------------- */


void Survival_Prob_Constants(std::string filename, double data[20]) {
    // Unpack data
    double U_Re[3][3];
    double U_Im[3][3];
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            U_Re[i][j] = data[3*i + j];
            U_Im[i][j] = data[9 + 3*i + j];

            std::cout << "U_" << i << j << " = " << U_Re[i][j] << "+ i" << U_Im[i][j] << std::endl;
        }
    }
    for (int k=0; k<3; ++k) {
        std::cout << "U_e" << k+1 << "U_mu" << k+1 << "^* = " << U_Re[0][k] * U_Re[1][k] + U_Im[0][k] * U_Im[1][k]
            << " + i" << U_Im[0][k] * U_Re[1][k] - U_Re[0][k] * U_Im[1][k] << std::endl;
    }
    double m21 = data[18];
    double m31 = data[19];

    std::cout << "m21 = " << m21 << std::endl;
    std::cout << "m31 = " << m31 << std::endl;

    // Compute vacuum constants
    double a0 = - (2.0/27.0) * (m21*m21*m21 + m31*m31*m31) + (1.0/9.0) * (m21*m21 * m31 + m21 * m31*m31);
    double a1 = (1.0/3.0) * (m21 * m31 - m21*m21 - m31*m31);

    // Input parameters
    double s12 = sqrt(0.297);
    double s13 = sqrt(0.0215);
    double s23;
    if(true){
        s23 = sqrt(0.545);
    } else {
        s23 = sqrt(0.547);
    }
    double c12 = sqrt(1.0 - s12 * s12);
    double c13 = sqrt(1.0 - s13 * s13);
    double c23 = sqrt(1.0 - s23 * s23);

    double H_ee = m21 * (s12*s12 * c13*c13 - (1.0/3.0)) + m31 * (s13*s13 - (1.0/3.0));

    double H_neq2 = c13*c13 * (m21*m21 * s12*s12 * (c12*c12 + s12*s12 * s13*s13) + m31*m31 * s13*s13
                    - 2.0 * m21 * m31 * s12*s12 * s13*s13);

    std::cout << H_neq2 << std::endl;

    double Y_ee = (2.0/3.0) * a1 + H_ee*H_ee + H_neq2;

    // Create new file, and print the needed constants to the file
    std::ofstream datafile(filename);
    datafile << H_ee << "\n";
    datafile << Y_ee << "\n";
    datafile << a0 << "\n";
    datafile << a1 << "\n";

    /*~~~~~~~~ mu -> e Extra stuff ~~~~~~~~~*/
    double delta = 1.36 * M_PI;
    double c_delta = cos(delta);
    double s_delta = sin(delta);

    double R_H_em = m21 * s12 * c13 * (c12 * c23 - s12 * s23 * s13 * c_delta)
                    + m31 * s13 * s23 * c13 * c_delta;
    double I_H_em = m21 * s12*s12 * s13 * s23 * c13 * s_delta - m31 * s13 * s23 * c13 * s_delta;

    double R_Y_em = Real_Yem(m21, m31, s12, s13, s23, delta);
    double I_Y_em = Im_Yem(m21, m31, s12, s13, s23, delta);

    datafile << R_H_em << "\n";
    datafile << I_H_em << "\n";
    datafile << R_Y_em << "\n";
    datafile << I_Y_em << "\n";
}

double Real_Yem(double m21, double m31, double s12, double s13, double s23, double delta) {

    double c12 = sqrt(1.0 - s12 * s12);
    double c13 = sqrt(1.0 - s13 * s13);
    double c23 = sqrt(1.0 - s23 * s23);
    double c_delta = cos(delta);

    double R_X2 = (1.0/3.0) * s12 * c13 * (c12 * c23 - s12 * s13 * s23 * c_delta);

    double R_X3 = (1.0 / 3.0) * s13 * s23 * c13 * c_delta;

    double R_X23 = - (2.0/3.0) * c12 * c13 * (s12 * c23 + s13 * s23 * c12 * c_delta);

    return m21*m21 * R_X2 + m31*m31 * R_X3 + m21 * m31 * R_X23;
}

double Im_Yem(double m21, double m31, double s12, double s13, double s23, double delta) {

    double c12 = sqrt(1.0 - s12 * s12);
    double c13 = sqrt(1.0 - s13 * s13);
    double c23 = sqrt(1.0 - s23 * s23);
    double s_delta = sin(delta);

    double I_X2 = (1.0/3.0) * s12*s12 * s13 * s23 * c13 * s_delta;

    double I_X3 = - (1.0 / 3.0) * s13 * s23 * c13 * s_delta;

    double I_X23 = (2.0/3.0) * s13 * s23 * c12*c12 * c13 * s_delta;
    
    return m21*m21 * I_X2 + m31*m31 * I_X3 + m21 * m31 * I_X23;
}