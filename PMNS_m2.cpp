/**
 * @file PMNSmatrix
 * @author James Page (jp643@sussex.ac.uk)
 * @brief 
 * @version 0.1
 * @date 2021-06-15
 * 
 * @copyright Copyright (c) 2021
 * 
 */
//include
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;
#include "cplx.hpp"


void PMNS(string filename, bool NO);

int main(){
    // create text file to store info
    string filename = "PMNS_m2_data.txt";
    ofstream datafile (filename);

    // Add PMNS matrix values of interst, for a chosen ordering
    PMNS(filename, true);

    return 0;
}



/**
 * @brief Compute PMNS matrix elements for a given ordering and print to file
 * 
 * @param NO true = Normal Ordering, false = Inverted Ordering
 * @return int 
 */
void PMNS(string filename, bool NO) {
    long double delta13 = 1.36 * M_PI;
    cplx phase1 = cplx::cplxPolar(1, -delta13);
    cplx phase2 = cplx::cplxPolar(1, delta13);

    // Input parameters
    long double s12 = sqrt(0.297);
    long double s13 = sqrt(0.0215);
    long double s23;
    if(NO){
        s23 = sqrt(0.545);
    } else {
        s23 = sqrt(0.547);
    }
    long double c12 = sqrt(1.0 - s12 * s12);
    long double c13 = sqrt(1.0 - s13 * s13);
    long double c23 = sqrt(1.0 - s23 * s23);

    // Define the three component matrices to multiply together
    cplx U23[3][3];
    cplx U13[3][3];
    cplx U12[3][3];
    
    U23[0][0] = 1;
    U23[1][1] = c23;
    U23[2][2] = c23;
    U23[1][2] = s23;
    U23[2][1] = -s23;

    U13[1][1] = 1;
    U13[0][0] = c13;
    U13[2][2] = c13;
    U13[0][2] = s13 * phase1;
    U13[2][0] = -s13 * phase2;

    U12[2][2] = 1;
    U12[0][0] = c12;
    U12[1][1] = c12;
    U12[0][1] = s12;
    U12[1][0] = -s12;

    // Compute each PMNS matrix component
    cplx U[3][3];
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            for(int k=0; k<3; ++k){
                for(int l=0; l<3; ++l){
                    U[i][j] += U23[i][k] * U13[k][l] * U12[l][j];
                }
            }
            cout << "U_" << i << j << " = " << U[i][j] << endl;
        }
    }

    // Print components to file in a list, real and then imaginary parts.
    ofstream datafile;
    datafile.open(filename);
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            datafile << U[i][j].real() << "\n";
        }
    }
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            datafile << U[i][j].imaginary() << "\n";
        }
    }

    // Print masses at the end too
    // Define mass differences
    long double m21 = 7.5288e-5; //(Delta m_21^2, in eV^2)
    long double m32;
    if(NO){
        m32 = 2.45301e-3; //(Delta m_32^2, in eV^2)
    } else {
        m32 = -2.546e-3; //(Delta m_32^2, in eV^2)
    }
    long double m31 = m32 + m21; //(Delta m_31^2, in eV^2)
    
    datafile << m21 << "\n";
    datafile << m31 << "\n";
    datafile << m32 << "\n";

    // Also save useful constants, for consistency
    long double GF = 1.663788e-14; //(eV^-2)
    long double hbar = 6.58211957e-16; //(eV.s)
    long double c = 299792458; //(m/s)
}
