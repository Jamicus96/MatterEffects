// Header file for cplx class

// header guard:
#ifndef complex
#define complex

//include
#include <iostream>
#include <assert.h>
#include <cmath>
#include <stdexcept>

class CplxToRealException : public std::runtime_error {
    public:
    /**
     * @brief Error from trying to assign cplx value to long double,
     * float or int.
     * 
     */
    CplxToRealException() : 
        runtime_error( "Assigned cplx with non-zero imaginary part to long double, float or int" ) {}
};


class cplx {
    private:
    long double re;
    long double im;

    public:

    /**
     * @brief Construct a new cplx object.
     * Input real and imaginary parts.
     * 
     * @param r 
     * @param i 
     */
    cplx(long double r, long double i) {re = r; im = i;}

    /**
     * @brief Construct a new cplx object.
     * Copy constructor.
     * 
     * @param z 
     */
    cplx(const cplx &z) {re = z.re; im = z.im;}

    /**
     * @brief Construct a new cplx object.
     * Default constructor: z = 0.
     * 
     */
    cplx() {re = 0.0; im = 0.0;}

    /**
     * @brief Explicit conversion: cplx to long double
     * (also works for float and int).
     * 
     * @return long double& 
     */
    operator long double&() {
        if (im == 0.0) {
            return re;
        } else {
            throw CplxToRealException();
        }
    }

/* ~~~~~~~ member functions ~~~~~~~~ */

    /**
     * @brief Magnitude.
     * 
     * @return long double 
     */
    long double mag() {return sqrt(re * re + im * im);}
    /**
     * @brief Magnitude squared.
     * 
     * @return long double 
     */
    long double mag2() {return re * re + im * im;}
    /**
     * @brief Complex phase.
     * Defaults to zero if magnitude is zero.
     * 
     * @return long double 
     */
    long double phase() {
        if (re == 0.0){
            if (im == 0.0) {
                return 0.0;
            } else {
                return 0.5 * M_PI;
            }
        } else {
            return atan(im / re);
        }
    }

    /**
     * @brief Real part Re(z).
     * 
     * @return long double& 
     */
    long double& real() {return re;}
    /**
     * @brief Imaginary part Im(z).
     * 
     * @return long double& 
     */
    long double& imaginary() {return im;}

    /**
     * @brief Complex conjugate.
     * 
     * @return cplx 
     */
    cplx cplxConj() {return cplx(re, -im);}

/* ~~~~~~~~~ operator overload ~~~~~~~~~ */

    /* ~~~~~~~~~ member functions ~~~~~~~~~ */

    void operator += (const cplx z) {re += z.re; im += z.im;}
    void operator += (long double s) {re += s;}

    void operator -= (const cplx z) {re -= z.re; im -= z.im;}
    void operator -= (long double s) {re -= s;}

    void operator *= (const cplx z) {
        long double temp = re;
        re = re * z.re - im * z.im;
        im = im * z.re + z.im * temp;
    }
    void operator *= (long double s) {re *= s;im *= s;}

    void operator /= (const cplx z) {
        long double denom = z.re * z.re + z.im * z.im;
        long double temp = re;
        re = (re * z.re + im * z.im) / denom;
        im = (im * z.re - z.im * temp) / denom;
    }
    void operator /= (long double s) {re /= s; im /= s;}

    /**
     * @brief Assign long double to cplx (also
     * works for float and int).
     * 
     * @param s 
     */
    void operator = (long double s){
        re = s;
        im = 0.0;
    }

/* ~~~~~~~~~~ friend functions ~~~~~~~~~~ */

    friend std::ostream& operator << (std::ostream& os, const cplx& z) {
        if (z.im == 0.0) {
            os << z.re;
        } else if (z.re == 0.0) {
            if (z.im == 1.0) {
                os << "i";
            } else if (z.im == - 1.0) {
                os << "-i";
            } else {
                os << z.im << "i";
            }
        } else if (z.im > 0.0) {
            if (z.im == 1.0) {
                os << z.re << "+i";
            } else {
                os << z.re << "+" << z.im << "i";
            }
        } else {
            if (z.im == -1.0) {
                os << z.re << "-i";
            } else {
                os << z.re << "-" << - z.im << "i";
            }
        }
        return os;
    }

    friend cplx operator + (const cplx& z, const cplx& w) {
        return cplx(z.re + w.re, z.im + w.im);
    }
    friend cplx operator + (const cplx& z, long double s) {
        return cplx(z.re + s, z.im);
    }
    friend cplx operator + (long double s, const cplx& z) {
        return cplx(s + z.re, z.im);
    }

    friend cplx operator - (const cplx& z, const cplx& w) {
        return cplx(z.re - w.re, z.im - w.im);
    }
    friend cplx operator - (const cplx& z, long double s) {
        return cplx(z.re - s, z.im);
    }
    friend cplx operator - (long double s, const cplx& z) {
        return cplx(s - z.re, - z.im);
    }
    friend cplx operator - (const cplx& z) {
        return cplx(- z.re, - z.im);
    }

    friend cplx operator * (const cplx& z, const cplx& w) {
        return cplx(z.re * w.re - z.im * w.im, z.im * w.re + w.im * z.re);
    }
    friend cplx operator * (const cplx& z, long double s) {
        return cplx(z.re * s, z.im * s);
    }
    friend cplx operator * (long double s, const cplx& z) {
        return cplx(s * z.re, s * z.im);
    }

    friend cplx operator / (const cplx& z, const cplx& w) {
        long double denom = w.re * w.re + w.im * w.im;
        return cplx((z.re * w.re + z.im * w.im) / denom, (z.im * w.re - w.im * z.re) / denom);
    }
    friend cplx operator / (const cplx& z, long double s) {
        return cplx(z.re / s, z.im / s);
    }
    friend cplx operator / (long double s, const cplx& z) {
        long double denom = z.re * z.re + z.im * z.im;
        return cplx(s * z.re / denom, - z.im * s / denom);
    }

    friend bool operator == (const cplx& z, const cplx& w) {
        return (z.re == w.re && z.im == w.im);
    }
    friend bool operator == (const cplx& z, long double s) {
        return (z.re == s && z.im == 0);
    }
    friend bool operator == (long double s, const cplx& z) {
        return (z.re == s && z.im == 0);
    }

    friend bool operator != (const cplx& z, const cplx& w) {
        return (z.re != w.re || z.im != w.im);
    }
    friend bool operator != (const cplx& z, long double s) {
        return (z.re != s || z.im != 0);
    }
    friend bool operator != (long double s, const cplx& z) {
        return (z.re != s || z.im != 0);
    }

/*~~~~~~~~~ end of operator overload ~~~~~~~~*/

    /**
     * @brief Create complex number from polar coordinates.
     * r is the magnitude and theta the phase.
     * 
     * @param r 
     * @param theta 
     * @return cplx 
     */
    static cplx cplxPolar(long double r, long double theta) {
        return cplx(r * std::cos(theta), r * std::sin(theta));
    }

    /**
     * @brief Magnitude.
     * 
     * @param z 
     * @return long double 
     */
    static long double mag(cplx& z) {return std::sqrt(z.re * z.re + z.im * z.im);}
    /**
     * @brief Magnitude.
     * 
     * @param s 
     * @return long double 
     */
    static long double mag(long double& s) {return std::fabs(s);}
    /**
     * @brief Magnitude squared.
     * 
     * @param z 
     * @return long double 
     */
    static long double mag2(cplx& z) {return z.re * z.re + z.im * z.im;}
    /**
     * @brief Magnitude squared.
     * 
     * @param s 
     * @return long double 
     */
    static long double mag2(long double& s) {return s * s;}

    /**
     * @brief Complex number to the power of complex number.
     * 
     * @param z 
     * @param w 
     * @return cplx 
     */
    static cplx pow(cplx& z, cplx& w) {
        long double arg = w.real() * z.phase() + w.imaginary() * std::log(w.mag());
        // z^w = |z|^Re(w) * e^(i * (Re(w) * phase(z) + Im(w) * ln(|w|)))
        return std::pow(z.mag(), w.real()) * cplx(std::cos(arg), std::sin(arg));
    }


};



//end header guard
#endif