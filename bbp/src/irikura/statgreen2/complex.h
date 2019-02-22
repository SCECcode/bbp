//=============================================================================
// complex.h
// Complex Class
// Function: =, ==, !=, +=, -=, *=, /=, +, -, *, /,
//           conjg(), abs(), sqrt(), exp(), cos(), sin(), tan(), log(), pow()
//=============================================================================

#if !defined(MSS_CCOMPLEX_1999_9_9_INCLUDED_)
#define MSS_CCOMPLEX_1999_9_9_INCLUDED_

#include <iostream>
#include <math.h>
using namespace std;

class CComplex {
private:
	double m_dRe; // real number part
	double m_dIm; // imaginary number part

public:
	CComplex(double r = 0.0, double i = 0.0) { m_dRe = r; m_dIm = i; } // constructor
	CComplex(const CComplex& x) { m_dRe = x.m_dRe; m_dIm = x.m_dIm; } // constructor

	friend double real(const CComplex& x) { return (x.m_dRe); } // return real number part
	friend double imag(const CComplex& x) { return (x.m_dIm); } // return imaginary number part

	CComplex operator+(void) const { return (*this); } // +
	CComplex operator-(void) const { return (CComplex(-m_dRe, -m_dIm)); } // -

// =
	CComplex& operator=(const CComplex& x) 
	{
		m_dRe = x.m_dRe;
		m_dIm = x.m_dIm;
		return (*this);
	}

	CComplex& operator=(const double x) 
	{
		m_dRe = x;
		m_dIm = 0.0;
		return (*this);
	}

// ==
	bool operator==(const CComplex& x) const 
	{
		return (((m_dRe == x.m_dRe) && (m_dIm == x.m_dIm)) ? true : false);
	}

	bool operator==(const double x) const 
	{
		return (((m_dRe == x) && (m_dIm == 0.0)) ? true : false);
	}

// !=
	bool operator!=(const CComplex& x) const 
	{
		return (!(*this == x));
	}

	bool operator!=(const double x) const 
	{
		return (!(*this == x));
	}

// +=
	CComplex& operator+=(const CComplex& x) 
	{
		m_dRe += x.m_dRe;
		m_dIm += x.m_dIm;
		return (*this);
	}

	CComplex& operator+=(const double x) 
	{
		m_dRe += x;
		return (*this);
	}

// -=
	CComplex& operator-=(const CComplex& x) 
	{
		m_dRe -= x.m_dRe;
		m_dIm -= x.m_dIm;
		return (*this);
	}

	CComplex& operator-=(const double x) 
	{
		m_dRe -= x;
		return (*this);
	}

// *=
	CComplex& operator*=(const CComplex& x) 
	{
		double dRe = m_dRe;
		double dIm = m_dIm;
		m_dRe = dRe * x.m_dRe - dIm * x.m_dIm;
		m_dIm = dIm * x.m_dRe + dRe * x.m_dIm;
		return (*this);
	}

	CComplex& operator*=(const double x) 
	{
		m_dRe *= x;
		m_dIm *= x;
		return (*this);
	}

// /=
	CComplex& operator/=(const CComplex& x) 
	{
		double dRe = m_dRe;
		double dIm = m_dIm;
		double r, den;
		double dReTmp = fabs(x.m_dRe);
		double dImTmp = fabs(x.m_dIm);
		if (dReTmp >= dImTmp) {
			r = x.m_dIm / x.m_dRe;
			den = x.m_dRe + r * x.m_dIm;
			m_dRe = (dRe + r * dIm) / den;
			m_dIm = (dIm - r * dRe) / den;
		}
		else {
			r = x.m_dRe / x.m_dIm;
			den = x.m_dIm + r * x.m_dRe;
			m_dRe = (dRe * r + dIm) / den;
			m_dIm = (dIm * r - dRe) / den;
		}
		return (*this);
	}

	CComplex& operator/=(const double x) 
	{
		m_dRe /= x;
		m_dIm /= x;
		return (*this);
	}

// plus
	friend CComplex operator+(const CComplex& x, const CComplex& y) 
	{
		return (CComplex(x.m_dRe + y.m_dRe, x.m_dIm + y.m_dIm));
	}

	friend CComplex operator+(const double x, const CComplex& y) 
	{
		return (CComplex(x + y.m_dRe, y.m_dIm));
	}

	friend CComplex operator+(const CComplex& x, const double y) 
	{
		return (CComplex(x.m_dRe + y, x.m_dIm));
	}

// minus
	friend CComplex operator-(const CComplex& x, const CComplex& y) 
	{
		return (CComplex(x.m_dRe - y.m_dRe, x.m_dIm - y.m_dIm));
	}

	friend CComplex operator-(const double x, const CComplex& y) 
	{
		return (CComplex(x - y.m_dRe, -y.m_dIm));
	}

	friend CComplex operator-(const CComplex& x, const double y) 
	{
		return (CComplex(x.m_dRe - y, x.m_dIm));
	}

// times
	friend CComplex operator*(const CComplex& x, const CComplex& y) 
	{
		return (CComplex(x.m_dRe * y.m_dRe - x.m_dIm * y.m_dIm, x.m_dIm * y.m_dRe + x.m_dRe * y.m_dIm));
	}

	friend CComplex operator*(const double x, const CComplex& y) 
	{
		return (CComplex(x * y.m_dRe, x * y.m_dIm));
	}

	friend CComplex operator*(const CComplex& x, const double y) 
	{
		return (CComplex(x.m_dRe * y, x.m_dIm * y));
	}

// divide
	friend CComplex operator/(const CComplex& x, const CComplex& y) 
	{
		double r, den;
		if (fabs(y.m_dRe) >= fabs(y.m_dIm)) {
			r = y.m_dIm / y.m_dRe;
			den = y.m_dRe + r * y.m_dIm;
			return (CComplex((x.m_dRe + r * x.m_dIm) / den, (x.m_dIm - r * x.m_dRe) / den));
		}
		else {
			r = y.m_dRe / y.m_dIm;
			den = y.m_dIm + r * y.m_dRe;
			return (CComplex((x.m_dRe * r + x.m_dIm) / den, (x.m_dIm * r - x.m_dRe) / den));
		}
	}

	friend CComplex operator/(const double x, const CComplex& y) 
	{
		double r, den;
		if (fabs(y.m_dRe) >= fabs(y.m_dIm)) {
			r = y.m_dIm / y.m_dRe;
			den = y.m_dRe + r * y.m_dIm;
			return (CComplex(x / den, -r * x / den));
		}
		else {
			r = y.m_dRe / y.m_dIm;
			den = y.m_dIm + r * y.m_dRe;
			return (CComplex(x * r / den, -x / den));
		}
	}

	friend CComplex operator/(const CComplex& x, const double y) 
	{
		return (CComplex(x.m_dRe / y, x.m_dIm / y));
	}

// conjugate
	friend CComplex conjg(const CComplex& x) 
	{
		return (CComplex(x.m_dRe, -x.m_dIm));
	}

// absolute
	friend double abs(const CComplex& z) 
	{
		double x, y, ans, temp;
		x = fabs(z.m_dRe);
		y = fabs(z.m_dIm);
		if (x == 0.0)
			ans = y;
		else if (y == 0.0)
			ans = x;
		else if (x > y) {
			temp = y / x;
			ans = x * sqrt(1.0 + temp * temp);
		}
		else {
			temp = x / y;
			ans = y * sqrt(1.0 + temp * temp);
		}
		return (ans);
	}

// square root
	friend CComplex sqrt(const CComplex& z) 
	{
		double x, y, w, r, tmp;
		if ((z.m_dRe == 0.0) && (z.m_dIm == 0.0)) {
			return (CComplex(0.0, 0.0));
		}
		else {
			x = fabs(z.m_dRe);
			y = fabs(z.m_dIm);
			if (x >= y) {
				r = y / x;
				w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
			}
			else {
				r = x / y;
				w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
			}
			if (z.m_dRe >= 0.0) {
				return (CComplex(w, z.m_dIm / (2.0 * w)));
			}
			else {
				tmp = (z.m_dIm >= 0.0) ? w : -w;
				return (CComplex(z.m_dIm / (2.0 * tmp), tmp));
			}
		}
	}

// exponential
	friend CComplex exp(const CComplex& z) 
	{
		return (CComplex(exp(z.m_dRe) * cos(z.m_dIm), exp(z.m_dRe) * sin(z.m_dIm)));
	}

// cos
	friend CComplex cos(const CComplex& z) 
	{
		return (CComplex(cos(z.m_dRe) * cosh(z.m_dIm), -sin(z.m_dRe) * sinh(z.m_dIm)));
	}

// sin
	friend CComplex sin(const CComplex& z) 
	{
		return (CComplex(sin(z.m_dRe) * cosh(z.m_dIm), cos(z.m_dRe) * sinh(z.m_dIm)));
	}

// tan
	friend CComplex tan(const CComplex& z)
	{
		return (sin(z) / cos(z));
	}

// log (n = 0)
	friend CComplex log(const CComplex& z)
	{
		double x, y, w, r;
		
		if (z.m_dIm == 0.0) {
			return (CComplex(log(z.m_dRe), 0.0));
		}
		else {
			x = fabs(z.m_dRe);
			y = fabs(z.m_dIm);
			if (x >= y) {
				r = y / x;
				w = 0.5 * log(x * x * (1.0 + r * r));
			}
			else {
				r = x / y;
				w = 0.5 * log(y * y * (1.0 + r * r));
			}
			return (CComplex(w, atan2(z.m_dIm, z.m_dRe)));
		}
	}

// pow (n = 0)
	friend CComplex pow(const CComplex& x, const CComplex& y)
	{
		if (x.m_dRe == 0.0 && x.m_dIm == 0.0) {
			if (y.m_dRe == 0.0 && y.m_dIm == 0.0) {
				return (CComplex(1.0, 0.0));
			}
			else {
				return (CComplex(0.0, 0.0));
			}
		}
		else {
			return (exp(y * log(x)));
		}
	}

	friend CComplex pow(const CComplex& x, const double y)
	{
		if (x.m_dRe == 0.0 && x.m_dIm == 0.0) {
			if (y == 0.0) {
				return (CComplex(1.0, 0.0));
			}
			else {
				return (CComplex(0.0, 0.0));
			}
		}
		else {
			return (exp(y * log(x)));
		}
	}

	friend CComplex pow(const double x, const CComplex& y)
	{
		if (x == 0.0) {
			if (y.m_dRe == 0.0 && y.m_dIm == 0.0) {
				return (CComplex(1.0, 0.0));
			}
			else {
				return (CComplex(0.0, 0.0));
			}
		}
		else {
			return (exp(y * log(x)));
		}
	}
};

inline ostream& operator<<(ostream& s,  const CComplex& x)
{
	return (s << '(' << real(x) << ", " << imag(x) << ')');
}

#endif
