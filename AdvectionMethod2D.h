#pragma once

//#ifndef AdvectionMethod2D_H
//#define AdvectionMethod2D_H
#include "CommonDef.h"
#include "Field2D.h"
#include "LevelSet2D.h"


template <class TT>
class AdvectionMethod2D
{
public:
	static double alpha; // = min(delta x, delta y) // This is a parameter for Heaviside function and delta function.

	AdvectionMethod2D();
	~AdvectionMethod2D();

	static TT Heaviside(const double& constant);
	static TT Heaviside2(const double& constant);
	static TT DeltaFt(const double& constant);
	static TT DeltaFt2(const double& constant);

	static TT sign(const TT& constant);
	static TT sign(const TT& constant, const double & epsilon);
	static TT Plus(const TT& constant);
	static TT Minus(const TT& constant);

	static void ENO3rdDerivation(const Field2D<TT>& ipField, Field2D<TT>& enoDxMinus, Field2D<TT>& enoDxPlus, Field2D<TT>& enoDyMinus, Field2D<TT>& enoDyPlus);
	static void ENO3rdDxMinus(const Field2D<TT>& ipField, Field2D<TT>& enoDxMinus);
	static void ENO3rdDxPlus(const Field2D<TT>& ipField, Field2D<TT>& enoDxPlus);
	static void ENO3rdDyMinus(const Field2D<TT>& ipField, Field2D<TT>& enoDyMinus);
	static void ENO3rdDyPlus(const Field2D<TT>& ipField, Field2D<TT>& enoDyPlus);
	static TT ENOD1x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD2x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD3x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD1y(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD2y(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD3y(const Field2D<TT>& ipField, const int& i, const int& j);

	static void WENO5thDerivation(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus);
	static void WENO5th(const TT& v1, const TT& v2, const TT& v3, const TT& v4, const TT& v5, TT& constant);
	static void WENO5thDxMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus);
	static void WENO5thDxPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus);
	static void WENO5thDyMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus);
	static void WENO5thDyPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus);

	static void LSReinitializationTVDRK3(LS& levelSet, const double& dt);
	static TT ReinitialGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& phi);
	//static void LSReinitializationFE(LS& levelSet, const double& dt, const double& cfl);
	//static void LSReinitializationGS(LS& levelSet, const double& dt);

	static void LSPropagatingTVDRK3(LS& levelSet, const double& dt);
	static void LSPropagatingTVDRK3(LS& levelSet, const FD& velocity, const double& dt);
	//static void LSPropagatingTVDRK3(LS& levelSet, const FV& velocity, const double& dt);
	static void LSPropagatingTVDRK3(LS& levelSet, const FD& velocityX, const FD& velocityY, const double& dt);
	static void LSPropagatingTVDRK3PeriodicX(LS& levelSet, const FD& velocityX, const FD& velocityY, const double& dt);

	static void LSPropagatingEuler(LS& levelSet, const FD& velocity, const double& dt);

	static TT PropagatingGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& sign);


	//// Local Level Set Functions
	static void LLSPropagatingTVDRK3(LS& levelSet, const double& dt);
	static void LLSPropagatingTVDRK3(LS& levelSet, const FD& velocity, const double& dt);
	static void LLSPropagatingTVDRK3(LS& levelSet, const FD& velocityX, const FD& velocityY, const double& dt);

	static void LLSWENO5thDerivation(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus);
	static void LLSWENO5thDxMinus(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus);
	static void LLSWENO5thDxPlus(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus);
	static void LLSWENO5thDyMinus(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus);
	static void LLSWENO5thDyPlus(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus);

	static void LLSReinitializationTVDRK3(LS& levelSet, const double& dt);
	static void LLSReinitializationTVDRK3(LS& levelSet, const double& dt, const int& iter);

	static void LLSQuantityExtension(const LS& levelSet, FD & ipField);

private:

};



//#endif // !AdvectionMethod2D


template<class TT>
AdvectionMethod2D<TT>::AdvectionMethod2D()
{
}

template<class TT>
AdvectionMethod2D<TT>::~AdvectionMethod2D()
{
}

template<class TT>
inline TT AdvectionMethod2D<TT>::Heaviside(const double & constant)
{
	if (constant > alpha)
	{
		return 1.0;
	}
	else if (constant < -alpha)
	{
		return 0.0;
	}
	else
	{
		return (1 + constant / alpha + sin(PI*constant / alpha) / PI) / 2.0;
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::Heaviside2(const double & constant)
{
	if (constant < -alpha)
	{
		return 0.0;
	}
	else if (constant < -0.5*alpha)
	{
		return -1.0 / (6.0*PI)*(1.0 + constant / alpha + sin(PI*constant / alpha));
	}
	else if (constant <= 0.5*alpha)
	{
		return -1.0 / (6.0*PI)*(1.0 + constant / alpha + sin(PI*constant / alpha)) + 1.0 / 3.0 *(2.0 + 4.0*constant / alpha + 2.0 / PI*sin(2.0*PI*constant / alpha));
	}
	else if (constant <= alpha)
	{
		return -1.0 / (6.0)*(1.0 + constant / alpha + sin(PI*constant / alpha)) + 4.0 / 3.0;
	}
	else if (constant > alpha)
	{
		return 1.0;
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::DeltaFt(const double & constant)
{
	if (abs(constant) > alpha)
	{
		return 0.0;
	}
	else
	{
		return (1 + cos(PI*constant / alpha)) / (2 * alpha);
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::DeltaFt2(const double & constant)
{
	if (abs(constant) > alpha)
	{
		return 0.0;
	}
	else if (abs(constant) < 0.5*alpha)
	{
		return -1.0 / (6.0*PI)*(1 + cos(PI*constant / alpha)) + 4.0 / (3.0*PI)*(1 + cos(2.0*PI*constant / alpha));
	}
	else
	{
		return -1.0 / (6.0*PI)*(1 + cos(PI*constant / alpha));
	}
}


template<class TT>
inline TT AdvectionMethod2D<TT>::sign(const TT & constant)
{
	return TT(constant / sqrt(constant*constant + DBL_EPSILON));
}

template<class TT>
inline TT AdvectionMethod2D<TT>::sign(const TT & constant, const double & epsilon)
{
	return TT(constant / sqrt(constant*constant + epsilon*epsilon));
}

template<class TT>
inline TT AdvectionMethod2D<TT>::Plus(const TT & constant)
{
	if (constant >= 0)
	{
		return TT(constant);
	}
	else
	{
		return TT(0);
	}

}

template<class TT>
inline TT AdvectionMethod2D<TT>::Minus(const TT & constant)
{
	if (constant <= 0)
	{
		return TT(constant);
	}
	else
	{
		return TT(0);
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDerivation(const Field2D<TT>& ipField, Field2D<TT>& enoDxMinus, Field2D<TT>& enoDxPlus, Field2D<TT>& enoDyMinus, Field2D<TT>& enoDyPlus)
{
	ENO3rdDxMinus(ipField, enoDxMinus);
	ENO3rdDxPlus(ipField, enoDxPlus);
	ENO3rdDyMinus(ipField, enoDyMinus);
	ENO3rdDyPlus(ipField, enoDyPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDxMinus(const Field2D<TT>& ipField, Field2D<TT>& enoDxMinus)
{
	double Q1, Q2, Q3, c1, c2;
	int k;
	double D2L, D2R;
	double D3L, D3R;
	double dx = ipField.dx;
#pragma omp parallel for private(Q1, Q2, Q3, c1, c2, k, D2L, D2R, D3L, D3R)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 3 || i > ipField.iEnd - 2)
			{
				enoDxMinus(i, j) = ipField.dxMinusPhi(i, j);
			}
			else
			{
				Q1 = ENOD1x(ipField, i - 1, j);
				D2L = ENOD2x(ipField, i - 1, j);
				D2R = ENOD2x(ipField, i, j);
				if (abs(D2L) <= abs(D2R))
				{
					c1 = D2L;
					k = i - 2;
				}
				else
				{
					c1 = D2R;
					k = i - 1;
				}
				Q2 = -c1*dx;

				D3L = ENOD3x(ipField, k, j);
				D3R = ENOD3x(ipField, k + 1, j);
				if (abs(D3L) <= abs(D3R))
				{
					c2 = D3L;
				}
				else
				{
					c2 = D3R;
				}
				Q3 = (double)c2*(3 * (i - k)*(i - k) - 6 * (i - k) + 2)*dx*dx;

				enoDxMinus(i, j) = Q1 + Q2 + Q3;
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDxPlus(const Field2D<TT>& ipField, Field2D<TT>& enoDxPlus)
{
	double Q1, Q2, Q3, c1, c2;
	int k;
	double D2L, D2R;
	double D3L, D3R;
	double dx = ipField.dx;
#pragma omp parallel for private(Q1, Q2, Q3, c1, c2, k, D2L, D2R, D3L, D3R)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{

			if (i < ipField.iStart + 2 || i > ipField.iEnd - 3)
			{
				enoDxPlus(i, j) = ipField.dxPlusPhi(i, j);
			}
			else
			{
				Q1 = ENOD1x(ipField, i, j);
				D2L = ENOD2x(ipField, i, j);
				D2R = ENOD2x(ipField, i + 1, j);
				if (abs(D2L) <= abs(D2R))
				{
					c1 = D2L;
					k = i - 1;
				}
				else
				{
					c1 = D2R;
					k = i;
				}
				Q2 = -c1*dx;

				D3L = ENOD3x(ipField, k, j);
				D3R = ENOD3x(ipField, k + 1, j);
				if (abs(D3L) <= abs(D3R))
				{
					c2 = D3L;
				}
				else
				{
					c2 = D3R;
				}
				Q3 = (double)c2*(3 * (i - k)*(i - k) - 6 * (i - k) + 2)*dx*dx;

				enoDxPlus(i, j) = Q1 + Q2 + Q3;
			}

		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDyMinus(const Field2D<TT>& ipField, Field2D<TT>& enoDyMinus)
{
	double Q1, Q2, Q3, c1, c2;
	int k;
	double D2L, D2R;
	double D3L, D3R;
	double dy = ipField.dy;
#pragma omp parallel for private(Q1, Q2, Q3, c1, c2, k, D2L, D2R, D3L, D3R)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 3 || j > ipField.jEnd - 2)
			{
				enoDyMinus(i, j) = ipField.dyMinusPhi(i, j);
			}
			else
			{
				Q1 = ENOD1y(ipField, i, j - 1);
				D2L = ENOD2y(ipField, i, j - 1);
				D2R = ENOD2y(ipField, i, j);
				if (abs(D2L) <= abs(D2R))
				{
					c1 = D2L;
					k = j - 2;
				}
				else
				{
					c1 = D2R;
					k = j - 1;
				}
				Q2 = -c1*dy;

				D3L = ENOD3y(ipField, i, k);
				D3R = ENOD3y(ipField, i, k + 1);
				if (abs(D3L) <= abs(D3R))
				{
					c2 = D3L;
				}
				else
				{
					c2 = D3R;
				}
				Q3 = (double)c2*(3 * (j - k)*(j - k) - 6 * (j - k) + 2)*dy*dy;

				enoDyMinus(i, j) = Q1 + Q2 + Q3;
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDyPlus(const Field2D<TT>& ipField, Field2D<TT>& enoDyPlus)
{
	double Q1, Q2, Q3, c1, c2;
	int k;
	double D2L, D2R;
	double D3L, D3R;
	double dy = ipField.dy;
#pragma omp parallel for private(Q1, Q2, Q3, c1, c2, k, D2L, D2R, D3L, D3R)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{

			if (j < ipField.jStart + 2 || j > ipField.jEnd - 3)
			{
				enoDyPlus(i, j) = ipField.dyPlusPhi(i, j);
			}
			else
			{
				Q1 = ENOD1y(ipField, i, j);
				D2L = ENOD2y(ipField, i, j);
				D2R = ENOD2y(ipField, i, j + 1);
				if (abs(D2L) <= abs(D2R))
				{
					c1 = D2L;
					k = j - 1;
				}
				else
				{
					c1 = D2R;
					k = j;
				}
				Q2 = -c1*dy;

				D3L = ENOD3y(ipField, i, k);
				D3R = ENOD3y(ipField, i, k + 1);
				if (abs(D3L) <= abs(D3R))
				{
					c2 = D3L;
				}
				else
				{
					c2 = D3R;
				}
				Q3 = (double)c2*(3 * (j - k)*(j - k) - 6 * (j - k) + 2)*dy*dy;

				enoDyPlus(i, j) = Q1 + Q2 + Q3;
			}
		}
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD1x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ipField(i + 1, j) - ipField(i, j)) / ipField.dx;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD2x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD1x(ipField, i, j) - ENOD1x(ipField, i - 1, j)) / (2 * ipField.dx);
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD3x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD2x(ipField, i + 1, j) - ENOD2x(ipField, i, j)) / (3 * ipField.dx);
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD1y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ipField(i, j + 1) - ipField(i, j)) / ipField.dy;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD2y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD1y(ipField, i, j) - ENOD1y(ipField, i, j - 1)) / (2 * ipField.dy);
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD3y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD2y(ipField, i, j + 1) - ENOD2y(ipField, i, j)) / (3 * ipField.dy);
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDerivation(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus)
{
	WENO5thDxMinus(ipField, wenoXMinus);
	WENO5thDxPlus(ipField, wenoXPlus);
	WENO5thDyMinus(ipField, wenoYMinus);
	WENO5thDyPlus(ipField, wenoYPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5th(const TT & v1, const TT & v2, const TT & v3, const TT & v4, const TT & v5, TT & constant)
{
	TT s1, s2, s3;
	TT a1, a2, a3;
	TT w1, w2, w3;

	s1 = 13.0 / 12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) + 1.0 / 4.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	s2 = 13.0 / 12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) + 1.0 / 4.0*(v2 - v4)*(v2 - v4);
	s3 = 13.0 / 12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) + 1.0 / 4.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);

	a1 = 1.0 / 10.0 * 1.0 / ((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	a2 = 6.0 / 10.0 * 1.0 / ((DBL_EPSILON + s2)*(DBL_EPSILON + s2));
	a3 = 3.0 / 10.0 * 1.0 / ((DBL_EPSILON + s3)*(DBL_EPSILON + s3));

	w1 = a1 / (a1 + a2 + a3);
	w2 = a2 / (a1 + a2 + a3);
	w3 = a3 / (a1 + a2 + a3);

	constant = w1*(1.0 / 3.0*v1 - 7.0 / 6.0*v2 + 11.0 / 6.0*v3) + w2*(-1.0 / 6.0*v2 + 5.0 / 6.0*v3 + 1.0 / 3.0*v4) + w3*(1.0 / 3.0*v3 + 5.0 / 6.0*v4 - 1.0 / 6.0*v5);
}


template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDxMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus)
{
	TT v1, v2, v3, v4, v5;
#pragma omp parallel for private(v1,v2,v3,v4,v5)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 3 || i > ipField.iEnd - 2)
			{
				wenoXMinus(i, j) = ipField.dxMinusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i - 2, j) - ipField(i - 3, j))*ipField.oneOverdx;
				v2 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;
				v3 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
				v4 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
				v5 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;

				WENO5th(v1, v2, v3, v4, v5, wenoXMinus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDxPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus)
{
	TT v1, v2, v3, v4, v5;
#pragma omp parallel for private(v1,v2,v3,v4,v5)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 2 || i > ipField.iEnd - 3)
			{
				wenoXPlus(i, j) = ipField.dxPlusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i + 3, j) - ipField(i + 2, j))*ipField.oneOverdx;
				v2 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;
				v3 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
				v4 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
				v5 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;

				WENO5th(v1, v2, v3, v4, v5, wenoXPlus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDyMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus)
{
	TT v1, v2, v3, v4, v5;
#pragma omp parallel for private(v1,v2,v3,v4,v5)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 3 || j > ipField.jEnd - 2)
			{
				wenoYMinus(i, j) = ipField.dyMinusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i, j - 2) - ipField(i, j - 3))*ipField.oneOverdy;
				v2 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;
				v3 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
				v4 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
				v5 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;

				WENO5th(v1, v2, v3, v4, v5, wenoYMinus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDyPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus)
{
	TT v1, v2, v3, v4, v5;
#pragma omp parallel for private(v1,v2,v3,v4,v5)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 2 || j > ipField.jEnd - 3)
			{
				wenoYPlus(i, j) = ipField.dyPlusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i, j + 3) - ipField(i, j + 2))*ipField.oneOverdy;
				v2 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;
				v3 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
				v4 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
				v5 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;

				WENO5th(v1, v2, v3, v4, v5, wenoYPlus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSReinitializationTVDRK3(LS& levelSet, const double& dt)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);


	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
		}
	}

}

template<class TT>
inline TT AdvectionMethod2D<TT>::ReinitialGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& phi)
{
	if (phi <= 0)
	{
		TT aPlus = max(dxPlus, 0.0);
		TT bMinus = min(dxMinus, 0.0);
		TT cPlus = max(dyPlus, 0.0);
		TT dMinus = min(dyMinus, 0.0);

		return TT(sqrt(max(aPlus*aPlus, bMinus*bMinus) + max(cPlus*cPlus, dMinus*dMinus)) - 1.0);
	}
	else
	{
		TT aMinus = min(dxPlus, 0.0);
		TT bPlus = max(dxMinus, 0.0);
		TT cMinus = min(dyPlus, 0.0);
		TT dPlus = max(dyMinus, 0.0);

		return TT(sqrt(max(aMinus*aMinus, bPlus*bPlus) + max(cMinus*cMinus, dPlus*dPlus)) - 1.0);
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSPropagatingTVDRK3(LS & levelSet, const double & dt)
{
	LS originLevelSet = levelSet;
	LS tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);
	clock_t before = clock();
	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);

		}
	}
	double  result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "1-step time : " << result << "\n";
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSPropagatingTVDRK3(LS & levelSet, const FD& velocity, const double & dt)
{
	LS originLevelSet = levelSet;
	LS tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSPropagatingTVDRK3(LS & levelSet, const FD& velocityX, const FD& velocityY, const double& dt)
{
	LS originLevelSet = levelSet;
	LS tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);



	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);

	double tempDxPhi, tempDyPhi;
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}

			k1(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}
	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}
			k2(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}
	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}
			k3(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSPropagatingTVDRK3PeriodicX(LS & levelSet, const FD& velocityX, const FD& velocityY, const double & dt)
{
	LS originLevelSet = levelSet;
	LS tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);



	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);

	double tempDxPhi, tempDyPhi;
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd - 1; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}

			k1(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}
#pragma omp parallel for
	for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
	{
		levelSet(levelSet.grid.iEnd, j) = levelSet(levelSet.grid.iStart, j);
	}

	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd - 1; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}
			k2(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}
#pragma omp parallel for
	for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
	{
		levelSet(levelSet.grid.iEnd, j) = levelSet(levelSet.grid.iStart, j);
	}

	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);
#pragma omp parallel for private(tempDxPhi, tempDyPhi)
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd - 1; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			if (velocityX(i, j) >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (velocityY(i, j) >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}
			k3(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
#pragma omp parallel for
	for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
	{
		levelSet(levelSet.grid.iEnd, j) = levelSet(levelSet.grid.iStart, j);
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSPropagatingEuler(LS & levelSet, const FD& velocity, const double & dt)
{
	LS tempLevelSet(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	WENO5thDxMinus(levelSet.phi, wenoXMinus);
	WENO5thDxPlus(levelSet.phi, wenoXPlus);
	WENO5thDyMinus(levelSet.phi, wenoYMinus);
	WENO5thDyPlus(levelSet.phi, wenoYPlus);

#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			levelSet(i, j) = levelSet(i, j) - dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		}
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::PropagatingGodunov(const TT & dxPlus, const TT & dxMinus, const TT & dyPlus, const TT & dyMinus, const TT & sign)
{
	if (sign <= 0)
	{
		TT aPlus = max(dxPlus, 0.0);
		TT bMinus = min(dxMinus, 0.0);
		TT cPlus = max(dyPlus, 0.0);
		TT dMinus = min(dyMinus, 0.0);

		return TT(sqrt(max(aPlus*aPlus, bMinus*bMinus) + max(cPlus*cPlus, dMinus*dMinus)));
	}
	else
	{
		TT aMinus = min(dxPlus, 0.0);
		TT bPlus = max(dxMinus, 0.0);
		TT cMinus = min(dyPlus, 0.0);
		TT dPlus = max(dyMinus, 0.0);

		return TT(sqrt(max(aMinus*aMinus, bPlus*bPlus) + max(cMinus*cMinus, dPlus*dPlus)));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3(LS & levelSet, const double & dt)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);
	int i, j;
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k1(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);

	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k2(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k3(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));
	}
	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy));
	LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);
	levelSet.UpdateInterface();
	levelSet.UpdateLLS();
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3(LS & levelSet, const FD& velocity, const double & dt)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);
	int i, j;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k1(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k2(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));

	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k < levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k3(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));

	}

	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy));
	LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);
	levelSet.UpdateInterface();
	levelSet.UpdateLLS();
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3(LS & levelSet, const FD& velocityX, const FD& velocityY, const double & dt)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);

	double tempDxPhi, tempDyPhi;
	int i, j;
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;

		if (velocityX(i, j) >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velocityY(i, j) >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k1(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;

		if (velocityX(i, j) >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velocityY(i, j) >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}
		k2(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		if (velocityX(i, j) >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velocityY(i, j) >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}
		k3(i, j) = -velocityX(i, j)*dt*tempDxPhi - velocityY(i, j)*dt*tempDyPhi;
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));
	}

	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy));
	LLSReinitializationTVDRK3(levelSet, dt,reinitialIter);
	levelSet.UpdateInterface();	
	levelSet.UpdateLLS();
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDerivation(const LS& levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus)
{
	LLSWENO5thDxMinus(levelSet, ipField, wenoXMinus);
	LLSWENO5thDxPlus(levelSet, ipField, wenoXPlus);
	LLSWENO5thDyMinus(levelSet, ipField, wenoYMinus);
	LLSWENO5thDyPlus(levelSet, ipField, wenoYPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDxMinus(const LS & levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		if (i < ipField.iStart + 3 || i > ipField.iEnd - 2)
		{
			wenoXMinus(i, j) = ipField.dxMinusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i - 2, j) - ipField(i - 3, j))*ipField.oneOverdx;
			v2 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;
			v3 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
			v4 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
			v5 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;

			WENO5th(v1, v2, v3, v4, v5, wenoXMinus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDxPlus(const LS & levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;

		if (i < ipField.iStart + 2 || i > ipField.iEnd - 3)
		{
			wenoXPlus(i, j) = ipField.dxPlusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i + 3, j) - ipField(i + 2, j))*ipField.oneOverdx;
			v2 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;
			v3 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
			v4 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
			v5 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;

			WENO5th(v1, v2, v3, v4, v5, wenoXPlus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDyMinus(const LS & levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;

		if (j < ipField.jStart + 3 || j > ipField.jEnd - 2)
		{
			wenoYMinus(i, j) = ipField.dyMinusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i, j - 2) - ipField(i, j - 3))*ipField.oneOverdy;
			v2 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;
			v3 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
			v4 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
			v5 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;

			WENO5th(v1, v2, v3, v4, v5, wenoYMinus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDyPlus(const LS & levelSet, const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;

		if (j < ipField.jStart + 2 || j > ipField.jEnd - 3)
		{
			wenoYPlus(i, j) = ipField.dyPlusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i, j + 3) - ipField(i, j + 2))*ipField.oneOverdy;
			v2 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;
			v3 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
			v4 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
			v5 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;

			WENO5th(v1, v2, v3, v4, v5, wenoYPlus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3(LS & levelSet, const double & dt)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	int i, j;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = levelSet.TubeIndex.iStart; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k1(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k2(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.TubeIndex(k).i;
		j = levelSet.TubeIndex(k).j;
		k3(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3(LS & levelSet, const double & dt, const int & iter)
{
	LS originLevelSet = levelSet;

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	int i, j;
	for (int l = 1; l <= iter; l++)
	{
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = levelSet.TubeIndex.iStart; k <= levelSet.numTube; k++)
		{
			i = levelSet.TubeIndex(k).i;
			j = levelSet.TubeIndex(k).j;
			k1(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.TubeIndex(k).i;
			j = levelSet.TubeIndex(k).j;
			k2(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
		}
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.TubeIndex(k).i;
			j = levelSet.TubeIndex(k).j;
			k3(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSQuantityExtension(const LS & levelSet, FD & ipField)
{
}
