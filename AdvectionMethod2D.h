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
	static double alpha; // = min(delta x, delta y) // This is a parameter for heaviside function and delta function.

	AdvectionMethod2D();
	~AdvectionMethod2D();

	static TT heaviside(const double& constant);
	static TT deltaFt(const double& constant);

	static TT sign(const TT& constant);


	static void WENO5th(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus);
	static void WENO5th(const TT& v1, const TT& v2, const TT& v3, const TT& v4, const TT& v5, TT& constant);
	static void WENO5thApproxXMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus);
	static void WENO5thApproxXPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus);
	static void WENO5thApproxYMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus);
	static void WENO5thApproxYPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus);

	static void levelSetReinitializationTVDRK3(LevelSet2D& levelSet, const double& dt);
	static TT reinitialGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& phi);
	static void levelSetReinitializationFE(LevelSet2D& levelSet, const double& dt, const double& cfl);
	static void levelSetReinitializationGS(LevelSet2D& levelSet, const double& dt);

	static void levelSetPropagatingTVDRK3(LevelSet2D& levelSet, const double& dt);
	static void levelSetPropagatingTVDRK3(LevelSet2D& levelSet, const Field2D<double>& velocity, const double& dt);
	//static void levelSetPropagatingTVDRK3(LevelSet2D& levelSet, const Field2D<Vector2D<double>>& velocity, const double& dt);
	static void levelSetPropagatingTVDRK3(LevelSet2D& levelSet, const Field2D<double>& velocityX, const Field2D<double>& velocityY, const double& dt);

	static void levelSetPropagatingEuler(LevelSet2D& levelSet, const Field2D<double>& velocity, const double& dt);

	static TT propagatingGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& sign);
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
inline TT AdvectionMethod2D<TT>::heaviside(const double & constant)
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
inline TT AdvectionMethod2D<TT>::deltaFt(const double & constant)
{
	if (abs(constant)>alpha)
	{
		return 0.0;
	}
	else
	{
		return (1 + cos(PI*constant / alpha)) / (2 * alpha);
	}
}


template<class TT>
inline TT AdvectionMethod2D<TT>::sign(const TT & constant)
{
	return TT(constant / sqrt(constant*constant + DBL_EPSILON));
}


template<class TT>
inline void AdvectionMethod2D<TT>::WENO5th(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus, Field2D<TT>& wenoXPlus, Field2D<TT>& wenoYMinus, Field2D<TT>& wenoYPlus)
{
	WENO5thApproxXMinus(ipField, wenoXMinus);
	WENO5thApproxXPlus(ipField, wenoXPlus);
	WENO5thApproxYMinus(ipField, wenoYMinus);
	WENO5thApproxYPlus(ipField, wenoYPlus);
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
inline void AdvectionMethod2D<TT>::WENO5thApproxXMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoXMinus)
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
inline void AdvectionMethod2D<TT>::WENO5thApproxXPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoXPlus)
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
inline void AdvectionMethod2D<TT>::WENO5thApproxYMinus(const Field2D<TT>& ipField, Field2D<TT>& wenoYMinus)
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
inline void AdvectionMethod2D<TT>::WENO5thApproxYPlus(const Field2D<TT>& ipField, Field2D<TT>& wenoYPlus)
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
inline void AdvectionMethod2D<TT>::levelSetReinitializationTVDRK3(LevelSet2D& levelSet, const double& dt)
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);


	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -sign(originLevelSet(i, j))*dt*reinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -sign(originLevelSet(i, j))*dt*reinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -sign(originLevelSet(i, j))*dt*reinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}

}

template<class TT>
inline TT AdvectionMethod2D<TT>::reinitialGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& phi)
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
inline void AdvectionMethod2D<TT>::levelSetPropagatingTVDRK3(LevelSet2D & levelSet, const double & dt)
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -dt*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -dt*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -dt*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);

		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::levelSetPropagatingTVDRK3(LevelSet2D & levelSet, const Field2D<double>& velocity, const double & dt)
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -dt*velocity(i, j)*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -dt*velocity(i, j)*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -dt*velocity(i, j)*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::levelSetPropagatingTVDRK3(LevelSet2D & levelSet, const Field2D<double>& velocityX, const Field2D<double>& velocityY, const double& dt)
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);



	WENO5thApproxXMinus(levelSet.phi, wenoXMinus);
	WENO5thApproxXPlus(levelSet.phi, wenoXPlus);
	WENO5thApproxYMinus(levelSet.phi, wenoYMinus);
	WENO5thApproxYPlus(levelSet.phi, wenoYPlus);

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
	WENO5thApproxXMinus(levelSet.phi, wenoXMinus);
	WENO5thApproxXPlus(levelSet.phi, wenoXPlus);
	WENO5thApproxYMinus(levelSet.phi, wenoYMinus);
	WENO5thApproxYPlus(levelSet.phi, wenoYPlus);
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
	WENO5thApproxXMinus(levelSet.phi, wenoXMinus);
	WENO5thApproxXPlus(levelSet.phi, wenoXPlus);
	WENO5thApproxYMinus(levelSet.phi, wenoYMinus);
	WENO5thApproxYPlus(levelSet.phi, wenoYPlus);
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
inline void AdvectionMethod2D<TT>::levelSetPropagatingEuler(LevelSet2D & levelSet, const Field2D<double>& velocity, const double & dt)
{
	LevelSet2D tempLevelSet(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

	WENO5thApproxXMinus(levelSet.phi, wenoXMinus);
	WENO5thApproxXPlus(levelSet.phi, wenoXPlus);
	WENO5thApproxYMinus(levelSet.phi, wenoYMinus);
	WENO5thApproxYPlus(levelSet.phi, wenoYPlus);

#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			levelSet(i, j) = levelSet(i, j) - dt*velocity(i, j)*propagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		}
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::propagatingGodunov(const TT & dxPlus, const TT & dxMinus, const TT & dyPlus, const TT & dyMinus, const TT & sign)
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
