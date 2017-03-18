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

	static double MinMod(const double & constant1, const double & constant2);
	static double MinAbs(const double & constant1, const double & constant2);

	static void ENO3rdDerivation(const Field2D<TT>& ipField, Array2D<TT>& enoDxMinus, Array2D<TT>& enoDxPlus, Array2D<TT>& enoDyMinus, Array2D<TT>& enoDyPlus);
	static void ENO3rdDxMinus(const Field2D<TT>& ipField, Array2D<TT>& enoDxMinus);
	static void ENO3rdDxPlus(const Field2D<TT>& ipField, Array2D<TT>& enoDxPlus);
	static void ENO3rdDyMinus(const Field2D<TT>& ipField, Array2D<TT>& enoDyMinus);
	static void ENO3rdDyPlus(const Field2D<TT>& ipField, Array2D<TT>& enoDyPlus);
	static TT ENOD1x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD2x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD3x(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD1y(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD2y(const Field2D<TT>& ipField, const int& i, const int& j);
	static TT ENOD3y(const Field2D<TT>& ipField, const int& i, const int& j);

	static void WENO3rdDerivation(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus);
	static void WENO3rd(const TT& v1, const TT& v2, const TT& v3, TT& constant);
	static void WENO3rdDxMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus);
	static void WENO3rdDxPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus);
	static void WENO3rdDyMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus);
	static void WENO3rdDyPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus);

	static void WENO5thDerivation(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus);
	static void WENO5th(const TT& v1, const TT& v2, const TT& v3, const TT& v4, const TT& v5, TT& constant);
	static void WENO5thDxMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus);
	static void WENO5thDxPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus);
	static void WENO5thDyMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus);
	static void WENO5thDyPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus);

	static void LSReinitializationTVDRK3(LS& levelSet, const double& dt);
	static void LSReinitializationTVDRK3(LS& levelSet, const double& dt, const int& iteration);
	static void LSReinitializationTVDRK3(LS& levelSet, LS& originLS, const double& dt, const int& spatialOrder);
	static void LSOneSidedDerivativesSubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LSOneSidedDerivativesDxSubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus);
	static void LSOneSidedDerivativesDySubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LSReinitializationTVDRK3SubcellFixSecondOrder(LS& levelSet, const double& dt, const int& iter);
	static TT ReinitialGodunov(const TT& dxPlus, const TT& dxMinus, const TT& dyPlus, const TT& dyMinus, const TT& phi);
	static TT ReinitialGodunov(LS& ipLS, const int& i, const int & j);
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
	static void LLSPropagatingTVDRK3MACGrid(LS& levelSet, const FD& velocityX, const FD& velocityY, const double& dt);
	static void LLSPropagatingTVDRK3MACGrid(LS& levelSet, const FD& velocityX, const FD& velocityY, const double& dt, const int& spatialOrder);

	static void LLSWENO3rdDerivation(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus);
	static void LLSWENO3rdDxMinus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus);
	static void LLSWENO3rdDxPlus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus);
	static void LLSWENO3rdDyMinus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus);
	static void LLSWENO3rdDyPlus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus);

	static void LLSWENO5thDerivation(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus);
	static void LLSWENO5thDxMinus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus);
	static void LLSWENO5thDxPlus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus);
	static void LLSWENO5thDyMinus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus);
	static void LLSWENO5thDyPlus(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus);

	static void LLSReinitializationTVDRK3(LS& levelSet, const double& dt);
	static void LLSReinitializationTVDRK3(LS& levelSet, const double& dt, const int& iter);
	static void LLSReinitializationTVDRK3(LS& levelSet, const double& dt, const int& iter, const int& spatialOrder);
	static void LLSReinitializationTVDRK3usingSubcellFix(LS& levelSet, const double& dt, const int& iter, const int& spatialOrder);

	static void LLSOneSidedDerivativesSubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LLSOneSidedDerivativesDxSubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus);
	static void LLSOneSidedDerivativesDySubcellFixSecondOrder(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LLSReinitializationTVDRK3SubcellFixSecondOrder(LS& levelSet, const double& dt, const int& iter);

	static void OneSideLeftDerivatives(const double& lsL, const double& lsC, const double& valueL, const double& valueC, const double& dx, double& val, double thetadx);
	static void OneSideRightDerivatives(const double& lsC, const double& lsR, const double& valueC, const double& valueR, const double& dx, double& val, double thetadx);
	static void LLSOneSidedDerivatives(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LLSOneSidedDerivativesDx(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus);
	static void LLSOneSidedDerivativesDy(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus);
	static void LLSReinitializationTVDRK3SubcellFixHighOrder(LS& levelSet, const double& dt, const int& iter);

	static void LLSQuantityExtension(LS& ipLS, FD& ipQuantity, const int & temporalOrder, const int & spatialOrder);
	static void LLSQuantityExtension(LS& ipLS, FD& ipQuantity, const int & temporalOrder, const int & spatialOrder, const int & iter);

	// Adaptive time step functions.
	static double AdaptiveTimeStep(const FD& velocity1, const double & cflCondition);
	static double AdaptiveTimeStep(const FD& velocity1, const FD& velocity2, const double & cflCondition);

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
		return 0.5 * (1 + constant / alpha + sin(PI*constant / alpha) / PI);
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
inline double AdvectionMethod2D<TT>::MinMod(const double & constant1, const double & constant2)
{
	if (constant1*constant2 <= 0)
	{
		return 0.0;
	}
	else if (abs(constant1) <= abs(constant2))
	{
		return constant1;
	}
	else
	{
		return constant2;
	}
}

template<class TT>
inline double AdvectionMethod2D<TT>::MinAbs(const double & constant1, const double & constant2)
{
	if (abs(constant1) <= abs(constant2))
	{
		return constant1;
	}
	else
	{
		return constant2;
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDerivation(const Field2D<TT>& ipField, Array2D<TT>& enoDxMinus, Array2D<TT>& enoDxPlus, Array2D<TT>& enoDyMinus, Array2D<TT>& enoDyPlus)
{
	ENO3rdDxMinus(ipField, enoDxMinus);
	ENO3rdDxPlus(ipField, enoDxPlus);
	ENO3rdDyMinus(ipField, enoDyMinus);
	ENO3rdDyPlus(ipField, enoDyPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDxMinus(const Field2D<TT>& ipField, Array2D<TT>& enoDxMinus)
{
	double dx(ipField.dx), dy(ipField.dy);
	double one_over_dx(ipField.oneOverdx), one_over_dy(ipField.oneOverdy), one_over_2dx(ipField.oneOver2dx), one_over_2dy(ipField.oneOver2dx);
	double one_over_3dx(one_over_dx*(double)1 / 3), one_over_3dy(one_over_dy*(double)1 / 3);

#pragma omp parallel for
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 3 || i > ipField.iEnd - 2)
			{
				enoDxMinus(i, j) = ipField.dxMinusPhi(i, j);
				continue;
			}
			double diff_1_x_n_3, diff_1_x_n_2, diff_1_x_n_1, diff_1_x_0, diff_1_x_p_1, diff_2_x_n_2, diff_2_x_n_1, diff_2_x_0, diff_2_x_p_1, diff_3_x_n_2, diff_3_x_n_1, diff_3_x_0;

			diff_1_x_n_3 = one_over_dx*(ipField(i - 2, j) - ipField(i - 3, j));
			diff_1_x_n_2 = one_over_dx*(ipField(i - 1, j) - ipField(i - 2, j));
			diff_1_x_n_1 = one_over_dx*(ipField(i, j) - ipField(i - 1, j));
			diff_1_x_0 = one_over_dx*(ipField(i + 1, j) - ipField(i, j));
			diff_1_x_p_1 = one_over_dx*(ipField(i + 2, j) - ipField(i + 1, j));

			diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
			diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
			diff_2_x_0 = one_over_2dx*(diff_1_x_0 - diff_1_x_n_1);
			diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0);

			diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
			diff_3_x_n_1 = one_over_3dx*(diff_2_x_0 - diff_2_x_n_1);
			diff_3_x_0 = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0);

			
			if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
			{
				if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
				{
					enoDxMinus(i, j) = diff_1_x_n_1 + diff_2_x_n_1*dx + 2 * diff_3_x_n_2*dx*dx;
				}
				else
				{
					enoDxMinus(i, j) = diff_1_x_n_1 + diff_2_x_n_1*dx + 2 * diff_3_x_n_1*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					enoDxMinus(i, j) = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					enoDxMinus(i, j) = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDxPlus(const Field2D<TT>& ipField, Array2D<TT>& enoDxPlus)
{
	double dx(ipField.dx), dy(ipField.dy);
	double one_over_dx(ipField.oneOverdx), one_over_dy(ipField.oneOverdy), one_over_2dx(ipField.oneOver2dx), one_over_2dy(ipField.oneOver2dx);
	double one_over_3dx(one_over_dx*(double)1 / 3), one_over_3dy(one_over_dy*(double)1 / 3);

#pragma omp parallel for
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{

			if (i < ipField.iStart + 2 || i > ipField.iEnd - 3)
			{
				enoDxPlus(i, j) = ipField.dxPlusPhi(i, j);
				continue;
			}
			double diff_1_x_n_2, diff_1_x_n_1, diff_1_x_0, diff_1_x_p_1, diff_1_x_p_2, diff_2_x_n_1, diff_2_x_0, diff_2_x_p_1, diff_2_x_p_2, diff_3_x_n_1, diff_3_x_0, diff_3_x_p_1;

			diff_1_x_n_2 = one_over_dx*(ipField(i - 1, j) - ipField(i - 2, j));
			diff_1_x_n_1 = one_over_dx*(ipField(i, j) - ipField(i - 1, j));
			diff_1_x_0 = one_over_dx*(ipField(i + 1, j) - ipField(i, j));
			diff_1_x_p_1 = one_over_dx*(ipField(i + 2, j) - ipField(i + 1, j));
			diff_1_x_p_2 = one_over_dx*(ipField(i + 3, j) - ipField(i + 2, j));

			diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
			diff_2_x_0 = one_over_2dx*(diff_1_x_0 - diff_1_x_n_1);
			diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0);
			diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

			diff_3_x_n_1 = one_over_3dx*(diff_2_x_0 - diff_2_x_n_1);
			diff_3_x_0 = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0);
			diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

			if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					enoDxPlus(i, j) = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					enoDxPlus(i, j) = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
				{
					enoDxPlus(i, j) = diff_1_x_0 - diff_2_x_p_1*dx + 2 * diff_3_x_0*dx*dx;
				}
				else
				{
					enoDxPlus(i, j) = diff_1_x_0 - diff_2_x_p_1*dx + 2 * diff_3_x_p_1*dx*dx;
				}
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDyMinus(const Field2D<TT>& ipField, Array2D<TT>& enoDyMinus)
{
	double dx(ipField.dx), dy(ipField.dy);
	double one_over_dx(ipField.oneOverdx), one_over_dy(ipField.oneOverdy), one_over_2dx(ipField.oneOver2dx), one_over_2dy(ipField.oneOver2dx);
	double one_over_3dx(one_over_dx*(double)1 / 3), one_over_3dy(one_over_dy*(double)1 / 3);

#pragma omp parallel for
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 3 || j > ipField.jEnd - 2)
			{
				enoDyMinus(i, j) = ipField.dyMinusPhi(i, j);
				continue;
			}
			double diff_1_y_n_3, diff_1_y_n_2, diff_1_y_n_1, diff_1_y_0, diff_1_y_p_1, diff_1_y_p_2, diff_2_y_n_2, diff_2_y_n_1, diff_2_y_0, diff_2_y_p_1, diff_2_y_p_2, diff_3_y_n_2, diff_3_y_n_1, diff_3_y_0, diff_3_y_p_1;

			diff_1_y_n_3 = one_over_dy*(ipField(i, j - 2) - ipField(i, j - 3));
			diff_1_y_n_2 = one_over_dy*(ipField(i, j - 1) - ipField(i, j - 2));
			diff_1_y_n_1 = one_over_dy*(ipField(i, j) - ipField(i, j - 1));
			diff_1_y_0 = one_over_dy*(ipField(i, j + 1) - ipField(i, j));
			diff_1_y_p_1 = one_over_dy*(ipField(i, j + 2) - ipField(i, j + 1));

			diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
			diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
			diff_2_y_0 = one_over_2dy*(diff_1_y_0 - diff_1_y_n_1);
			diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0);

			diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
			diff_3_y_n_1 = one_over_3dy*(diff_2_y_0 - diff_2_y_n_1);
			diff_3_y_0 = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0);

			if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
			{
				if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
				{
					enoDyMinus(i, j) = diff_1_y_n_1 + diff_2_y_n_1*dy + 2 * diff_3_y_n_2*dy*dy;
				}
				else
				{
					enoDyMinus(i, j) = diff_1_y_n_1 + diff_2_y_n_1*dy + 2 * diff_3_y_n_1*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					enoDyMinus(i, j) = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					enoDyMinus(i, j) = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}
			
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::ENO3rdDyPlus(const Field2D<TT>& ipField, Array2D<TT>& enoDyPlus)
{
	double dx(ipField.dx), dy(ipField.dy);
	double one_over_dx(ipField.oneOverdx), one_over_dy(ipField.oneOverdy), one_over_2dx(ipField.oneOver2dx), one_over_2dy(ipField.oneOver2dx);
	double one_over_3dx(one_over_dx*(double)1 / 3), one_over_3dy(one_over_dy*(double)1 / 3);

#pragma omp parallel for
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 2 || j > ipField.jEnd - 3)
			{
				enoDyPlus(i, j) = ipField.dyPlusPhi(i, j);
				continue;
			}
			double diff_1_y_n_3, diff_1_y_n_2, diff_1_y_n_1, diff_1_y_0, diff_1_y_p_1, diff_1_y_p_2, diff_2_y_n_2, diff_2_y_n_1, diff_2_y_0, diff_2_y_p_1, diff_2_y_p_2, diff_3_y_n_2, diff_3_y_n_1, diff_3_y_0, diff_3_y_p_1;

			diff_1_y_n_2 = one_over_dy*(ipField(i, j - 1) - ipField(i, j - 2));
			diff_1_y_n_1 = one_over_dy*(ipField(i, j) - ipField(i, j - 1));
			diff_1_y_0 = one_over_dy*(ipField(i, j + 1) - ipField(i, j));
			diff_1_y_p_1 = one_over_dy*(ipField(i, j + 2) - ipField(i, j + 1));
			diff_1_y_p_2 = one_over_dy*(ipField(i, j + 3) - ipField(i, j + 2));

			diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
			diff_2_y_0 = one_over_2dy*(diff_1_y_0 - diff_1_y_n_1);
			diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0);
			diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

			diff_3_y_n_1 = one_over_3dy*(diff_2_y_0 - diff_2_y_n_1);
			diff_3_y_0 = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0);
			diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);
			
			if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					enoDyPlus(i, j) = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					enoDyPlus(i, j) = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
				{
					enoDyPlus(i, j) = diff_1_y_0 - diff_2_y_p_1*dy + 2 * diff_3_y_0*dy*dy;
				}
				else
				{
					enoDyPlus(i, j) = diff_1_y_0 - diff_2_y_p_1*dy + 2 * diff_3_y_p_1*dy*dy;
				}
			}

		}
	}
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD1x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ipField(i + 1, j) - ipField(i, j)) * ipField.oneOverdx;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD2x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD1x(ipField, i, j) - ENOD1x(ipField, i - 1, j)) * ipField.oneOverdx;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD3x(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD2x(ipField, i + 1, j) - ENOD2x(ipField, i, j)) * ipField.oneOverdx;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD1y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ipField(i, j + 1) - ipField(i, j)) * ipField.oneOverdy;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD2y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD1y(ipField, i, j) - ENOD1y(ipField, i, j - 1)) * ipField.oneOverdy;
}

template<class TT>
inline TT AdvectionMethod2D<TT>::ENOD3y(const Field2D<TT>& ipField, const int& i, const int& j)
{
	return (ENOD2y(ipField, i, j + 1) - ENOD2y(ipField, i, j)) * ipField.oneOverdy;
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rdDerivation(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus)
{
	WENO3rdDxMinus(ipField, wenoXMinus);
	WENO3rdDxPlus(ipField, wenoXPlus);
	WENO3rdDyMinus(ipField, wenoYMinus);
	WENO3rdDyPlus(ipField, wenoYPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rd(const TT & v1, const TT & v2, const TT & v3, TT & constant)
{
	TT s1, s2;
	TT a1, a2;
	TT w1, w2;

	s1 = (v1 - v2)*(v1 - v2);
	s2 = (v2 - v3)*(v2 - v3);

	a1 = 1.0 / 3.0 / ((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	a2 = 2.0 / 3.0 / ((DBL_EPSILON + s2)*(DBL_EPSILON + s2));

	w1 = a1 / (a1 + a2);
	w2 = a2 / (a1 + a2);

	constant = w1*(-1.0 / 2.0*v1 + 3.0 / 2.0*v2) + w2*(1.0 / 2.0*v2 + 1.0 / 2.0*v3);
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rdDxMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 2 || i > ipField.iEnd - 1)
			{
				wenoXMinus(i, j) = ipField.dxMinusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;
				v2 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
				v3 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;

				WENO3rd(v1, v2, v3, wenoXMinus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rdDxPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (i < ipField.iStart + 1 || i > ipField.iEnd - 2)
			{
				wenoXPlus(i, j) = ipField.dxPlusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;
				v2 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
				v3 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;

				WENO3rd(v1, v2, v3, wenoXPlus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rdDyMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (j < ipField.jStart + 2 || j > ipField.jEnd - 1)
			{
				wenoYMinus(i, j) = ipField.dyMinusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;
				v2 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
				v3 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;

				WENO3rd(v1, v2, v3, wenoYMinus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO3rdDyPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{

			if (j < ipField.jStart + 1 || j > ipField.jEnd - 2)
			{
				wenoYPlus(i, j) = ipField.dyPlusPhi(i, j);
			}
			else
			{
				v1 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;
				v2 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
				v3 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;

				WENO3rd(v1, v2, v3, wenoYPlus(i, j));
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::WENO5thDerivation(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus)
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
inline void AdvectionMethod2D<TT>::WENO5thDxMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus)
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
inline void AdvectionMethod2D<TT>::WENO5thDxPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus)
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
inline void AdvectionMethod2D<TT>::WENO5thDyMinus(const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus)
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
inline void AdvectionMethod2D<TT>::WENO5thDyPlus(const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus)
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
	levelSet.phi.SaveOld();
	Array2D<TT>& OldLS = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;


	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -sign(OldLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), OldLS(i, j));
			levelSet(i, j) = OldLS(i, j) + k1(i, j);
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -sign(OldLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), OldLS(i, j));
			levelSet(i, j) = 3.0 / 4.0*OldLS(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
		}
	}

	WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -sign(OldLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), OldLS(i, j));
			levelSet(i, j) = 1.0 / 3.0*OldLS(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
		}
	}

}

template<class TT>
inline void AdvectionMethod2D<TT>::LSReinitializationTVDRK3(LS & levelSet, const double & dt, const int & iteration)
{

	levelSet.phi.SaveOld();
	Array2D<TT>& OldLS = levelSet.phi.dataArrayOld;
	Array2D<TT> originLS = levelSet.phi.dataArray;
	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	for (int k = 0; k < iteration; k++)
	{
		levelSet.phi.SaveOld();

		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				k1(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
				levelSet(i, j) = OldLS(i, j) + k1(i, j);
			}
		}

		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				k2(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
				levelSet(i, j) = 3.0 / 4.0*OldLS(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
			}
		}

		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				k3(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
				levelSet(i, j) = 1.0 / 3.0*OldLS(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
			}
		}
	}

}

template<class TT>
inline void AdvectionMethod2D<TT>::LSReinitializationTVDRK3(LS & levelSet, LS & originLS, const double & dt, const int& spatialOrder)
{
	levelSet.phi.SaveOld();
	Array2D<TT>& OldLS = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;


	if (spatialOrder == 3)
	{
		WENO3rdDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else
	{
		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
			levelSet(i, j) = OldLS(i, j) + k1(i, j);
		}
	}

	if (spatialOrder == 3)
	{
		WENO3rdDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else
	{
		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
			levelSet(i, j) = 3.0 / 4.0*OldLS(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
		}
	}

	if (spatialOrder == 3)
	{
		WENO3rdDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else
	{
		WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = -sign(originLS(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLS(i, j));
			levelSet(i, j) = 1.0 / 3.0*OldLS(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSOneSidedDerivativesSubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
{
	LSOneSidedDerivativesDxSubcellFixSecondOrder(levelSet, ipField, DxMinus, DxPlus);
	LSOneSidedDerivativesDySubcellFixSecondOrder(levelSet, ipField, DyMinus, DyPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSOneSidedDerivativesDxSubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus)
{
	int iStart = ipField.grid.iStart;
	int iEnd = ipField.grid.iEnd;
	double dx = ipField.dx;
	double dx2 = dx*dx;
	double thresL = pow(10.0, -10.0);
	double thresU = 2 / dx;

#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{

			int iLL = max(i - 2, iStart), iL = max(i - 1, iStart), iR = min(i + 1, iEnd), iRR = min(i + 2, iEnd);
			double lsLL = ipField(iLL, j), lsL = ipField(iL, j), ls = ipField(i, j), lsR = ipField(iR, j), lsRR = ipField(iRR, j);
			double d2L = ipField.dxxPhi(iL, j), d2C = ipField.dxxPhi(i, j), d2R = ipField.dxxPhi(iR, j);

			if (i == iStart || i == iEnd)
			{
				DxMinus(i, j) = ipField.dxPhi(i, j) + 0.5*dx*d2C;
				DxPlus(i, j) = ipField.dxPhi(i, j) - 0.5*dx*d2C;
				continue;
			}

			double dxL = dx;
			double phixx = MinMod(d2C*dx2, d2L*dx2);
			if (ls*lsL < 0 && abs(ls)>10 * DBL_EPSILON)
			{
				if (abs(phixx) < thresL)
				{
					dxL = dx*(ls) / (ls - lsL);
				}
				else
				{
					double D = (0.5*phixx - ls - lsL)*(0.5*phixx - ls - lsL) - 4 * ls*lsL;
					dxL = dx*(0.5 + (ls - lsL - sign(ls - lsL)*sqrt(D)) / phixx);
				}
				DxMinus(i, j) = ls / dxL + 0.5*dx * phixx / dx2;
			}
			else
			{
				DxMinus(i, j) = (ls - lsL) / dx + 0.5*dx * phixx / dx2;
			}

			double dxR = dx;
			phixx = MinMod(d2C*dx2, d2R*dx2);
			if (ls*lsR < 0 && abs(ls)>10 * DBL_EPSILON)
			{
				if (abs(phixx) < thresL)
				{
					dxR = dx*(ls) / (ls - lsR);
				}
				else
				{
					double D = (0.5*phixx - ls - lsR)*(0.5*phixx - ls - lsR) - 4 * ls*lsR;
					dxR = dx*(0.5 + (ls - lsR - sign(ls - lsR)*sqrt(D)) / phixx);
				}
				DxPlus(i, j) = -ls / dxR - 0.5*dx * phixx / dx2;
			}
			else
			{
				DxPlus(i, j) = (lsR - ls) / dx - 0.5*dx * phixx / dx2;
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSOneSidedDerivativesDySubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
{
	int jStart = ipField.grid.jStart;
	int jEnd = ipField.grid.jEnd;
	double dy = ipField.dy;
	double dy2 = dy*dy;
	double thresL = pow(10.0, -10.0);
	double thresU = 2 / dy;

#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{

			int jBB = max(j - 2, jStart), jB = max(j - 1, jStart), jT = min(j + 1, jEnd), jTT = min(j + 2, jEnd);
			double lsBB = ipField(i, jBB), lsB = ipField(i, jB), ls = ipField(i, j), lsT = ipField(i, jT), lsTT = ipField(i, jTT);
			double d2B = ipField.dyyPhi(i, jB), d2C = ipField.dyyPhi(i, j), d2T = ipField.dyPhi(i, jT);

			if (j == jStart || j == jEnd)
			{
				DyMinus(i, j) = ipField.dyPhi(i, j) + 0.5*dy*d2C;
				DyPlus(i, j) = ipField.dyPhi(i, j) - 0.5*dy*d2C;
				continue;
			}

			double dyB = dy;
			double phiyy = MinMod(d2C*dy2, d2B*dy2);
			if (ls*lsB < 0 && abs(ls)>10 * DBL_EPSILON)
			{
				if (abs(phiyy) < thresL)
				{
					dyB = dy*(ls) / (ls - lsB);
				}
				else
				{
					double D = (0.5*phiyy - ls - lsB)*(0.5*phiyy - ls - lsB) - 4 * ls*lsB;
					dyB = dy*(0.5 + (ls - lsB - sign(ls - lsB)*sqrt(D)) / phiyy);
				}
				DyMinus(i, j) = ls / dyB + 0.5*dy * phiyy / dy2;
			}
			else
			{
				DyMinus(i, j) = (ls - lsB) / dy + 0.5*dy * phiyy / dy2;
			}

			double dyT = dy;
			phiyy = MinMod(d2C*dy2, d2T*dy2);
			if (ls*lsT < 0 && abs(ls)>10 * DBL_EPSILON)
			{
				if (abs(phiyy) < thresL)
				{
					dyT = dy*(ls) / (ls - lsT);
				}
				else
				{
					double D = (0.5*phiyy - ls - lsT)*(0.5*phiyy - ls - lsT) - 4 * ls*lsT;
					dyT = dy*(0.5 + (ls - lsT - sign(ls - lsT)*sqrt(D)) / phiyy);
				}
				DyPlus(i, j) = -ls / dyT - 0.5*dy * phiyy / dy2;
			}
			else
			{
				DyPlus(i, j) = (lsT - ls) / dy - 0.5*dy * phiyy / dy2;
			}
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LSReinitializationTVDRK3SubcellFixSecondOrder(LS & levelSet, const double & dt, const int & iter)
{

	Array2D<TT>& DxMinus = levelSet.phi.dfdxM;
	Array2D<TT>& DxPlus = levelSet.phi.dfdxP;
	Array2D<TT>& DyMinus = levelSet.phi.dfdyM;
	Array2D<TT>& DyPlus = levelSet.phi.dfdyP;

	Array2D<TT> origin = levelSet.phi.dataArray;
	Array2D<TT>& oldLevelSet = levelSet.phi.dataArrayOld;
	for (int l = 1; l <= iter; l++)
	{
		levelSet.phi.SaveOld();

		AdvectionMethod2D<double>::LSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				double originL = origin(i, j);

				double k1 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
				levelSet(i, j) = oldLevelSet(i, j) + k1;
			}
		}

		AdvectionMethod2D<double>::LSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				double originL = origin(i, j);

				double k2 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
				levelSet(i, j) = 3.0 / 4.0*oldLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2);
			}
		}

		AdvectionMethod2D<double>::LSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
		{
			for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
			{
				double originL = origin(i, j);

				double k3 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
				levelSet(i, j) = 1.0 / 3.0*oldLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3);
			}
		}
		//levelSet.phi.Variable("phi");
		//MATLAB.Command("plot(0,0);hold on,contour(X, Y, phi0,[0 0],'b');contour(X, Y, phi, [0 0],'r');grid on,axis equal");
		//MATLAB.Command("hold off");
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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;



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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;



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
	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;
	int i, j;
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k1(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);

	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k2(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k3(i, j) = -dt*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), 1);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3(LS & levelSet, const FD& velocity, const double & dt)
{
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;
	int i, j;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k1(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k2(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));

	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k3(i, j) = -dt*velocity(i, j)*PropagatingGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), -1);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));

	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3(LS & levelSet, const FD& velocityX, const FD& velocityY, const double & dt)
{
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);

	double tempDxPhi, tempDyPhi;
	int i, j;
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

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

		k1(i, j) = -dt*(velocityX(i, j)*tempDxPhi + velocityY(i, j)*tempDyPhi);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

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
		k2(i, j) = -dt*(velocityX(i, j)*tempDxPhi + velocityY(i, j)*tempDyPhi);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
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
		k3(i, j) = -dt*(velocityX(i, j)*tempDxPhi + velocityY(i, j)*tempDyPhi);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3MACGrid(LS & levelSet, const FD & velocityX, const FD & velocityY, const double & dt)
{
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);

	double tempDxPhi, tempDyPhi, velX, velY;
	int i, j;
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}
		
		k1(i, j) = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1(i, j);
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k2(i, j) = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2(i, j));
	}

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k3(i, j) = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3(i, j));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSPropagatingTVDRK3MACGrid(LS & levelSet, const FD & velocityX, const FD & velocityY, const double & dt, const int & spatialOrder)
{
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	//Array2D<TT>& k1 = levelSet.phi.K1;
	//Array2D<TT>& k2 = levelSet.phi.K2;
	//Array2D<TT>& k3 = levelSet.phi.K3;
	double k1, k2, k3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	if (spatialOrder == 3)
	{
		LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else if (spatialOrder == 5)
	{
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}

	double tempDxPhi, tempDyPhi, velX, velY;
	int i, j;
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY, k1)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k1 = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = originLevelSet(i, j) + levelSet.Cutoff(i, j)*k1;
	}

	if (spatialOrder == 3)
	{
		LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else if (spatialOrder == 5)
	{
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY, k2)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k2 = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k2);
	}

	if (spatialOrder == 3)
	{
		LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else if (spatialOrder == 5)
	{
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for private(i, j, tempDxPhi, tempDyPhi, velX, velY, k3)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		velX = 0.5*(velocityX(i + 1, j) + velocityX(i, j));
		velY = 0.5*(velocityY(i, j + 1) + velocityY(i, j));
		if (velX >= 0)
		{
			tempDxPhi = wenoXMinus(i, j);
		}
		else
		{
			tempDxPhi = wenoXPlus(i, j);
		}
		if (velY >= 0)
		{
			tempDyPhi = wenoYMinus(i, j);
		}
		else
		{
			tempDyPhi = wenoYPlus(i, j);
		}

		k3  = -dt*(velX*tempDxPhi + velY*tempDyPhi);
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + levelSet.Cutoff(i, j)*k3);
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO3rdDerivation(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus)
{
	LLSWENO3rdDxMinus(levelSet, ipField, wenoXMinus);
	LLSWENO3rdDxPlus(levelSet, ipField, wenoXPlus);
	LLSWENO3rdDyMinus(levelSet, ipField, wenoYMinus);
	LLSWENO3rdDyPlus(levelSet, ipField, wenoYPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO3rdDxMinus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		if (i < ipField.iStart + 2 || i > ipField.iEnd - 1)
		{
			wenoXMinus(i, j) = ipField.dxMinusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i - 1, j) - ipField(i - 2, j))*ipField.oneOverdx;
			v2 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;
			v3 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;

			WENO3rd(v1, v2, v3, wenoXMinus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO3rdDxPlus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

		if (i < ipField.iStart + 1 || i > ipField.iEnd - 2)
		{
			wenoXPlus(i, j) = ipField.dxPlusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i + 2, j) - ipField(i + 1, j))*ipField.oneOverdx;
			v2 = (ipField(i + 1, j) - ipField(i, j))*ipField.oneOverdx;
			v3 = (ipField(i, j) - ipField(i - 1, j))*ipField.oneOverdx;

			WENO3rd(v1, v2, v3, wenoXPlus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO3rdDyMinus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

		if (j < ipField.jStart + 2 || j > ipField.jEnd - 1)
		{
			wenoYMinus(i, j) = ipField.dyMinusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i, j - 1) - ipField(i, j - 2))*ipField.oneOverdy;
			v2 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;
			v3 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;

			WENO3rd(v1, v2, v3, wenoYMinus(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO3rdDyPlus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus)
{
	TT v1, v2, v3;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

		if (j < ipField.jStart + 1 || j > ipField.jEnd - 2)
		{
			wenoYPlus(i, j) = ipField.dyPlusPhi(i, j);
		}
		else
		{
			v1 = (ipField(i, j + 2) - ipField(i, j + 1))*ipField.oneOverdy;
			v2 = (ipField(i, j + 1) - ipField(i, j))*ipField.oneOverdy;
			v3 = (ipField(i, j) - ipField(i, j - 1))*ipField.oneOverdy;

			WENO3rd(v1, v2, v3, wenoYPlus(i, j));
		}
	}
}


template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDerivation(const LS& levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus, Array2D<TT>& wenoXPlus, Array2D<TT>& wenoYMinus, Array2D<TT>& wenoYPlus)
{
	LLSWENO5thDxMinus(levelSet, ipField, wenoXMinus);
	LLSWENO5thDxPlus(levelSet, ipField, wenoXPlus);
	LLSWENO5thDyMinus(levelSet, ipField, wenoYMinus);
	LLSWENO5thDyPlus(levelSet, ipField, wenoYPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSWENO5thDxMinus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXMinus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
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
inline void AdvectionMethod2D<TT>::LLSWENO5thDxPlus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoXPlus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

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
inline void AdvectionMethod2D<TT>::LLSWENO5thDyMinus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYMinus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

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
inline void AdvectionMethod2D<TT>::LLSWENO5thDyPlus(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& wenoYPlus)
{
	TT v1, v2, v3, v4, v5;
	int i, j;
#pragma omp parallel for private(i, j, v1, v2, v3, v4, v5)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;

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
	levelSet.phi.SaveOld();
	Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	int i, j;

	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k1(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k2(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
	}
	LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		k3(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
		levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3(LS & levelSet, const double & dt, const int & iter)
{
	Array2D<TT>& k1 = levelSet.phi.K1;
	Array2D<TT>& k2 = levelSet.phi.K2;
	Array2D<TT>& k3 = levelSet.phi.K3;

	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	int i, j;
	for (int l = 1; l <= iter; l++)
	{
		levelSet.phi.SaveOld();
		Array2D<TT>& originLevelSet = levelSet.phi.dataArrayOld;

		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			k1(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			k2(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2(i, j));
		}
		LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(i,j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			k3(i, j) = -sign(originLevelSet(i, j))*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originLevelSet(i, j));
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3(i, j));
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3(LS & levelSet, const double & dt, const int & iter, const int & spatialOrder)
{
	//Array2D<TT>& k1 = levelSet.phi.K1;
	//Array2D<TT>& k2 = levelSet.phi.K2;
	//Array2D<TT>& k3 = levelSet.phi.K3;
	double dx = levelSet.grid.dx;
	double k1, k2, k3;
	double originL;
	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;
	
	Array2D<TT> origin = levelSet.phi.dataArray;
	Array2D<TT>& oldLevelSet = levelSet.phi.dataArrayOld;
	int i, j;
	for (int l = 1; l <= iter; l++)
	{
		levelSet.phi.SaveOld();
		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i, j, k1, originL)
		for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			originL = origin(i, j);
			k1 = -sign(originL, dx)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			levelSet(i, j) = oldLevelSet(i, j) + k1;
		}

		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i,j, k2, originL)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			originL = origin(i, j);
			k2  = -sign(originL, dx)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			levelSet(i, j) = 3.0 / 4.0*oldLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2);
		}

		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i,j, k3, originL)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			originL = origin(i, j);
			k3 = -sign(originL, dx)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			levelSet(i, j) = 1.0 / 3.0*oldLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3);
		}
	}
}

// Peter Smereka, Subcell fixed reinitialization
template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3usingSubcellFix(LS & levelSet, const double & dt, const int & iter, const int & spatialOrder)
{
	double k1, k2, k3;
	Array2D<TT>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<TT>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<TT>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<TT>& wenoYPlus = levelSet.phi.dfdyP;

	int iStart = levelSet.grid.iStart, iEnd = levelSet.grid.iEnd, jStart = levelSet.grid.jStart, jEnd = levelSet.grid.jEnd;
	double dx = levelSet.grid.dx, dy = levelSet.grid.dy;
	double oneOverdx = levelSet.grid.oneOverdx, oneOverdy = levelSet.grid.oneOverdy;
	int i, j;
	int iL, iR, jB, jT;
	double Dij, lsL, lsR, lsB, lsT, lsC;
	double originL;
	double D1, D2;
	double a, b, c, d, theta;
	Array2D<TT> DDD(levelSet.grid);
	Array2D<TT> origin = levelSet.phi.dataArray;
	Array2D<TT> currentLevelSet(levelSet.grid);
	Array2D<TT>& oldLevelSet = levelSet.phi.dataArrayOld;
	for (int l = 1; l <= iter; l++)
	{
		levelSet.phi.SaveOld();
		currentLevelSet = levelSet.phi.dataArray;
		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i, j, k1, iL, iR, jB, jT, Dij, D1, D2, lsL, lsR, lsB, lsT, lsC, originL)
		for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;

			iL = max(i - 1, iStart), iR = min(i + 1, iEnd);
			jB = max(j - 1, jStart), jT = min(j + 1, jEnd);
			lsL = currentLevelSet(iL, j), lsR = currentLevelSet(iR, j);
			lsB = currentLevelSet(i, jB), lsT = currentLevelSet(i, jT);
			lsC = currentLevelSet(i, j);

			originL = origin(i, j);
			k1 = 0, D1 = 0, D2 = 0, Dij = 0;

			if (lsL*lsC < 0 || lsR*lsC < 0 || lsB*lsC < 0 || lsT*lsC < 0)
			{
				//if (lsL*lsR < 0)	D1 = abs(lsR - lsL) / 2;
				//else				D1 = max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON);
				//if (lsB*lsT < 0)	D2 = abs(lsT - lsB) / 2;
				//else				D2 = max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON);
				D1 = max(abs(lsR - lsL) / 2, max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON));
				D2 = max(abs(lsT - lsB) / 2, max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON));
				Dij = lsC / sqrt(D1*D1 + D2*D2);

				//k1 = -dt*(sign(originL)*abs(levelSet(i, j)) * oneOverdx - Dij);
				k1 = -dt*(sign(originL)*abs(levelSet(i, j)) / sqrt(D1*D1 + D2*D2) - Dij);
				k1 = 0;
			}
			else
			{
				k1 = -sign(originL)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			}
			levelSet(i, j) = oldLevelSet(i, j) + k1;
			
		}
	

		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}

		currentLevelSet = levelSet.phi.dataArray;
#pragma omp parallel for private(i, j, k2, iL, iR, jB, jT, Dij, D1, D2, lsL, lsR, lsB, lsT, lsC, originL)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;

			iL = max(i - 1, iStart), iR = min(i + 1, iEnd);
			jB = max(j - 1, jStart), jT = min(j + 1, jEnd);
			lsL = currentLevelSet(iL, j), lsR = currentLevelSet(iR, j);
			lsB = currentLevelSet(i, jB), lsT = currentLevelSet(i, jT);
			lsC = currentLevelSet(i, j);

			originL = origin(i, j);
			k2 = 0, D1 = 0, D2 = 0, Dij = 0;

			//if ((lsL*lsC < 0 || lsR*lsC < 0 || lsB*lsC < 0 || lsT*lsC < 0) && (lsL*lsR < 0 || lsB*lsT < 0))
			if (lsL*lsC < 0 || lsR*lsC < 0 || lsB*lsC < 0 || lsT*lsC < 0)
			{
				//if (lsL*lsR < 0)	D1 = abs(lsR - lsL) / 2;
				//else				D1 = max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON);
				//if (lsB*lsT < 0)	D2 = abs(lsT - lsB) / 2;
				//else				D2 = max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON);
				D1 = max(abs(lsR - lsL) / 2, max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON));
				D2 = max(abs(lsT - lsB) / 2, max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON));
				Dij = lsC / sqrt(D1*D1 + D2*D2);

				//k2 = -dt*(sign(originL, dx)*abs(levelSet(i, j)) * oneOverdx - Dij);
				k2 = -dt*(sign(originL)*abs(levelSet(i, j)) / sqrt(D1*D1 + D2*D2) - Dij);
			}
			else
			{
				k2 = -sign(originL, dx)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			}
			levelSet(i, j) = 3.0 / 4.0*oldLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2);
		}
		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(levelSet, levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}

		currentLevelSet = levelSet.phi.dataArray;
#pragma omp parallel for private(i, j, k3, iL, iR, jB, jT, Dij, D1, D2, lsL, lsR, lsB, lsT, lsC, originL)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			i = levelSet.tubeIndex(k).i;
			j = levelSet.tubeIndex(k).j;
			
			iL = max(i - 1, iStart), iR = min(i + 1, iEnd);
			jB = max(j - 1, jStart), jT = min(j + 1, jEnd);
			lsL = currentLevelSet(iL, j), lsR = currentLevelSet(iR, j);
			lsB = currentLevelSet(i, jB), lsT = currentLevelSet(i, jT);
			lsC = currentLevelSet(i, j);

			originL = origin(i, j);
			k3 = 0, D1 = 0, D2 = 0, Dij = 0;

			if (lsL*lsC < 0 || lsR*lsC < 0 || lsB*lsC < 0 || lsT*lsC < 0)
			{
				//if (lsL*lsR < 0)	D1 = abs(lsR - lsL) / 2;
				//else				D1 = max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON);
				//if (lsB*lsT < 0)	D2 = abs(lsT - lsB) / 2;
				//else				D2 = max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON);
				D1 = max(abs(lsR - lsL) / 2, max(max(abs(lsL - lsC), abs(lsR - lsC)), DBL_EPSILON));
				D2 = max(abs(lsT - lsB) / 2, max(max(abs(lsB - lsC), abs(lsT - lsC)), DBL_EPSILON));
				Dij = lsC / sqrt(D1*D1 + D2*D2);

				//k3 = -dt*(sign(originL, dx)*abs(levelSet(i, j)) * oneOverdx - Dij);
				k3 = -dt*(sign(originL)*abs(levelSet(i, j)) / sqrt(D1*D1 + D2*D2) - Dij);
			}
			else
			{
				k3 = -sign(originL, dx)*dt*ReinitialGodunov(wenoXPlus(i, j), wenoXMinus(i, j), wenoYPlus(i, j), wenoYMinus(i, j), originL);
			}
			levelSet(i, j) = 1.0 / 3.0*oldLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3);
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSOneSidedDerivativesSubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
{
	LLSOneSidedDerivativesDxSubcellFixSecondOrder(levelSet, ipField, DxMinus, DxPlus);
	LLSOneSidedDerivativesDySubcellFixSecondOrder(levelSet, ipField, DyMinus, DyPlus);
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSOneSidedDerivativesDxSubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus)
{
	int iStart = ipField.grid.iStart;
	int iEnd = ipField.grid.iEnd;
	double dx = ipField.dx;
	double dx2 = dx*dx;
	double thresL = pow(10.0, -10.0);
	double thresU = 2 / dx;

#pragma omp parallel for
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		int i = levelSet.tubeIndex(k).i;
		int j = levelSet.tubeIndex(k).j;

		int iLL = max(i - 2, iStart), iL = max(i - 1, iStart), iR = min(i + 1, iEnd), iRR = min(i + 2, iEnd);
		double lsLL = ipField(iLL, j), lsL = ipField(iL, j), ls = ipField(i, j), lsR = ipField(iR, j), lsRR = ipField(iRR, j);
		double d2L = ipField.dxxPhi(iL, j), d2C = ipField.dxxPhi(i, j), d2R = ipField.dxxPhi(iR, j);

		if (i == iStart || i == iEnd)
		{
			DxMinus(i, j) = ipField.dxPhi(i, j) + 0.5*dx*d2C;
			DxPlus(i, j) = ipField.dxPhi(i, j) - 0.5*dx*d2C;
			continue;
		}

		double dxL = dx;
		double phixx = MinMod(d2C*dx2, d2L*dx2);
		if (ls*lsL < 0 && abs(ls)>10*DBL_EPSILON)
		{
			if (abs(phixx) < thresL)
			{
				dxL = dx*(ls) / (ls - lsL);
			}
			else
			{
				double D = (0.5*phixx - ls - lsL)*(0.5*phixx - ls - lsL) - 4 * ls*lsL;
				dxL = dx*(0.5 + (ls - lsL - sign(ls - lsL)*sqrt(D)) / phixx);
			}
			DxMinus(i, j) = ls / dxL + 0.5*dx * phixx / dx2;
		}
		else
		{
			DxMinus(i, j) = (ls - lsL) / dx + 0.5*dx * phixx / dx2;
		}

		double dxR = dx;
		phixx = MinMod(d2C*dx2, d2R*dx2);
		if (ls*lsR < 0 && abs(ls)>10 * DBL_EPSILON)
		{
			if (abs(phixx) < thresL)
			{
				dxR = dx*(ls) / (ls - lsR);
			}
			else
			{
				double D = (0.5*phixx - ls - lsR)*(0.5*phixx - ls - lsR) - 4 * ls*lsR;
				dxR = dx*(0.5 + (ls - lsR - sign(ls - lsR)*sqrt(D)) / phixx);
			}
			DxPlus(i, j) = -ls / dxR - 0.5*dx * phixx / dx2;
		}
		else
		{
			DxPlus(i, j) = (lsR - ls) / dx - 0.5*dx * phixx / dx2;
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSOneSidedDerivativesDySubcellFixSecondOrder(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
{
	int jStart = ipField.grid.jStart;
	int jEnd = ipField.grid.jEnd;
	double dy = ipField.dy;
	double dy2 = dy*dy;
	double thresL = pow(10.0, -10.0);
	double thresU = 2 / dy;

#pragma omp parallel for
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		int i = levelSet.tubeIndex(k).i;
		int j = levelSet.tubeIndex(k).j;

		int jBB = max(j - 2, jStart), jB = max(j - 1, jStart), jT = min(j + 1, jEnd), jTT = min(j + 2, jEnd);
		double lsBB = ipField(i, jBB), lsB = ipField(i, jB), ls = ipField(i, j), lsT = ipField(i, jT), lsTT = ipField(i, jTT);
		double d2B = ipField.dyyPhi(i, jB), d2C = ipField.dyyPhi(i, j), d2T = ipField.dyPhi(i, jT);

		if (j == jStart || j == jEnd)
		{
			DyMinus(i, j) = ipField.dyPhi(i, j) + 0.5*dy*d2C;
			DyPlus(i, j) = ipField.dyPhi(i, j) - 0.5*dy*d2C;
			continue;
		}

		double dyB = dy;
		double phiyy = MinMod(d2C*dy2, d2B*dy2);
		if (ls*lsB < 0 && abs(ls)>10 * DBL_EPSILON)
		{
			if (abs(phiyy) < thresL)
			{
				dyB = dy*(ls) / (ls - lsB);
			}
			else
			{
				double D = (0.5*phiyy - ls - lsB)*(0.5*phiyy - ls - lsB) - 4 * ls*lsB;
				dyB = dy*(0.5 + (ls - lsB - sign(ls - lsB)*sqrt(D)) / phiyy);
			}
			DyMinus(i, j) = ls / dyB + 0.5*dy * phiyy / dy2;
		}
		else
		{
			DyMinus(i, j) = (ls - lsB) / dy + 0.5*dy * phiyy / dy2;
		}

		double dyT = dy;
		phiyy = MinMod(d2C*dy2, d2T*dy2);
		if (ls*lsT < 0 && abs(ls)>10 * DBL_EPSILON)
		{
			if (abs(phiyy) < thresL)
			{
				dyT = dy*(ls) / (ls - lsT);
			}
			else
			{
				double D = (0.5*phiyy - ls - lsT)*(0.5*phiyy - ls - lsT) - 4 * ls*lsT;
				dyT = dy*(0.5 + (ls - lsT - sign(ls - lsT)*sqrt(D)) / phiyy);
			}
			DyPlus(i, j) = -ls / dyT - 0.5*dy * phiyy / dy2;
		}
		else
		{
			DyPlus(i, j) = (lsT - ls) / dy - 0.5*dy * phiyy / dy2;
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3SubcellFixSecondOrder(LS & levelSet, const double & dt, const int & iter)
{
	Array2D<TT>& DxMinus = levelSet.phi.dfdxM;
	Array2D<TT>& DxPlus = levelSet.phi.dfdxP;
	Array2D<TT>& DyMinus = levelSet.phi.dfdyM;
	Array2D<TT>& DyPlus = levelSet.phi.dfdyP;

	Array2D<TT> origin = levelSet.phi.dataArray;
	Array2D<TT>& oldLevelSet = levelSet.phi.dataArrayOld;
	for (int l = 1; l <= iter; l++)
	{
		levelSet.phi.SaveOld();

		AdvectionMethod2D<double>::LLSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
		{
			int i = levelSet.tubeIndex(k).i;
			int j = levelSet.tubeIndex(k).j;
			double originL = origin(i, j);

			double k1 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
			levelSet(i, j) = oldLevelSet(i, j) + k1;
		}

		AdvectionMethod2D<double>::LLSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			int i = levelSet.tubeIndex(k).i;
			int j = levelSet.tubeIndex(k).j;
			double originL = origin(i, j);

			double k2 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
			levelSet(i, j) = 3.0 / 4.0*oldLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2);
		}

		AdvectionMethod2D<double>::LLSOneSidedDerivativesSubcellFixSecondOrder(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
#pragma omp parallel for
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			int i = levelSet.tubeIndex(k).i;
			int j = levelSet.tubeIndex(k).j;
			double originL = origin(i, j);

			double k3 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
			levelSet(i, j) = 1.0 / 3.0*oldLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3);
		}
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::OneSideLeftDerivatives(const double & lsL, const double & lsC, const double & valueL, const double & valueC, const double & dx, double & val, double thetadx)
{
	if (lsL*lsC >= 0)
	{
		thetadx = dx;
		val = (valueC - valueL) / thetadx;
	}
	else
	{
		thetadx = dx*abs(lsC) / (abs(lsC) + abs(lsL));
		val = (valueC - 0) / thetadx;
	}
}

template<class TT>
inline void AdvectionMethod2D<TT>::OneSideRightDerivatives(const double & lsC, const double & lsR, const double & valueC, const double & valueR, const double & dx, double & val, double thetadx)
{
	if (lsR*lsC >= 0)
	{
		thetadx = dx;
		val = (valueR - valueC) / thetadx;
	}
	else
	{
		thetadx = dx*abs(lsC) / (abs(lsC) + abs(lsR));
		val = (0 - valueC) / thetadx;
	}
}

//template<class TT>
//inline void AdvectionMethod2D<TT>::LLSOneSidedDerivatives(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
//{
//	LLSOneSidedDerivativesDx(levelSet, ipField, DxMinus, DxPlus);
//	LLSOneSidedDerivativesDy(levelSet, ipField, DyMinus, DyPlus);
//}
//
//template<class TT>
//inline void AdvectionMethod2D<TT>::LLSOneSidedDerivativesDx(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DxMinus, Array2D<TT>& DxPlus)
//{
//	int computedDomaion = 1;
//	int iStart = ipField.grid.iStart;
//	int iEnd = ipField.grid.iEnd;
//	double dx = ipField.grid.dx;
//	double oneOverdx = ipField.grid.oneOverdx;
//
//#pragma omp parallel for
//	for (int k = 1; k <= levelSet.numTube; k++)
//	{
//		int i = levelSet.tubeIndex(k).i;
//		int j = levelSet.tubeIndex(k).j;
//
//		int iLLL = max(i - 3, iStart), iLL = max(i - 2, iStart), iL = max(i - 1, iStart);
//		int iRRR = min(i + 3, iEnd), iRR = min(i + 2, iEnd), iR = min(i + 1, iEnd);
//		double dx1,dx2,dx3,dx
//		double dfdx1M3, dfdx1M2, dfdx1M1, dfdx1P1, dfdx1P2;
//		double dfdx2M2, dfdx2M1, dfdx2, dfdx2P1;
//		double dfdx3M2, dfdx3M1, dfdx3P1;
//
//		double D2M = D2(iL, j);
//		double D20 = D2(i, j);
//		if (i == iStart || i == iEnd)
//		{
//			DxMinus(i, j) = ipField.dxMinusPhi(i, j);
//		}
//		else
//		{
//			if (levelSet.tube(i,j) <= computedDomaion)
//			{
//				double a;
//				if (abs(D2M) < abs(D20))
//				{
//					a = 2 * dx*dx*MinAbs(D3(iLL, j), D3(iL, j));
//				}
//				else
//				{
//					a = -dx*dx*MinAbs(D3(iL, j), D3(i, j));
//				}
//				DxMinus(i, j) = ipField.dxMinusPhi(i, j) + MinMod(D2M, D20)*dx + a;
//			}
//		}
//		
//		double D2P = D2(iR, j);
//		if (i == iStart || i == iEnd)
//		{
//			DxPlus(i, j) = ipField.dxPlusPhi(i, j);
//		}
//		else
//		{
//			if (levelSet.tube(i, j) <= computedDomaion)
//			{
//				double b;
//				if (abs(D20) < abs(D2P))
//				{
//					b = -dx*dx*MinAbs(D3(iL, j), D3(i, j));
//				}
//				else
//				{
//					b = 2 * dx*dx*MinAbs(D3(i, j), D3(iR, j));
//				}
//				DxPlus(i, j) = ipField.dxPlusPhi(i, j) - MinMod(D20, D2P)*dx + b;
//			}
//
//		}
//	}
//
//}
//
//template<class TT>
//inline void AdvectionMethod2D<TT>::LLSOneSidedDerivativesDy(const LS & levelSet, const Field2D<TT>& ipField, Array2D<TT>& DyMinus, Array2D<TT>& DyPlus)
//{
//	int computedDomaion = 1;
//	int jStart = ipField.grid.jStart;
//	int jEnd = ipField.grid.jEnd;
//	double dy = ipField.grid.dy;
//	double oneOverdy = ipField.grid.oneOverdy;
//
//#pragma omp parallel for
//	for (int k = 1; k <= levelSet.numTube; k++)
//	{
//		int i = levelSet.tubeIndex(k).i;
//		int j = levelSet.tubeIndex(k).j;
//
//		int jBB = j - 2;
//		int jB = j - 1;
//		double D2M = D2(i, jB);
//		double D20 = D2(i, j);
//		if (j < jStart + 2 || j > jEnd - 1)
//		{
//			DyMinus(i, j) = ipField.dyMinusPhi(i, j) + MinMod(D2M, D20)*dy;
//		}
//		else
//		{
//			if (levelSet.tube(i, j) <= computedDomaion)
//			{
//
//				double a;
//				if (abs(D2M) < abs(D20))
//				{
//					a = 2 * dy*dy*MinAbs(D3(i, jBB), D3(i, jB));
//				}
//				else
//				{
//					a = -dy*dy*MinAbs(D3(i, jB), D3(i, j));
//				}
//				DyMinus(i, j) = ipField.dyMinusPhi(i, j) + MinMod(D2M, D20)*dy + a;
//			}
//		}
//
//		int jT = j + 1;
//		double D2P = D2(i, jT);
//		if (j < jStart + 1 || j > jEnd - 2)
//		{
//			DyPlus(i, j) = ipField.dyPlusPhi(i, j) - MinMod(D20, D2P)*dy;
//		}
//		else
//		{
//			if (levelSet.tube(i, j) <= computedDomaion)
//			{
//				double b;
//				if (abs(D20) < abs(D2P))
//				{
//					b = -dy*dy*MinAbs(D3(i, jB), D3(i, j));
//				}
//				else
//				{
//					b = 2 * dy*dy*MinAbs(D3(i, j), D3(i, jT));
//				}
//				DyPlus(i, j) = ipField.dyPlusPhi(i, j) - MinMod(D20, D2P)*dy + b;
//			}
//		}
//	}
//}
//
//// Chene, Min, Gibou, High-order reinitialization
//template<class TT>
//inline void AdvectionMethod2D<TT>::LLSReinitializationTVDRK3SubcellFixHighOrder(LS & levelSet, const double & dt, const int & iter)
//{
//	Array2D<TT>& DxMinus = levelSet.phi.dfdxM;
//	Array2D<TT>& DxPlus = levelSet.phi.dfdxP;
//	Array2D<TT>& DyMinus = levelSet.phi.dfdyM;
//	Array2D<TT>& DyPlus = levelSet.phi.dfdyP;
//	
//	Array2D<TT> origin = levelSet.phi.dataArray;
//	Array2D<TT>& oldLevelSet = levelSet.phi.dataArrayOld;
//	for (int l = 1; l <= iter; l++)
//	{
//		levelSet.phi.SaveOld();
//
//		AdvectionMethod2D<double>::LLSOneSidedDerivatives(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
//		DxMinus.Variable("DxM");
//		DxPlus.Variable("DxP");
//		DyMinus.Variable("DyM");
//		DyPlus.Variable("DyP");
//#pragma omp parallel for
//		for (int k = levelSet.tubeIndex.iStart; k <= levelSet.numTube; k++)
//		{
//			int i = levelSet.tubeIndex(k).i;
//			int j = levelSet.tubeIndex(k).j;
//			double originL = origin(i, j);
//
//			double k1 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
//			levelSet(i, j) = oldLevelSet(i, j) + k1;
//		}
//
//		AdvectionMethod2D<double>::LLSOneSidedDerivatives(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
//#pragma omp parallel for
//		for (int k = 1; k <= levelSet.numTube; k++)
//		{
//			int i = levelSet.tubeIndex(k).i;
//			int j = levelSet.tubeIndex(k).j;
//			double originL = origin(i, j);
//
//			double k2 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
//			levelSet(i, j) = 3.0 / 4.0*oldLevelSet(i, j) + 1.0 / 4.0*(levelSet(i, j) + k2);
//		}
//
//		AdvectionMethod2D<double>::LLSOneSidedDerivatives(levelSet, levelSet.phi, DxMinus, DxPlus, DyMinus, DyPlus);
//#pragma omp parallel for
//		for (int k = 1; k <= levelSet.numTube; k++)
//		{
//			int i = levelSet.tubeIndex(k).i;
//			int j = levelSet.tubeIndex(k).j;
//			double originL = origin(i, j);
//
//			double k3 = -sign(originL)*dt*ReinitialGodunov(DxPlus(i, j), DxMinus(i, j), DyPlus(i, j), DyMinus(i, j), originL);
//			levelSet(i, j) = 1.0 / 3.0*oldLevelSet(i, j) + 2.0 / 3.0*(levelSet(i, j) + k3);
//		}
//	}
//}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSQuantityExtension(LS & ipLS, FD & ipQuantity, const int & temporalOrder, const int & spatialOrder)
{
	ipQuantity.SaveOld();
	Array2D<double>& originQuantity = ipQuantity.dataArrayOld;

	//Array2D<double>& k1 = ipQuantity.K1;
	//Array2D<double>& k2 = ipQuantity.K2;
	//Array2D<double>& k3 = ipQuantity.K3;
	double k1, k2, k3;

	Array2D<double>& wenoXMinus = ipQuantity.dfdxM;
	Array2D<double>& wenoXPlus = ipQuantity.dfdxP;
	Array2D<double>& wenoYMinus = ipQuantity.dfdyM;
	Array2D<double>& wenoYPlus = ipQuantity.dfdyP;

	ipLS.LComputeUnitNormal();
	//ArrayVec2DVariable("unitnormal", ipLS.unitNormal.dataArray);

	VT normal;
	double signPhi;
	double tempDxPhi, tempDyPhi;
	int i, j;
	int updatedRegion = 2;
	double cflCondition = 0.2;
	double dt = cflCondition * min(ipLS.grid.dx, ipLS.grid.dy);

	if (spatialOrder == 3)
	{
		AdvectionMethod2D<double>::LLSWENO3rdDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
	else if (spatialOrder == 5)
	{
		AdvectionMethod2D<double>::LLSWENO5thDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
	}
#pragma omp parallel for private(i, j, normal, signPhi, tempDxPhi, tempDyPhi, k1)
	for (int k = 1; k <= ipLS.numTube; k++)
	{
		ipLS.TubeIndex(k, i, j);
		if (ipLS.tube(i, j) <= updatedRegion)
		{
			normal = ipLS.unitNormal(i, j);
			signPhi = AdvectionMethod2D<double>::sign(ipLS(i, j));
			if (signPhi*normal.i >= 0)
			{
				tempDxPhi = wenoXMinus(i, j);
			}
			else
			{
				tempDxPhi = wenoXPlus(i, j);
			}
			if (signPhi*normal.j >= 0)
			{
				tempDyPhi = wenoYMinus(i, j);
			}
			else
			{
				tempDyPhi = wenoYPlus(i, j);
			}
			k1 = -dt*signPhi*(normal.x*tempDxPhi + normal.y*tempDyPhi);
			ipQuantity(i, j) = originQuantity(i, j) + k1;
		}
	}

	if (temporalOrder == 3)
	{
		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i, j, normal, signPhi, tempDxPhi, tempDyPhi, k2)
		for (int k = 1; k <= ipLS.numTube; k++)
		{
			ipLS.TubeIndex(k, i, j);
			if (ipLS.tube(i, j) <= updatedRegion)
			{
				normal = ipLS.unitNormal(i, j);
				signPhi = AdvectionMethod2D<double>::sign(ipLS(i, j));
				if (signPhi*normal.i >= 0)
				{
					tempDxPhi = wenoXMinus(i, j);
				}
				else
				{
					tempDxPhi = wenoXPlus(i, j);
				}
				if (signPhi*normal.j >= 0)
				{
					tempDyPhi = wenoYMinus(i, j);
				}
				else
				{
					tempDyPhi = wenoYPlus(i, j);
				}
				k2 = -dt*signPhi*(normal.x*tempDxPhi + normal.y*tempDyPhi);
				ipQuantity(i, j) = 1.0 / 4.0 * (3 * originQuantity(i, j) + (ipQuantity(i, j) + k2));
			}
		}

		if (spatialOrder == 3)
		{
			AdvectionMethod2D<double>::LLSWENO3rdDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
		else if (spatialOrder == 5)
		{
			AdvectionMethod2D<double>::LLSWENO5thDerivation(ipLS, ipQuantity, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
		}
#pragma omp parallel for private(i, j, normal, signPhi, tempDxPhi, tempDyPhi, k3)
		for (int k = 1; k <= ipLS.numTube; k++)
		{
			ipLS.TubeIndex(k, i, j);
			if (ipLS.tube(i, j) <= updatedRegion)
			{
				normal = ipLS.unitNormal(i, j);
				signPhi = AdvectionMethod2D<double>::sign(ipLS(i, j));
				if (signPhi*normal.i >= 0)
				{
					tempDxPhi = wenoXMinus(i, j);
				}
				else
				{
					tempDxPhi = wenoXPlus(i, j);
				}
				if (signPhi*normal.j >= 0)
				{
					tempDyPhi = wenoYMinus(i, j);
				}
				else
				{
					tempDyPhi = wenoYPlus(i, j);
				}
				k3 = -dt*signPhi*(normal.x*tempDxPhi + normal.y*tempDyPhi);
				ipQuantity(i, j) = 1.0 / 3.0 * (originQuantity(i, j) + 2.0 *(ipQuantity(i, j) + k3));
			}
		}
	}
	
}

template<class TT>
inline void AdvectionMethod2D<TT>::LLSQuantityExtension(LS & ipLS, FD & ipQuantity, const int & temporalOrder, const int & spatialOrder, const int & iter)
{
	int i, j;
	int updatedRegion = 2;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= ipLS.numTube; k++)
	{
		ipLS.TubeIndex(k, i, j);
		if (ipLS.tube(i, j) > updatedRegion)
		{
			ipQuantity(i, j) = 0;
		}
	}

	for (int m = 1;  m <= iter;  m++)
	{
		LLSQuantityExtension(ipLS, ipQuantity, temporalOrder, spatialOrder);
	}
}


template<class TT>
inline double AdvectionMethod2D<TT>::AdaptiveTimeStep(const FD & velocity1, const double & cflCondition)
{
	double maxVel1 = 0;

	for (int i = velocity1.grid.iStart; i <= velocity1.grid.iEnd; i++)
	{
		for (int j = velocity1.grid.jStart; j <= velocity1.grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
		}
	}
	return cflCondition*max(velocity1.grid.dx, velocity1.grid.dy) / (maxVel1 + DBL_EPSILON);
}

template<class TT>
inline double AdvectionMethod2D<TT>::AdaptiveTimeStep(const FD & velocity1, const FD & velocity2, const double & cflCondition)
{
	double maxVel1 = 0;
	double maxVel2 = 0;

	for (int i = velocity1.grid.iStart; i <= velocity1.grid.iEnd; i++)
	{
		for (int j = velocity1.grid.jStart; j <= velocity1.grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
			if (abs(velocity2(i, j)) > maxVel2)
			{
				maxVel2 = abs(velocity2(i, j));
			}
		}
	}
	return cflCondition*min(velocity1.grid.dx / (maxVel1 + DBL_EPSILON),velocity1.grid.dy / (maxVel2 + DBL_EPSILON));
}