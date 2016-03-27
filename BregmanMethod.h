#pragma once

#include "CommonDef.h"
#include "Field2D.h"
#include "LevelSet2D.h"

template <class TT>
class BregmanMethod
{
public:
	BregmanMethod();
	~BregmanMethod();

	
	static Vector2D<TT> Shrink(const Vector2D<TT>& input, const TT& lambda);
	static Field2D<Vector2D<TT>> Shrink(const Field2D<Vector2D<TT>>& input, const TT& lambda);

	double testFunction(const double& a);


private:

};

template<class TT>
inline BregmanMethod<TT>::BregmanMethod()
{
}

template<class TT>
inline BregmanMethod<TT>::~BregmanMethod()
{
}


template<class TT>
inline Vector2D<TT> BregmanMethod<TT>::Shrink(const Vector2D<TT>& input, const TT & lambda)
{
	return max(input.magnitude() - lambda, 0)*input / input.magnitude();
}

template<class TT>
inline Field2D<Vector2D<TT>> BregmanMethod<TT>::Shrink(const Field2D<Vector2D<TT>>& input, const TT & lambda)
{
	Field2D<Vector2D<TT>> temp = Grid2D(input.grid);
	for (int i = temp.iStart; i <= temp.iEnd; i++)
	{
		for (int j = temp.jStart; j <= temp.jEnd; j++)
		{
			temp(i, j) = Shrink(input(i, j), lambda);
		}
	}
	return temp;
}








