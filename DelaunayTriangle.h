#pragma once
#include "Polygon.h"


class DelaunayTriangulization
{
public:
	DelaunayTriangulization();
	~DelaunayTriangulization();

	static bool CircumCircle(const Vector2D<double> & p, const Vector2D<double> & p1, const Vector2D<double> & p2, const Vector2D<double> & p3, Vector2D<double> & center, double &r);
	static void CircumCircle(const Vector2D<double> & p1, const Vector2D<double> & p2, const Vector2D<double> & p3, Vector2D<double> & center, double &r);
	static void DelaunayTriangulate(const VectorND<Vector2D<double>> & points, VectorND<Polygon2D> & Triangles, VectorND<Vector2D<double>> & center, VectorND<double> & radius);
private:

};

DelaunayTriangulization::DelaunayTriangulization()
{
}

DelaunayTriangulization::~DelaunayTriangulization()
{
}

////////////////////////////////////////////////////////////////////////
// CircumCircle() :
//   Return true if a point (xp,yp) is inside the circumcircle made up
//   of the points (x1,y1), (x2,y2), (x3,y3)
//   The circumcircle centre is returned in (xc,yc) and the radius r
//   Note : A point on the edge is inside the circumcircle
////////////////////////////////////////////////////////////////////////
inline bool DelaunayTriangulization::CircumCircle(const Vector2D<double>& p, const Vector2D<double>& p1, const Vector2D<double>& p2, const Vector2D<double>& p3, Vector2D<double>& center, double & r)
{
	double xp = p.x, yp = p.y;
	double x1 = p1.x, y1 = p1.y;
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;

	double xc, yc;

	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;
	double EPSILON = 0.000000001;
	/* Check for coincident points */
	if (abs(y1 - y2) < EPSILON  && abs(y2 - y3) < EPSILON)
	{
		return(false);
	}

	if (abs(y2 - y1) < EPSILON)
	{
		m2 = -(x3 - x2) / (y3 - y2);
		mx2 = (x2 + x3) / 2.0;
		my2 = (y2 + y3) / 2.0;

		center.x = (x2 + x1) / 2.0;
		center.y = m2 * (center.x - mx2) + my2;
	}
	else if (abs(y3 - y2) < EPSILON)
	{
		m1 = -(x2 - x1) / (y2 - y1);
		mx1 = (x1 + x2) / 2.0;
		my1 = (y1 + y2) / 2.0;

		center.x = (x3 + x2) / 2.0;
		center.y = m1 * (center.x - mx1) + my1;
	}
	else
	{
		m1 = -(x2 - x1) / (y2 - y1);
		m2 = -(x3 - x2) / (y3 - y2);
		mx1 = (x1 + x2) / 2.0;
		mx2 = (x2 + x3) / 2.0;
		my1 = (y1 + y2) / 2.0;
		my2 = (y2 + y3) / 2.0;

		center.x = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		center.y = m1 * (center.x - mx1) + my1;
	}

	dx = x2 - center.x;
	dy = y2 - center.y;
	rsqr = dx * dx + dy * dy;
	r = sqrt(rsqr);
	dx = xp - center.x;
	dy = yp - center.y;
	drsqr = dx * dx + dy * dy;

	if (drsqr <= rsqr)
	{
		return true;
	}
	else
	{
		return false;
	}
}

////////////////////////////////////////////////////////////////////////
// CircumCircle() :
//   Return  circumcircle centre (xc,yc)  and the circumcircle radius r
//   made up of the points (x1,y1), (x2,y2), (x3,y3).
//   Note : A point on the edge is inside the circumcircle
////////////////////////////////////////////////////////////////////////
inline void DelaunayTriangulization::CircumCircle(const Vector2D<double>& p1, const Vector2D<double>& p2, const Vector2D<double>& p3, Vector2D<double>& center, double & r)
{
	double x1 = p1.x, y1 = p1.y;
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;

	double xc, yc;

	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;
	double EPSILON = 0.000001;
	/* Check for coincident points */
	if (abs(y1 - y2) < EPSILON  && abs(y2 - y3) < EPSILON)
	{
		assert(true);
	}

	if (abs(y2 - y1) < EPSILON)
	{
		m2 = -(x3 - x2) / (y3 - y2);
		mx2 = (x2 + x3) / 2.0;
		my2 = (y2 + y3) / 2.0;

		center.x = (x2 + x1) / 2.0;
		center.y = m2 * (center.x - mx2) + my2;
	}
	else if (abs(y3 - y2) < EPSILON)
	{
		m1 = -(x2 - x1) / (y2 - y1);
		mx1 = (x1 + x2) / 2.0;
		my1 = (y1 + y2) / 2.0;

		center.x = (x3 + x2) / 2.0;
		center.y = m1 * (center.x - mx1) + my1;
	}
	else
	{
		m1 = -(x2 - x1) / (y2 - y1);
		m2 = -(x3 - x2) / (y3 - y2);
		mx1 = (x1 + x2) / 2.0;
		mx2 = (x2 + x3) / 2.0;
		my1 = (y1 + y2) / 2.0;
		my2 = (y2 + y3) / 2.0;

		center.x = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		center.y = m1 * (center.x - mx1) + my1;
	}
	dx = x2 - center.x;
	dy = y2 - center.y;
	rsqr = dx * dx + dy * dy;
	r = sqrt(rsqr);
}

inline void DelaunayTriangulization::DelaunayTriangulate(const VectorND<Vector2D<double>>& points, VectorND<Polygon2D> & Triangles, VectorND<Vector2D<double>> & center, VectorND<double> & radius)
{
	int triMax = points.iLength * 4;
	int edgeMax = 100;

	VectorND<Vector2D<int>> Edge(1, edgeMax);
	VectorND<Vector2D<int>> tempEdge(1, 2*edgeMax);

	VectorND<VectorND<int>> triangle(1, triMax);
#pragma omp parallel for
	for (int i = triangle.iStart; i <= triangle.iEnd; i++)
	{
		triangle(i) = VectorND<int>(1, 3);
	}


	VectorND<bool> complete(points.iStart, triMax);


	/*
	Find the maximum and minimum vertex bounds.
	This is to allow calculation of the bounding triangle
	*/
	double xMin = points(points.iStart).x;
	double xMax = points(points.iStart).x;
	double yMin = points(points.iStart).y;
	double yMax = points(points.iStart).y;

	for (int i = points.iStart; i <= points.iEnd; i++)
	{
		if (points(i).x<xMin)
		{
			xMin = points(i).x;
		}
		if (points(i).x>xMax)
		{
			xMax = points(i).x;
		}
		if (points(i).y<yMin)
		{
			yMin = points(i).y;
		}
		if (points(i).y>yMax)
		{
			yMax = points(i).y;
		}
	}

	double dx = xMax - xMin;
	double dy = yMax - yMin;
	double dMax = (dx > dy) ? dx : dy;
	double xMid = (xMax + xMin) / 2.0;
	double yMid = (yMax + yMin) / 2.0;


	/*
	Set up the supertriangle
	his is a triangle which encompasses all the sample points.
	The supertriangle coordinates are added to the end of the
	vertex list. The supertriangle is the first triangle in
	the triangle list.
	*/
	VectorND<Vector2D<double>>  newPoints(points.iStart, points.iLength + 3);
#pragma omp parallel for
	for (int i = points.iStart; i <= points.iEnd; i++)
	{
		newPoints(i) = points(i);
	}
	newPoints(points.iEnd + 1).x = xMid - 20 * dMax;
	newPoints(points.iEnd + 1).y = yMid - dMax;
	newPoints(points.iEnd + 2).x = xMid;
	newPoints(points.iEnd + 2).y = yMid + 20 * dMax;
	newPoints(points.iEnd + 3).x = xMid + 20 * dMax;
	newPoints(points.iEnd + 3).y = yMid - dMax;

	int numTri = 1;
	triangle(numTri)(1) = points.iEnd + 1;
	triangle(numTri)(2) = points.iEnd + 2;
	triangle(numTri)(3) = points.iEnd + 3;
	complete(points.iStart) = false;




	/*
	Include each point one at a time into the existing mesh
	*/
	int nbhdEdge;
	Vector2D<double> currentPoint;
	Vector2D<double> P1, P2, P3;
	center = VectorND<Vector2D<double>>(1, triMax);
	radius = VectorND<double>(1, triMax);
	bool inside;
	for (int i = points.iStart; i <= points.iEnd; i++)
	{
		currentPoint = newPoints(i);
		nbhdEdge = 0;
		//currentPoint.Variable("currentPoint");
		/*
		Set up the edge buffer.
		If the point (xp,yp) lies inside the circumcircle then the
		three edges of that triangle are added to the edge buffer
		and that triangle is removed.
		*/
		for (int j = 1; j <= numTri; j++)
		{
			if (complete(j))
			{
				continue;
			}
			P1 = newPoints(triangle(j)(1));
			P2 = newPoints(triangle(j)(2));
			P3 = newPoints(triangle(j)(3));
			//P1.Variable("P1");
			//P2.Variable("P2");
			//P3.Variable("P3");
			//MATLAB.Command("plot([P1(1) P2(1) P3(1) P1(1)],[P1(2) P2(2) P3(2) P1(2)]),hold on, plot(currentPoint(1),currentPoint(2),'ro'), hold off");

			inside = CircumCircle(currentPoint, P1, P2, P3, center(j), radius(j));

			// Suggested
			// if (xc + r + EPSILON < xp)
			if (center(j).x + radius(j)<currentPoint.x)
			{
				complete(j) = true;
			}

			if (inside)
			{
				/* Check that we haven't exceeded the edge list size */
				if (nbhdEdge + 3 >= edgeMax)
				{
					edgeMax += 100;
					tempEdge = VectorND<Vector2D<int>>(edgeMax);
					#pragma omp parallel for
					for (int i = Edge.iStart; i <= Edge.iEnd; i++)
					{
						tempEdge(i) = Edge(i);
					}
					Edge = tempEdge;
				}
				nbhdEdge++;
				Edge(nbhdEdge) = Vector2D<int>(triangle(j)(1), triangle(j)(2));
				nbhdEdge++;
				Edge(nbhdEdge) = Vector2D<int>(triangle(j)(2), triangle(j)(3));
				nbhdEdge++;
				Edge(nbhdEdge) = Vector2D<int>(triangle(j)(3), triangle(j)(1));

				triangle(j) = triangle(numTri);
				complete(j) = complete(numTri);
				center(j) = center(numTri);
				radius(j) = radius(numTri);

				numTri--;
				j--;
			}
		}

		/*
		Tag multiple edges
		Note: if all triangles are specified anticlockwise then all
		interior edges are opposite pointing in direction.
		*/
#pragma omp parallel for
		for (int j = Edge.iStart; j <= nbhdEdge - 1; j++)
		{
			for (int k = j + 1; k <= nbhdEdge; k++)
			{
				if ((Edge(j).i == Edge(k).j) && (Edge(j).j == Edge(k).i))
				{
					Edge(j) = -1;
					Edge(k) = -1;
				}
				/* Shouldn't need the following, see note above */
				if ((Edge(j).i == Edge(k).i) && (Edge(j).j == Edge(k).j))
				{
					Edge(j) = -1;
					Edge(k) = -1;
				}
			}
		}

		/*
		Form new triangles for the current point
		Skipping over any tagged edges.
		All edges are arranged in clockwise order.
		*/
		for (int j = Edge.iStart; j <= nbhdEdge; j++)
		{
			if (Edge(j).x < 0 || Edge(j).y < 0)
			{
				continue;
			}

			numTri++;
			triangle(numTri)(1) = Edge(j).i;
			triangle(numTri)(2) = Edge(j).j;
			triangle(numTri)(3) = i;
			CircumCircle(newPoints(triangle(numTri)(1)), newPoints(triangle(numTri)(2)), newPoints(triangle(numTri)(3)), center(numTri), radius(numTri));
			complete(numTri) = false;
			//newPoints(triangle(numTri)(1)).Variable("P1");
			//newPoints(triangle(numTri)(2)).Variable("P2");
			//newPoints(triangle(numTri)(3)).Variable("P3");
			//MATLAB.Variable("radius", radius(numTri));
			//MATLAB.Command("plot([P1(1) P2(1) P3(1) P1(1)],[P1(2) P2(2) P3(2) P1(2)])");
			//triangle(numTri).Variable("triangleN");
		}
	}

	/*
	Remove triangles with supertriangle vertices
	These are triangles which have a vertex number greater than nv
	*/
	for (int i = triangle.iStart; i <= numTri; i++)
	{
		if (triangle(i)(1)>points.iEnd || triangle(i)(2)>points.iEnd || triangle(i)(3)>points.iEnd)
		{
			triangle(i) = triangle(numTri);
			center(i) = center(numTri);
			radius(i) = radius(numTri);
			numTri--;
			i--;
		}
	}

	VectorND<Vector2D<double>> tempCenter = center;
	VectorND<double> tempRadius = radius;

	center = VectorND<Vector2D<double>>(1, numTri);
	radius = VectorND<double>(1, numTri);
	Triangles = VectorND<Polygon2D> (1, numTri);

	#pragma omp parallel for
	for (int i = 1; i <= numTri; i++)
	{
		Triangles(i) = Polygon2D(3);
		/*
		triangle is clock-wise
		*/
		Triangles(i).Index(1) = triangle(i).index(1);
		Triangles(i).Index(2) = triangle(i).index(3);
		Triangles(i).Index(3) = triangle(i).index(2);
		Triangles(i)(1) = points(triangle(i)(1));
		Triangles(i)(2) = points(triangle(i)(3));
		Triangles(i)(3) = points(triangle(i)(2));
		center(i) = tempCenter(i);
		radius(i) = tempRadius(i);
	}



	bool show = false;
	if (show)
	{
		string str;
		const char* cmd;

		VecND2DVariable("points", points);
		MATLAB.Command("theta = 0:2*pi/100:2*pi+2*pi/100;");
		MATLAB.Command("plot(points(:,1),points(:,2),'ro'),hold on");
		for (int i = 1; i <= numTri; i++)
		{
			Triangles(i).Plot();
			MATLAB.Command("axis equal");
			str = string("title(['triangle : ', num2str(") + to_string(i) + string("),'/', num2str(") + to_string(numTri) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Variable("i", i);
			center(i).Variable("center");
			radius.Variable("radius");
			MATLAB.Command("plot(radius(i)*cos(theta)+center(1),radius(i)*sin(theta)+center(2)), axis equal");
			str = string("text(sum(tri(:,1))/3,sum(tri(:,2))/3,num2str(i))");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}
	}

}
