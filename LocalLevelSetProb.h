#pragma once

#include "AdvectionMethod2D.h"

class LocalLevelSetAdvection
{
public:
	Grid2D grid;
	LS LLS;
	double cflCondition;
	double dt;

	bool isVelocity;
	bool needReinitial;
	FD velocityX;
	FD velocityY;

	FD quantity;

	int reinitialIter;
	int maxIteration;
	int writeIter;

	// Zero level set point.
	int givenPointNum;
	VectorND<VT> givenPoint;
	FD distance;

	LocalLevelSetAdvection();
	~LocalLevelSetAdvection();

	void InitialCondition(const int & example);
	void AdvectionSolver(const int & example);
	
	void QuantityExtensionSolver(const int & example);

	double Distance2Data(const int& i, const int& j);
	void ExactDistance();
private:

};

LocalLevelSetAdvection::LocalLevelSetAdvection()
{
}

LocalLevelSetAdvection::~LocalLevelSetAdvection()
{
}

inline void LocalLevelSetAdvection::InitialCondition(const int & example)
{
	if (example == 1)
	{
		isVelocity = true;
		needReinitial = false;

		cout << "Local Local Level set advection Test - Rotating a circle" << endl;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);

		LLS = LS(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}
		LLS.InitialTube();

		velocityX = FD(grid);
		velocityY = FD(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				velocityX(i, j) = -grid(i, j).y;
				velocityY(i, j) = grid(i, j).x;
			}
		}

		isVelocity = true;
		cflCondition = 0.5;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 100;
	}
	else if (example == 2)
	{
		cout << "Local Level set advection Test - Unit circle" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		LLS = LS(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - 0.5;
			}
		}
		LLS.InitialTube();
		cflCondition = 0.5;

		//dt = grid.dx*grid.dy;
		maxIteration = 100;
		writeIter = 10;
	}
	else if (example == 3)
	{
		cout << "Quantity Extension using Local LS - Unit circle" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		LLS = LS(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - 0.5;
			}
		}
		LLS.InitialTube();

		quantity = FD(grid);

		int i, j;
		double x0, y0;
#pragma omp parallel for private(i, j, x0, y0)
		for (int k = 1; k <= LLS.numTube; k++)
		{
			i = LLS.tubeIndex(k).i;
			j = LLS.tubeIndex(k).j;
			x0 = grid(i, j).x;
			y0 = grid(i, j).y;
			if (LLS.tube(i,j) ==1)
			{
				quantity(i, j) = 1; sin(y0 / sqrt(x0*x0 + y0*y0)); x0; x0*x0;
			}
		}

		cflCondition = 0.2;
		dt = cflCondition * min(grid.dx, grid.dy);
		//dt = grid.dx*grid.dy;
		maxIteration = (int)ceil((LLS.gamma2 - LLS.gamma1) / dt);
		writeIter = 10;
	}
}


inline void LocalLevelSetAdvection::AdvectionSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;

	InitialCondition(example);

	grid.Variable();
	LLS.phi.Variable("phi0");
	LLS.phi.Variable("phi");
	LLS.tube.Variable("Tube");
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("subplot(1,3,1)");
	MATLAB.Command("surf(X, Y, Tube);grid on;axis([-1 1 -1 1]);axis equal;");
	MATLAB.Command("subplot(1,3,2)");
	//MATLAB.Command("contour(X, Y, Tube);hold on;grid on;axis([-1 1 -1 1]);axis equal;");
	MATLAB.Command("surf(X,Y,phi);hold off;");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	MATLAB.Command("subplot(1,3,3)");
	MATLAB.Command("contour(X, Y, phi0, [0 0],'b'),grid on;axis([-1 1 -1 1]);axis equal;");

	if (writeFile)
	{
		str = "phi0";
		LLS.phi.WriteFile(str);

		if (isVelocity)
		{
			str = "velocityX";
			velocityX.WriteFile(str);
			str = "velocityY";
			velocityY.WriteFile(str);
		}
	}

	//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	//MATLAB.Command("v = VideoWriter('newfile.avi');");
	//MATLAB.Command("open(v)");
	clock_t before;
	double  result;
	int reinitialIter;
	double totalT = 0;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << endl;
		cout << "********************************" << endl;
		cout << "Level set advection : " << i << endl;

		if (isVelocity)
		{
			dt = AdvectionMethod2D<double>::AdaptiveTimeStep(velocityX, velocityY, cflCondition);
			totalT += dt;
			before = clock();
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, velocityX, velocityY, dt);
			reinitialIter = int(LLS.gamma1 / min(LLS.phi.dx, LLS.phi.dy));
			AdvectionMethod2D<double>::LLSReinitializationTVDRK3(LLS, dt, reinitialIter);
			LLS.UpdateInterface();
			LLS.UpdateLLS();

			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << result << endl;
			LLS.tube.Variable("Tube");
			LLS.phi.Variable("phi");
			MATLAB.Command("subplot(1,3,1)");
			MATLAB.Command("surf(X, Y, Tube);grid on;axis([-1 1 -1 1]);axis equal;set(gca,'fontsize',20)");
			MATLAB.Command("subplot(1,3,2)");
			//MATLAB.Command("contour(X, Y, Tube);hold on;grid on;axis([-1 1 -1 1]);axis equal;");
			MATLAB.Command("surf(X,Y,phi);hold off;set(gca,'fontsize',20)");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("subplot(1,3,3)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;contour(X, Y, phi,[0 0],'r');axis([-1 1 -1 1]);axis equal;hold off;set(gca,'fontsize',20)");

		}
		else
		{
			dt = cflCondition*min(grid.dx, grid.dy);
			before = clock();
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, dt);
			reinitialIter = int(LLS.gamma1 / min(LLS.phi.dx, LLS.phi.dy));
			AdvectionMethod2D<double>::LLSReinitializationTVDRK3(LLS, dt, reinitialIter);
			LLS.UpdateInterface();
			LLS.UpdateLLS();

			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << result << endl;
			LLS.phi.Variable("phi");
			MATLAB.Command("subplot(1,2,1)");
			MATLAB.Command("surf(X,Y,phi);");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("subplot(1,2,2)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;contour(X, Y, phi,[0 0],'r');axis([-1 1 -1 1]);axis equal;hold off;");
			//MATLAB.Command("F=getframe;");
			//MATLAB.Command("writeVideo(v,F)");
		}


		if (needReinitial)
		{
			for (int j = 0; j < reinitialIter; j++)
			{
				cout << "Reinitialization : " << i << "-" << j + 1 << endl;
				dt = cflCondition*min(grid.dx, grid.dy);
				AdvectionMethod2D<double>::LLSReinitializationTVDRK3(LLS, dt);
			}
		}

		if (i%writeIter == 0 && writeFile)
		{
			str = "phi" + to_string(i);
			LLS.phi.WriteFile(str);
		}
	}
	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << result << endl;
	//MATLAB.Command("close(v)");

}

inline void LocalLevelSetAdvection::QuantityExtensionSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;
	int timeAdvectionOrder = 3; // 1 or 3
	int spacialOrder = 3; // 3 or 5
	InitialCondition(example);

	grid.Variable();
	LLS.phi.Variable("phi");
	LLS.tube.Variable("Tube");
	quantity.Variable("quantity0");
	
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("surf(X,Y,quantity0);hold on;contour(X, Y, Tube);hold off,");
	clock_t before = clock();
	double  result;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << endl;
		cout << "********************************" << endl;
		cout << "Quantity Extension : " << i << endl;
		AdvectionMethod2D<double>::LLSQuantityExtension(LLS, quantity, timeAdvectionOrder, spacialOrder);
		quantity.Variable("quantity");
		//MATLAB.Command("subplot(1,2,1)");
		MATLAB.Command("surf(X,Y,quantity);%hold on;%contour(X, Y, Tube);");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		//MATLAB.Command("subplot(1,2,2)");
		//MATLAB.Command("surf(X,Y,quantity-quantity0);%hold on;contour(X, Y, Tube);");
	}
	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << result << endl;
}

inline double LocalLevelSetAdvection::Distance2Data(const int & i, const int & j)
{
	double distance = 100;
	double tempDist;

	for (int k = 0; k < givenPointNum; k++)
	{
		tempDist = (givenPoint(k) - grid(i, j)).magnitude();
		if (distance > tempDist)
		{
			distance = tempDist;
		}
	}

	return distance;
}

inline void LocalLevelSetAdvection::ExactDistance()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			distance(i, j) = Distance2Data(i, j);
		}
	}
}


