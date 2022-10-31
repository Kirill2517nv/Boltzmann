#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <numbers>

using namespace std;

class Solver
{
public:
	Solver(int Nx, int Ny, vector<vector<double>>& rho);
	void SaveVTKFile(int tStep);
	void oneStepWithoutWalls();
	void oneStepWithWalls();
	void check_rho();
	void calculate_r_p1_p2();
	void check_angle();



private:
	double wettability = 1.0;
	double max_rho = 0;
	double min_rho = 2;
	const double A = -0.152;
	int mNx;
	int mNy;
	const int tau = 1;
	const double temperature = 0.6;
	vector<vector<double>> rho;
	vector<vector<double>> effrho;
	vector<vector<double>> sqr_effrho;
	vector<vector<double>> ux;
	vector<vector<double>> uy;
	vector<vector<double>> dux_force;
	vector<vector<double>> duy_force;
	vector<vector<vector<double>>> f;
	const int dx[9] = { 0,   1, 0, -1,  0,   1, -1, -1,  1 };
	const int dy[9] = { 0,   0, 1,  0, -1,   1,  1, -1, -1 };
	
	double PressureVanDerVaals(double& rho,const double& temperature);
	double PressurePengRobinson(double& rho, const double& temperature);
	void set_border_conditions();
	void set_border_conditions_2();
	void eq_func(double rho, double ux, double uy, double* f_eq);
	void movement();
	void set_rho_u();
	void set_effrho();
	void set_effrho_2();
	void set_du_force();
	void collision();

	
};



