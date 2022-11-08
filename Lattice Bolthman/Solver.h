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
	Solver(int Nx, int Ny, vector<vector<vector<double>>>& rho, 
		int numberspec, vector<double> crit_temp, vector<double> crit_rho, 
		vector<double> molmass, vector<double> omega1);
	void SaveVTKFile(int tStep);
	void oneStepWithoutWalls();
	void oneStepWithoutWallsMulticomponent();
	void oneStepWithWalls();
	void check_rho();
	void calculate_r_p1_p2();
	void check_angle();



private:
	double wettability = 1.0;
	double max_rho = 0;
	double min_rho = 10;
	const double A = -0.152;
	const double A2 = -0.0456; /*Пенг-Робинсон*/
	const int tau = 1;
	const double temperature = 0.9;
	int mNx;
	int mNy;
	int number_of_species;
	vector<double> omega;
	vector<double> critical_temperatures;
	vector<double> critical_rho;
	vector<double> molarmass;
	vector<vector<double>> rho;
	vector<vector<double>> effrho;
	vector<vector<double>> sqr_effrho;
	vector<vector<double>> ux;
	vector<vector<double>> uy;
	vector<vector<double>> dux_force;
	vector<vector<double>> duy_force;
	vector<vector<vector<double>>> f;
	vector<vector<vector<double>>> rhomulticomponent;
	vector<vector<vector<vector<double>>>> fmulticomponent;
	vector<vector<double>> pressure;
	const int dx[9] = { 0,   1, 0, -1,  0,   1, -1, -1,  1 };
	const int dy[9] = { 0,   0, 1,  0, -1,   1,  1, -1, -1 };
	
	double PressureVanDerVaals(double& rho,const double& temperature);
	double PressurePengRobinson(double& rho, const double& temperature, double omega);
	void PressurePengRobinsonMulticomponent();
	void set_border_conditions();
	void set_border_conditions_2();
	void set_border_conditions_multicomponent();
	void eq_func(double rho, double ux, double uy, double* f_eq);
	void movement();
	void movement_multicomponent();
	void set_rho_u();
	void set_rho_u_multicomponent();
	void set_effrho();
	void set_effrho_multicomponent();
	void set_effrho_2();
	void set_du_force();
	void collision();
	void collision_multicomponent();

	
};



