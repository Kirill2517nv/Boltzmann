#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

class Solver
{
public:
	Solver(int Nx, int Ny, vector<vector<double>>& rho);
	void SaveVTKFile(int tStep);
	void oneStep();
	void check_rho();
	void calculate_r_p1_p2();

private:
	const double A = -0.152;
	int mNx;
	int mNy;
	const int tau = 1;
	const double temperature = 0.9;
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



