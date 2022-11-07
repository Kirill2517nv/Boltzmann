#include "Solver.h"

int main()
{
    system("mkdir VTK");
    system("mkdir DATA");
    int nx = 200;
    int ny = 200;
    int max_time = 40000;
    int numberspec = 2;
    vector<double> critical_temperatures = {3.435, 3.435}; /*[k]*/
    vector<double> critical_rho = { 469,65, 469.65 }; /*[kg/m^3]*/
    vector<double> molar_mass = { 0.07215, 0.07215 }; /*[kg/mol]*/
    vector<vector<vector<double>>> rho (numberspec, vector<vector<double>>(nx + 2, vector<double>(ny + 2, 0.96))); 
    /*for (int i = 1; i < nx + 1; i++)
    {
        for (int j = 1; j < ny + 1; j++)
        {
            if ((i - 100) * (i - 100)  + (j - 100) * (j - 100) <= 2500)
            {
                rho[i][j] = 2.201;
            }
        }
    }*/
    for (int numspec = 0; numspec < numberspec; numspec++)
    {
        for (int i = 1; i <= nx; i++)
        {
            for (int j = 1; j <= ny; j++)
            {
                rho[numspec][i][j] += (rand() % 10) / 100.0;
            }
        }
    }

    Solver A(nx, ny, rho, numberspec, critical_temperatures, critical_rho, molar_mass);
    A.SaveVTKFile(0);
    for (int t = 1; t <= max_time; t++) 
    {
        
        A.oneStepWithoutWalls();
        if (t % 100  == 0) 
        {
            A.SaveVTKFile(t);
            A.check_rho();
            /*A.calculate_r_p1_p2();*/
            /*A.check_angle();*/
        }
        /*if (t % 1000 == 0)
        {
            A.calculate_r_p1_p2();
        }*/
    }
}


