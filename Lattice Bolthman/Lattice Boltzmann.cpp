#include "Solver.h"

int main()
{
    system("mkdir VTK");
    system("mkdir DATA");
    int nx = 200;
    int ny = 200;
    int max_time = 40000;
    int numberspec = 2;
    int x_centre = 100;
    int y_centre = 100;
    int radius = 30;
    int thickness = 8;
    double rho_l = 1.0;
    double rho_v = 0.2;
    vector<double> omega = { 0.2514, 0.1522 };
    vector<double> critical_temperatures = { 469.65, 190.6 }; /*[k] пентан(omega = 0.2514), метан (omega = 0.1522) */
    vector<double> critical_rho = { 232, 162 }; /*[kg/m^3]*/
    vector<double> molar_mass = { 0.07215, 0.016 }; /*[kg/mol]*/
    vector<double> gamma = {1, 1};
    vector<vector<vector<double>>> rho (numberspec, vector<vector<double>>(nx + 2, vector<double>(ny + 2, 0.46))); 
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

    /*for (int i = 1; i < nx + 1; i++)
    {
        for (int j = 1; j < ny + 1; j++)
        {
            rho[0][i][j] = rho_l / 2 - rho_l / 2 * 
                tanh(2 * (sqrt((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) - radius) / thickness);
            rho[1][i][j] = rho_v / 2 + rho_v / 2 *
                tanh(2 * (sqrt((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) - radius) / thickness);
        }
    }*/


    Solver A(nx, ny, rho, numberspec, critical_temperatures, critical_rho, molar_mass, omega, gamma);
    A.SaveVTKFile(0);
    for (int t = 1; t <= max_time; t++) 
    {
        
        A.oneStepWithoutWallsMulticomponent();
        if (t % 50  == 0) 
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


