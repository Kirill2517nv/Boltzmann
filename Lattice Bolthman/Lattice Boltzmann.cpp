#include "Solver.h"

int main()
{
    system("mkdir VTK");
    system("mkdir DATA");
    int nx = 200;
    int ny = 200;
    int max_time = 20000;
    vector<vector<double>> rho(nx + 2, vector<double>(ny + 2, 0.00825)); 
    for (int i = 1; i < nx + 1; i++)
    {
        for (int j = 1; j < ny + 1; j++)
        {
            if ((i - 100) * (i - 100)  + (j - 100) * (j - 100) <= 1225)
            {
                rho[i][j] = 3.253;
            }
        }
    }
    /*for (int i = 1; i <= nx ; i++)
    {
        for (int j = 1; j <= ny; j++)
        {
            rho[i][j] += (rand() % 10) / 100.0;
        }
    }*/

    Solver A(nx, ny, rho);
    A.SaveVTKFile(0);
    for (int t = 1; t <= max_time; t++) 
    {
        
        A.oneStepWithWalls();
        if (t % 100  == 0) 
        {
            A.SaveVTKFile(t);
            /*A.check_angle();*/
            A.check_rho();
            A.calculate_r_p1_p2();
        }
        /*if (t % 1000 == 0)
        {
            A.calculate_r_p1_p2();
        }*/
    }
}


