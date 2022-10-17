#include "Solver.h"

int main()
{
    system("mkdir VTK");
    system("mkdir DATA");
    int nx = 400;
    int ny = 100;
    int max_time = 10000;
    double ful_rho = 0;
    vector<vector<double>> rho(nx + 2, vector<double>(ny + 2, 0.42)); /*при t=0.8 rho_min = 0.24, rho_max = 1.93 при t = 0.9 rho_min = 0.42 rho_max = 1.68*/
    for (int i = 1; i < nx + 1; i++)
    {
        for (int j = 1; j < ny + 1; j++)
        {
            if ((i - 200) * (i - 200)  + (j) * (j) <= 1600)
            {
                rho[i][j] = 1.68;
            }
            ful_rho += rho[i][j];
        }
    }
    cout << ful_rho << endl;
   /* for (int i = 1; i <= nx ; i++)
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
            A.check_angle();
            A.calculate_r_p1_p2();
        }
        /*if (t % 1000 == 0)
        {
            A.calculate_r_p1_p2();
        }*/
    }



}


