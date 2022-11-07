#include "Solver.h"



Solver::Solver(int Nx, int Ny, vector<vector<double>>& rho) :
    mNx(Nx),
    mNy(Ny),
    rho(rho),
    ux(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    uy(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    dux_force(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    duy_force(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    effrho(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    sqr_effrho(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    f(9, vector<vector<double>>(Nx + 2, vector<double>(Ny + 2)))
{
    double f_eq[9];
    for (int i = 1; i <= mNx; i++)
    {
        for (int j = 1; j <= mNy; j++)
        {
            eq_func(rho[i][j], ux[i][j], uy[i][j], f_eq);
            for (int k = 0; k < 9; k++)
            {
                f[k][i][j] =  f_eq[k];
            }
        }
    }
}

Solver::Solver(int Nx, int Ny, vector<vector<vector<double>>>& rho, 
    int numberspec, vector<double> crit_temp, vector<double> crit_rho, 
    vector<double> molmass) :
    mNx(Nx),
    mNy(Ny),
    rhomulticomponent(rho),
    number_of_species(numberspec),
    critical_temperatures(crit_temp),
    critical_rho(crit_rho),
    molarmass(molmass),
    pressure(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    ux(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    uy(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    dux_force(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    duy_force(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    effrho(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    sqr_effrho(vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))),
    fmulticomponent(numberspec, vector<vector<vector<double>>>(9, vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))))
{
    double f_eq[9];
    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        for (int i = 1; i <= mNx; i++)
        {
            for (int j = 1; j <= mNy; j++)
            {
                eq_func(rhomulticomponent[numspec][i][j], ux[i][j], uy[i][j], f_eq);
                for (int k = 0; k < 9; k++)
                {
                    fmulticomponent[numspec][k][i][j] = f_eq[k];
                }
            }
        }

    }
}

void Solver::SaveVTKFile(int tStep)
{
    stringstream fname;
    fname << "VTK/adv_";
    if (tStep < 10) fname << "0";
    if (tStep < 100) fname << "0";
    if (tStep < 1000) fname << "0";
    if (tStep < 10000) fname << "0";
    if (tStep < 100000) fname << "0";
    if (tStep < 1000000) fname << "0";
    if (tStep < 10000000) fname << "0";
    fname << tStep << ".vtk";
    ofstream vtk_file(fname.str().c_str());
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Immiscible displacement\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << mNx << " " << mNy << " 1\n";
    vtk_file << "X_COORDINATES " << mNx << " double\n";
    for (int i = 1; i <= mNx; i++) vtk_file << i << " ";
    vtk_file << endl;
    vtk_file << "Y_COORDINATES " << mNy << " double\n";
    for (int i = 1; i <= mNy; i++) vtk_file << i << " ";
    vtk_file << endl;
    vtk_file << "Z_COORDINATES 1 double\n0\n";
    vtk_file << "POINT_DATA " << mNx * mNy << endl;
  /*  vtk_file << "SCALARS delta_rho double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++)
            if (mask[i][j] == 0) vtk_file << rho_spec[0][i][j] - rho_spec[1][i][j] << " ";
            else vtk_file << -2.0 << " ";
    vtk_file << endl;*/
    vtk_file << "SCALARS rho double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++) vtk_file << rho[i][j] << " ";
    vtk_file << endl;
    vtk_file << "VECTORS uflow double\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++) vtk_file << ux[i][j] << "  " << uy[i][j] << "  0.0" << " ";
    vtk_file << endl;
    vtk_file.close();

    cout << endl << "File " << fname.str() << " written" << endl << endl;
}


void Solver::oneStepWithWalls()
{
    set_border_conditions_2();
    movement();
    set_rho_u();
    set_effrho_2();
    set_du_force();
    collision();
}

void Solver::oneStepWithoutWalls()
{
    set_border_conditions();
    movement();
    set_rho_u();
    set_effrho();
    set_du_force();
    collision();
}

void Solver::oneStepWithoutWallsMulticomponent()
{
    set_border_conditions_multicomponent();
    movement_multicomponent();
    set_rho_u_multicomponent();
    set_effrho_multicomponent();
    set_du_force();
    collision_multicomponent();
}

double Solver::PressureVanDerVaals(double& rho, const double& temperature) 
{
    double pressure = 8 * temperature * rho / (3 - rho) - 3 * rho * rho;
    return pressure;
}
 /* Peng*/
double a(const double& temperature)
{
    double omega = 0.2514;
    double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    double a = pow((1 + m * (1 - sqrt(temperature))), 2);
    return a;
}

double Solver::PressurePengRobinson(double& rho, const double& temperature)
{
    
    double pressure = 1 / 0.307 * (temperature / (1. / rho - 0.253) - 
        1.487 * a(temperature) / (1. / rho / rho + 2 * 0.253 / rho - 0.253 * 0.253));
    return pressure;
}

void Solver::PressurePengRobinsonMulticomponent()
{
    double Z = 0.307;
    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            double D = 0;
            double B = 0;
            double A = 0;
            for (int numspec = 0; numspec < number_of_species; numspec++)
            {
                D += rhomulticomponent[numspec][i][j] * molarmass[0] / molarmass[numspec];
                B += rhomulticomponent[numspec][i][j] * critical_rho[0] / critical_rho[numspec];
                for (int numspec2 = 0; numspec2 < number_of_species; numspec2++)
                {
                    A += rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec2][i][j] *
                        molarmass[0] / molarmass[numspec] * molarmass[0] / molarmass[numspec2] *
                        sqrt(
                            (critical_temperatures[numspec] / critical_temperatures[0] *
                            molarmass[0] / molarmass[numspec] *
                            critical_rho[0] / critical_rho[numspec]) *
                            (critical_temperatures[numspec2] / critical_temperatures[0] *
                            molarmass[0] / molarmass[numspec2] *
                            critical_rho[0] / critical_rho[numspec2]) *
                            a(temperature * critical_temperatures[0] / critical_temperatures[numspec]) *
                            a(temperature * critical_temperatures[0] / critical_temperatures[numspec2])
                                );
                }

            }
            A *= 1.487;
            B *= 0.253;
            pressure[i][j] = 1. / Z * (D * temperature / (1 - B) - A / (1 + 2 * B - B * B));
        }
    }
    
}

void Solver::set_border_conditions()
{
    f[5][0][0] = f[5][mNx][mNy];
    f[6][mNx + 1][0] = f[6][1][mNy];
    f[7][mNx + 1][mNy + 1] = f[7][1][1];
    f[8][0][mNy + 1] = f[8][mNx][1];

    for (int i = 1; i < mNy + 1; i++) 
    {
        f[1][0][i] = f[1][mNx][i];
        f[3][mNx + 1][i] = f[3][1][i];
        f[5][0][i] = f[5][mNx][i];
        f[8][0][i] = f[8][mNx][i];
        f[6][mNx + 1][i] = f[6][1][i];
        f[7][mNx + 1][i] = f[7][1][i];
    }
    for (int i = 1; i < mNx + 1; i++) 
    {
        f[2][i][0] = f[2][i][mNy];
        f[4][i][mNy + 1] = f[4][i][1];
        f[5][i][0] = f[5][i][mNy];
        f[6][i][0] = f[6][i][mNy];
        f[7][i][mNy + 1] = f[7][i][1];
        f[8][i][mNy + 1] = f[8][i][1];
    }
}

void Solver::set_border_conditions_multicomponent()
{
    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        fmulticomponent[numspec][5][0][0] = fmulticomponent[numspec][5][mNx][mNy];
        fmulticomponent[numspec][6][mNx + 1][0] = fmulticomponent[numspec][6][1][mNy];
        fmulticomponent[numspec][7][mNx + 1][mNy + 1] = fmulticomponent[numspec][7][1][1];
        fmulticomponent[numspec][8][0][mNy + 1] = fmulticomponent[numspec][8][mNx][1];

        for (int i = 1; i < mNy + 1; i++)
        {
            fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][1][mNx][i];
            fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][3][1][i];
            fmulticomponent[numspec][5][0][i] = fmulticomponent[numspec][5][mNx][i];
            fmulticomponent[numspec][8][0][i] = fmulticomponent[numspec][8][mNx][i];
            fmulticomponent[numspec][6][mNx + 1][i] = fmulticomponent[numspec][6][1][i];
            fmulticomponent[numspec][7][mNx + 1][i] = fmulticomponent[numspec][7][1][i];
        }
        for (int i = 1; i < mNx + 1; i++)
        {
            fmulticomponent[numspec][2][i][0] = fmulticomponent[numspec][2][i][mNy];
            fmulticomponent[numspec][4][i][mNy + 1] = fmulticomponent[numspec][4][i][1];
            fmulticomponent[numspec][5][i][0] = fmulticomponent[numspec][5][i][mNy];
            fmulticomponent[numspec][6][i][0] = fmulticomponent[numspec][6][i][mNy];
            fmulticomponent[numspec][7][i][mNy + 1] = fmulticomponent[numspec][7][i][1];
            fmulticomponent[numspec][8][i][mNy + 1] = fmulticomponent[numspec][8][i][1];
        }

    }
}

void Solver::eq_func(double rho, double ux, double uy, double* f_eq)
{
    double du2 = 1.0 - 1.5 * (ux * ux + uy * uy);
    f_eq[0] = 4.0 / 9.0 * rho * du2;
    f_eq[1] = rho / 9.0 * (du2 + ux * (3.0 + 4.5 * ux));
    f_eq[2] = rho / 9.0 * (du2 + uy * (3.0 + 4.5 * uy));
    f_eq[3] = rho / 9.0 * (du2 - ux * (3.0 - 4.5 * ux));
    f_eq[4] = rho / 9.0 * (du2 - uy * (3.0 - 4.5 * uy));
    f_eq[5] = rho / 36.0 * (du2 + (ux + uy) * (3.0 + 4.5 * (ux + uy)));
    f_eq[6] = rho / 36.0 * (du2 + (uy - ux) * (3.0 + 4.5 * (uy - ux)));
    f_eq[7] = rho / 36.0 * (du2 - (ux + uy) * (3.0 - 4.5 * (ux + uy)));
    f_eq[8] = rho / 36.0 * (du2 + (ux - uy) * (3.0 + 4.5 * (ux - uy)));
}

void Solver::movement()
{
    vector<vector<double>> f_temp(mNx + 2, vector<double>(mNy + 2));
    for (int k = 1; k < 9; k++)
    {
        for (int i = 1; i <= mNx; i++)
        {
            for (int j = 1; j <= mNy; j++)
            {
            int i_adj = i - dx[k];
            int j_adj = j - dy[k];
            f_temp[i][j] = f[k][i_adj][j_adj];

            }
        }
        f[k].swap(f_temp);
    }
}

void Solver::movement_multicomponent()
{
    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        vector<vector<double>> f_temp(mNx + 2, vector<double>(mNy + 2));
        for (int k = 1; k < 9; k++)
        {
            for (int i = 1; i <= mNx; i++)
            {
                for (int j = 1; j <= mNy; j++)
                {
                    int i_adj = i - dx[k];
                    int j_adj = j - dy[k];
                    f_temp[i][j] = fmulticomponent[numspec][k][i_adj][j_adj];

                }
            }
            fmulticomponent[numspec][k].swap(f_temp);
        }

    }
}

void Solver::set_rho_u()
{
    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            rho[i][j] = 0;
            for (int k = 0; k < 9; k++) {
                rho[i][j] += f[k][i][j];
            }
            ux[i][j] = (f[1][i][j] + f[5][i][j] + f[8][i][j] - f[3][i][j] - f[6][i][j] - f[7][i][j]) / rho[i][j];
            uy[i][j] = (f[2][i][j] + f[5][i][j] + f[6][i][j] - f[4][i][j] - f[7][i][j] - f[8][i][j]) / rho[i][j];
            
        }
    }
}

void Solver::set_rho_u_multicomponent()
{
    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            rho[i][j] = 0;
            for (int numspec = 0; numspec < number_of_species; numspec++)
            {
                rhomulticomponent[numspec][i][j] = fmulticomponent[numspec][0][i][j];
                for (int k = 1; k < 9; k++) 
                {
                    rhomulticomponent[numspec][i][j] += fmulticomponent[numspec][k][i][j];
                }
                rho[i][j] += rhomulticomponent[numspec][i][j];
                ux[i][j] = (fmulticomponent[numspec][1][i][j] + fmulticomponent[numspec][5][i][j] +
                    fmulticomponent[numspec][8][i][j] - fmulticomponent[numspec][3][i][j] -
                    fmulticomponent[numspec][6][i][j] - fmulticomponent[numspec][7][i][j]);
                uy[i][j] = (fmulticomponent[numspec][2][i][j] + fmulticomponent[numspec][5][i][j] + 
                    fmulticomponent[numspec][6][i][j] - fmulticomponent[numspec][4][i][j] - 
                    fmulticomponent[numspec][7][i][j] - fmulticomponent[numspec][8][i][j]);
                
            }
            ux[i][j] = ux[i][j] / rho[i][j];
            uy[i][j] = uy[i][j] / rho[i][j];
            
        }
    }    
    
}


void Solver::set_effrho()
{
    
    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            sqr_effrho[i][j] = rho[i][j] / 3.0 - 0.01 * PressurePengRobinson(rho[i][j], temperature);
            effrho[i][j] = sqrt(sqr_effrho[i][j]);
        }
    }
    effrho[0][mNy + 1] = effrho[mNx][1];
    effrho[mNx + 1][mNy + 1] = effrho[1][1];
    effrho[0][0] = effrho[mNx][mNy];
    effrho[mNx + 1][0] = effrho[1][mNy];
    sqr_effrho[0][0] = effrho[0][0] * effrho[0][0];
    sqr_effrho[0][mNy + 1] = effrho[0][mNy + 1] * effrho[0][mNy + 1];
    sqr_effrho[mNx + 1][0] = effrho[mNx + 1][0] * effrho[mNx + 1][0];
    sqr_effrho[mNx + 1][mNy + 1] = effrho[mNx + 1][mNy + 1] * effrho[mNx + 1][mNy + 1];
    for (int i = 1; i <= mNx; i++) {
        effrho[i][mNy + 1] = effrho[i][1];
        effrho[i][0] = effrho[i][mNy];
        sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
        sqr_effrho[i][0] = effrho[i][0] * effrho[i][0];
    }
    for (int i = 1; i <= mNy; i++) {
        effrho[0][i] = effrho[mNx][i];
        effrho[mNx + 1][i] = effrho[1][i];
        sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
        sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
    }
}

void Solver::set_effrho_multicomponent()
{

    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            sqr_effrho[i][j] = rho[i][j] / 3.0 - 0.01 * pressure[i][j];
            effrho[i][j] = sqrt(sqr_effrho[i][j]);
        }
    }
    effrho[0][mNy + 1] = effrho[mNx][1];
    effrho[mNx + 1][mNy + 1] = effrho[1][1];
    effrho[0][0] = effrho[mNx][mNy];
    effrho[mNx + 1][0] = effrho[1][mNy];
    sqr_effrho[0][0] = effrho[0][0] * effrho[0][0];
    sqr_effrho[0][mNy + 1] = effrho[0][mNy + 1] * effrho[0][mNy + 1];
    sqr_effrho[mNx + 1][0] = effrho[mNx + 1][0] * effrho[mNx + 1][0];
    sqr_effrho[mNx + 1][mNy + 1] = effrho[mNx + 1][mNy + 1] * effrho[mNx + 1][mNy + 1];
    for (int i = 1; i <= mNx; i++) {
        effrho[i][mNy + 1] = effrho[i][1];
        effrho[i][0] = effrho[i][mNy];
        sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
        sqr_effrho[i][0] = effrho[i][0] * effrho[i][0];
    }
    for (int i = 1; i <= mNy; i++) {
        effrho[0][i] = effrho[mNx][i];
        effrho[mNx + 1][i] = effrho[1][i];
        sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
        sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
    }
}

void Solver::set_du_force()
{
    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            double force_x = 0;
            force_x = 2.0 / 3 * ((1 - 2 * A2) * effrho[i][j] * (effrho[i + 1][j] - effrho[i - 1][j] + 
            0.25 * (effrho[i + 1][j + 1] + effrho[i + 1][j - 1] - effrho[i - 1][j + 1] - effrho[i - 1][j - 1])) +
            A2 * (sqr_effrho[i + 1][j] - sqr_effrho[i - 1][j] + 0.25 *
            (sqr_effrho[i + 1][j + 1] + sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j + 1] - sqr_effrho[i - 1][j - 1])));
            dux_force[i][j] = force_x / rho[i][j];
            double force_y = 0;
            force_y = 2.0 / 3 * ((1 - 2 * A2) * effrho[i][j] * (effrho[i][j + 1] - effrho[i][j - 1] +
                0.25 * (effrho[i + 1][j + 1] + effrho[i - 1][j + 1] - effrho[i + 1][j - 1] - effrho[i - 1][j - 1])) +
                A2 * (sqr_effrho[i][j + 1] - sqr_effrho[i][j - 1] + 0.25 *
                    (sqr_effrho[i + 1][j + 1] + sqr_effrho[i - 1][j + 1] - sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j - 1])));
            duy_force[i][j] = force_y / rho[i][j];
        } 
    }
}


void Solver::collision()
{
    double f_eq[9], f_eq_n[9];
    for (int i = 1; i <= mNx; i++)
    {
        for (int j = 1; j <= mNy; j++)
        {
            eq_func(rho[i][j], ux[i][j], uy[i][j], f_eq);
            eq_func(rho[i][j], ux[i][j] + dux_force[i][j], uy[i][j] + duy_force[i][j], f_eq_n);
            for (int k = 0; k < 9; k++)
            {
                f[k][i][j] += f_eq_n[k] - f_eq[k] + (f_eq[k] - f[k][i][j]) / tau;
                /*f[k][i][j] += (f_eq[k] - f[k][i][j]) / tau;   only collision*/ 
            }
        }
    }

}

void Solver::collision_multicomponent()
{
    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        double f_eq[9], f_eq_n[9];
        for (int i = 1; i <= mNx; i++)
        {
            for (int j = 1; j <= mNy; j++)
            {
                eq_func(rhomulticomponent[numspec][i][j], ux[i][j], uy[i][j], f_eq);
                eq_func(rhomulticomponent[numspec][i][j], ux[i][j] + dux_force[i][j], uy[i][j] + duy_force[i][j], f_eq_n);
                for (int k = 0; k < 9; k++)
                {
                    fmulticomponent[numspec][k][i][j] += f_eq_n[k] - f_eq[k] + (f_eq[k] - fmulticomponent[numspec][k][i][j]) / tau;
                    /*f[k][i][j] += (f_eq[k] - f[k][i][j]) / tau;   only collision*/
                }
            }
        }
    }

}

void Solver::check_rho()
{
    double full_rho = 0;
    double max_rho_temp = 0;
    double min_rho_temp = 10;
    int max_rho_x = 0;
    int max_rho_y = 0;
    int min_rho_x = 0;
    int min_rho_y = 0;
    for (int i = 1; i <= mNx; i++)
    {
        for (int j = 1; j <= mNy; j++)
        {
            full_rho += rho[i][j];
            if (rho[i][j] < min_rho_temp)
            {
                min_rho = rho[i][j];
                min_rho_temp = rho[i][j];
            }
            if (rho[i][j] > max_rho_temp)
            {
                max_rho = rho[i][j];
                max_rho_temp = rho[i][j];
            }
            if (max_rho_temp == rho[i][j])
            {
                max_rho_x = i;
                max_rho_y = j;
            }
            if (min_rho_temp == rho[i][j])
            {
                min_rho_x = i;
                min_rho_y = j;
            }
        }
    }
    cout << "full rho = " << setprecision(6) << full_rho << endl;
    cout << "max rho = " << setprecision(6) << max_rho << "\t" << "min rho = " << min_rho << endl;
    cout << "max rho(x) = " << setprecision(6) << max_rho_x << "\t" << "max rho(y) = " << max_rho_y << endl;
    cout << "min rho(x) = " << setprecision(6) << min_rho_x << "\t" << "min rho(y) = " << min_rho_y << endl;
    cout << "rho[100][100] = " << setprecision(6) << rho[100][100] << "\t" << "rho[5][5] = " << rho[5][5] << endl;
}


void Solver::calculate_r_p1_p2()
{
    
    double r = 0;
    double s = 0;
    double p1 = 0;
    double p2 = 0;

    
    for (int i = 1; i <= mNx; i++)
    {
        for (int j = 1; j <= mNy; j++)
        {
            if (rho[i][j] > (max_rho + min_rho) / 2.)
            {
                s++;
            }
        }
    }
    r = sqrt(s / M_PI);
    /*p1 = PressurePengRobinson(rho[100][100], temperature);
    p2 = PressurePengRobinson(rho[5][5], temperature);*/

    p1 = PressurePengRobinson(rho[100][100], temperature);
    p2 = PressurePengRobinson(rho[5][5], temperature);

    cout << "p1 = " << setprecision(6) << p1 << "\tp2 = " << p2 << "\tdelta p = " << p1 - p2 << endl;
    stringstream fname;
    fname << "DATA/T = 0.9 r = 50.txt";
    ofstream vtk_file(fname.str().c_str(), ios_base::app);
    vtk_file << "Radius = " << setprecision(6) << r << "\t";
    vtk_file << "P in drop = " << setprecision(6) << p1 << "\t";
    vtk_file << "P external = " << setprecision(6) << p2 << "\t";
    vtk_file << "Delta P = " << setprecision(6) << p1 - p2 << "\n";
    vtk_file.close();

    cout << endl << "File " << fname.str() << " written" << endl << endl;
}

void Solver::set_border_conditions_2()
{

    for (int i = 1; i < mNy + 1; i++)
    {
        f[3][mNx + 1][i] = f[1][mNx][i];
        f[1][0][i] = f[3][1][i];
        f[7][mNx + 1][i + 1] = f[5][mNx][i];
        f[6][mNx + 1][i - 1] = f[8][mNx][i];
        f[8][0][i + 1] = f[6][1][i];
        f[5][0][i - 1] = f[7][1][i];
    }
    for (int i = 1; i < mNx + 1; i++)
    {
        f[4][i][mNy + 1] = f[2][i][mNy];
        f[2][i][0] = f[4][i][1];
        f[7][i + 1][mNy + 1] = f[5][i][mNy];
        f[8][i - 1][mNy + 1] = f[6][i][mNy];
        f[5][i - 1][0] = f[7][i][1];
        f[6][i + 1][0] = f[8][i][1];
    }
}

void Solver::set_effrho_2()
{

    for (int i = 1; i <= mNx; i++) {
        for (int j = 1; j <= mNy; j++) {
            sqr_effrho[i][j] = rho[i][j] / 3.0 - 0.01 * PressurePengRobinson(rho[i][j], temperature);
            effrho[i][j] = sqrt(sqr_effrho[i][j]);
        }
    }
    effrho[0][mNy + 1] = effrho[1][mNy];
    effrho[mNx + 1][mNy + 1] = effrho[mNx][mNy];
    effrho[0][0] = wettability * effrho[1][1];
    effrho[mNx + 1][0] = wettability * effrho[mNx][1];
    sqr_effrho[0][0] = wettability * wettability * effrho[0][0] * effrho[0][0];
    sqr_effrho[0][mNy + 1] = effrho[0][mNy + 1] * effrho[0][mNy + 1];
    sqr_effrho[mNx + 1][0] = wettability * wettability * effrho[mNx + 1][0] * effrho[mNx + 1][0];
    sqr_effrho[mNx + 1][mNy + 1] = effrho[mNx + 1][mNy + 1] * effrho[mNx + 1][mNy + 1];
    for (int i = 1; i <= mNx; i++) {
        effrho[i][mNy + 1] = effrho[i][mNy];
        effrho[i][0] = wettability * effrho[i][1];
        sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
        sqr_effrho[i][0] = wettability * wettability * effrho[i][0] * effrho[i][0];
    }
    for (int i = 1; i <= mNy; i++) {
        effrho[0][i] = effrho[1][i];
        effrho[mNx + 1][i] = effrho[mNx][i];
        sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
        sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
    }
}

void Solver::check_angle()
{
    double h = 0;
    double l = 0;
    for (int i = 200; i <= mNx; i++)
    {
        if (rho[i][1] > (max_rho + min_rho) / 2) {
            l++;
        }
    }
    l = l + (rho[200 + l][1] - (max_rho + min_rho) / 2) / (rho[200 + l][1] - rho[200 + l + 1][1]);
    for (int i = 1; i <= mNy; i++)
    {
        if (rho[200][i] > (max_rho + min_rho) / 2) {
            h++;
        }
    }
    h = h + (rho[200][h] - (max_rho + min_rho) / 2) / (rho[200][h] - rho[200][h + 1]);
    cout << "h = " << h << "\tl = " << l << "\talpha = " << 2*atan(h / l) * 180.0 / M_PI << endl;
    stringstream fname;
    fname << "DATA/T = 0.7 wett = 1.3.txt";
    ofstream vtk_file(fname.str().c_str(), ios_base::app);
    vtk_file << "h = " << h << "\t";
    vtk_file << "l = " << l << "\t";
    vtk_file << "alpha = " << 2 * atan(h / l) * 180.0 / M_PI << "\n";
    vtk_file.close();

    cout << endl << "File " << fname.str() << " written" << endl << endl;
    
}
