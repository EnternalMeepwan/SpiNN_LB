#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define Q 9
// define the problem size
#define xdim 512
#define ydim 512

// Define velocity vectors
const int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

// fi
double fi[xdim + 2][ydim + 2][Q];
double feq[xdim + 2][ydim + 2][Q];
double fi_star[xdim + 2][ydim + 2][Q] = {0.0};

// Wi
const double wi[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

// macro fluid density rho
double rho[xdim + 2][ydim + 2];

double u_x[xdim + 2][ydim + 2] = {0.001};
double u_y[xdim + 2][ydim + 2] = {1.0};

// parameters for check
double mass = 0.0;

// parameter for the checking problem
double U_0 = 0.01;
double K = 30.0;
double delta = 0.05;
double cVel = 1 / sqrt(3); // c_s

/*
 * Init the velocity u, v = (u_x, u_y)
 * 
 * ydim = xdim = N = 128
 * K = 30 / N
 * delta = 0.05
 * 
 * u = tanh( K (y - 0.25 * ydim) ) for y <= 0.5 * ydim 
 * u = U_0 tanh( K (0.75 * ydim - y) ) for y >  0.5 * ydim
 * 
 * v = delta sin( 2pi *x / xdim )
*/
void initRho_N_U(double (*u_x)[ydim + 2], double (*u_y)[ydim + 2], double (*rho)[ydim + 2])
{
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            double x_temp, y_temp;
            x_temp = 1.0 * (double)(i - 1) / (double)xdim;
            y_temp = 1.0 * (double)(j - 1) / (double)ydim; // scale the i and j to make 0 < i,j < 1

            if (y_temp <= 0.5)
            {
                u_x[i][j] = U_0 * tanh(K * (y_temp - 0.25));
            }
            else
            {
                u_x[i][j] = U_0 * tanh(K * (0.75 - y_temp));
            }

            u_y[i][j] = U_0 * delta * sin(2.0 * M_PI * (x_temp + 0.25));
            rho[i][j] = 1.0;
        }
    }
}

double findBig(double (*u)[ydim + 2])
{
    double big = 0.0;
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            if (big < u[i][j])
                big = u[i][j];
        }
    }
    return big;
}
double findSmall(double (*u)[ydim + 2])
{
    double small = 0.0;
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            if (small > u[i][j])
                small = u[i][j];
        }
    }
    return small;
}

void printux()
{
    printf("the init ux = ");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("%f ", u_x[i][j]);
        }
        printf("\n");
    }
}

int getDelta(int a, int b)
{
    if (a == b)
        return 1;
    return 0;
}

/*
 * init the distribution function fi 
*/
void initFi(double (*fi)[ydim + 2][Q], double (*rho)[ydim + 2])
{

    // for check problem
    double A = 1 / (cVel * cVel);
    double B = 1 / 2 * cVel * cVel * cVel * cVel;
    for (int i = 0; i <= xdim + 1; i++)
    {
        for (int j = 0; j <= ydim + 1; j++)
        {
            for (int k = 0; k < Q; k++)
            {
                fi[i][j][k] = wi[k] * rho[i][j] *
                              (1.0 + A * (ex[k] * u_x[i][j] + ey[k] * u_y[i][j]) //  A c_ia u_a =  A (c_ix u_x + c_iy u_y)
                               + B * ((ex[k] * ex[k] - 1.0 / 3.0)                // Q_ixx = (c_ix c_ix - 1/3 c_s^2)
                                      + (ex[k] * ey[k])                          // Q_ixy = (c_ix c_iy)
                                      + (ey[k] * ex[k])                          // Q_iyx = (c_iy c_ix)
                                      + (ey[k] * ey[k] - 1.0 / 3.0))             //Q_iyy = (c_iy c_iy - 1/3 c_s^2)
                              );
            }
        }
    }
}

void printfi(int index)
{
    printf("the fi is ");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("%f ", fi[i][j][index]);
        }
        printf("\n");
    }
}

void printfeq()
{
    printf("the feq is ");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("%f ", feq[i][j][4]);
        }
        printf("\n");
    }
}

/*
 *  Assuming zero velocity as init for every lattice
 *  Use computeFeq() instead !!
 */

void initFeq(double (*rho)[ydim + 2], double (*feq)[ydim + 2][Q])
{
    double c_s = 1 / sqrt(3); // the speed of sound in the lattice
    for (int i = 0; i <= xdim + 1; i++)
    {
        for (int j = 0; j <= ydim + 1; j++)
        {
            for (int k = 0; k < Q; k++)
            {
                //feq_i=ω_i*ρ*(1+(e_i•u_i)/(cs^2)+((e_i•u_i)^2)/(2*(cs^4))+(u_i^2)/(2*(cs^2))
                // feq[i][j][k] = wi[k] * rho[i][j] * ( 1 + (ex[k]*u_x[i][j]/(cVel*cVel)) + (()) );
            }
        }
    }
}

// move fi --> fi_star in the direction of ei
void streamStep(double (*fi_star)[ydim + 2][Q], double (*fi)[ydim + 2][Q])
{

    // periodic condition

    for (int j = 1; j < ydim + 1; j++)
    {
        for (int k = 0; k < Q; k++)
        {
            // copy the extra line in the top to the bottom
            fi_star[0][j][k] = fi_star[xdim][j][k];
            fi_star[xdim + 1][j][k] = fi_star[1][j][k];
        }
    }

    for (int i = 0; i < xdim + 2; i++)
    {
        for (int k = 0; k < Q; k++)
        {
            // copy the extra line in the right to the beginning
            fi_star[i][0][k] = fi_star[i][ydim][k];
            fi_star[i][ydim + 1][k] = fi_star[i][1][k];
        }
    }

    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            for (int k = 0; k < Q; k++)
            {
                int in = i - ex[k];
                int jn = j - ey[k];
                fi[i][j][k] = fi_star[in][jn][k];
            }
        }
    }
}

void printfi_star()
{
    printf("the fi_star is \n");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("%f ", fi_star[i][j][4]);
        }
        printf("\n");
    }
}

/*
 * compute the macroscopic density rho and velocity u, u = (ux, uy) for every lattice
*/

void computeRho_N_U(double (*rho)[ydim + 2], double (*u_x)[ydim + 2], double (*u_y)[ydim + 2], double (*fi)[ydim + 2][Q])
{
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            rho[i][j] = 0.0;
            u_x[i][j] = 0.0;
            u_y[i][j] = 0.0;
            for (int k = 0; k < Q; k++)
            {
                rho[i][j] += fi[i][j][k];
                u_x[i][j] += ex[k] * fi[i][j][k];
                u_y[i][j] += ey[k] * fi[i][j][k];
            }
            u_x[i][j] /= rho[i][j];
            u_y[i][j] /= rho[i][j];
        }
    }
}

/*
* Compute the equilibrium distribution function f_a^eq for every lattice
* int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
* int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
* double **u_x, double **u_y
*/
void computeFeq(double (*rho)[ydim + 2], double (*u_x)[ydim + 2], double (*u_y)[ydim + 2], double (*feq)[ydim + 2][Q])
{
    double f1 = 3.0;
    double f2 = 9.0 / 2.0;
    double f3 = 3.0 / 2.0;
    double uuDot, euDot, si;

    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            uuDot = u_x[i][j] * u_x[i][j] + u_y[i][j] * u_y[i][j];
            for (int k = 0; k < Q; k++)
            {
                euDot = ex[k] * u_x[i][j] + ey[k] * u_y[i][j];
                si = f1 * euDot + f2 * euDot * euDot - f3 * uuDot; // cVel^2 = 1/3
                // feq[i][j][k] = wi[k] * rho[i][j] + wi[k] * rho[i][j] * si;
                feq[i][j][k] = (1.0 + si) * wi[k] * rho[i][j];
            }
        }
    }
}

/*
* Collide Step : calculate the update distribution function fi = fi_Star - (fi_star - feq) / tau
*/

void collideStep(double (*fi)[ydim + 2][Q], double (*feq)[ydim + 2][Q], double (*fi_star)[ydim + 2][Q], double tau)
{
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            for (int k = 0; k < Q; k++)
            {
                fi_star[i][j][k] = fi[i][j][k] - (fi[i][j][k] - feq[i][j][k]) / tau; // fi_star is the post
            }
        }
    }
}

void printRho(double (*rho)[ydim + 2])
{
    printf("\nthe rho metics is: \n");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("%f ", rho[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void calWholeDensity(double (*fi)[ydim + 2][Q], double mass)
{
    mass = 0.0;
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                mass += fi[i][j][k];
            }
        }
    }
    printf("The total mass is %f \n", mass);
}

void calWholeMomentum(double (*rho)[ydim + 2], double (*u_x)[ydim + 2], double (*u_y)[ydim + 2])
{
    double momentum[2] = {0.0};
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            momentum[0] += rho[i][j] * u_x[i][j];
            momentum[1] += rho[i][j] * u_y[i][j];
        }
    }
    printf("the Momentum in X is %f, in y is %f \n", momentum[0], momentum[1]);
}

void checkCorrectness(int iteration, double (*rho)[ydim + 2], double (*u_x)[ydim + 2], double (*u_y)[ydim + 2], double (*fi)[ydim + 2][Q], double mass)
{
    if (iteration % 1000 == 0)
    {
        calWholeMomentum(rho, u_x, u_y);
        calWholeDensity(fi, mass);
    }
}

/// make the calcualation
///

void calVorticity(double (*u_x)[ydim + 2], double (*u_y)[ydim + 2])
{
    FILE *fp;
    fp = fopen("Vorticity.csv", "w");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            double voriticity = 0.5 * (u_y[i + 1][j] - u_y[i - 1][j]) + 0.5 * (u_x[i][j + 1] - u_x[i][j - 1]);
            fprintf(fp, "%.7f,", voriticity);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void printUxAndUy()
{   
    printf("\n\n");
    for (int i = 1; i < xdim + 1; i++)
    {
        for (int j = 1; j < ydim + 1; j++)
        {
            printf("[%.10f][%.10f] ", u_x[i][j], u_y[i][j]);
        }
        printf("\n");
    }
}

void showTheParameters(int maxiteration, double nu, double tau)
{
    printf("The size of the simulation is %d * %d,\nThe max iteration is %d\nU_0 = 0.01\nK = 30.0\ndelta = 0.05\n", xdim, ydim, maxiteration);
    printf("Simulation Start:\n");
}

float getTime(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec / (float) 1.0e6;
    
}

int main(int argc, char const *argv[])
{

    double big, small;
    int maxiterations = 12000;
    int iterations = 0;
    double nu = 0.000128;
    double tau = (1 + nu * 6) / 2;
    float tstart, tend;
    tstart = getTime();
    initRho_N_U(u_x, u_y, rho);
    initFi(fi, rho);
    computeFeq(rho, u_x, u_y, feq);
    showTheParameters(maxiterations, nu, tau);
    while (iterations < maxiterations)
    {
        computeRho_N_U(rho, u_x, u_y, fi);
        computeFeq(rho, u_x, u_y, feq);
        collideStep(fi, feq, fi_star, tau);
        streamStep(fi_star, fi);
        checkCorrectness(iterations, rho, u_x, u_y, fi, mass);
        iterations++;
    }
    tend = getTime();
    printf("The running time is %f", tend - tstart);
    calVorticity(u_x, u_y);
    return 0;
}
