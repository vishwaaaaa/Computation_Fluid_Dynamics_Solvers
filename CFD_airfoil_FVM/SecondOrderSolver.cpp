//AERO623 Project 2
//First-order Solver

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace Eigen;

MatrixXd ArraytoMatrix(double* A) {
    int n = 4;
    MatrixXd M(1,n);
    for (int i=0; i<n; i++) {
        M(0,i) = (double) A[i];
    }
    return M;
}

double* MatrixtoArray(MatrixXd M) {
    int n = M.cols();
    double* A = new double[n];
    for (int i=0; i<n; i++) {
        A[i] = M(0,i);
    }
    return A;
}

/*double[][] MatrixtoArray2(MatrixXd M, int row, int col) {
    const int m = row;
    const int n = col;
    double A[m][n];
    for (int i=0; i < m; i++) {
        for (int j=0; j<n; j++) {
            A[i][j] = M(i,j);
        }
    }
    return A;
}*/

MatrixXd FileRead(string filename, int nrow, int ncol) {
    
  // Define the matrix to store the data
  MatrixXd matrix(nrow,ncol);

  // Open the input file
  ifstream file(filename);
  if (!file.is_open()) {
    cout << "Error opening file" << endl;
    return matrix;
  }

  // Read the file line by line
  string line;
  int row = 0;
  while (getline(file, line)) {
    // Split the line into separate values
    istringstream iss(line);
    double value;
    int col = 0;
    while (iss >> value) {
      // Resize the matrix if necessary
      /*if (row == matrix.rows()) {
        matrix.conservativeResize(matrix.rows() + 1, Eigen::NoChange);
      }
      if (col == matrix.cols()) {
        matrix.conservativeResize(Eigen::NoChange, matrix.cols() + 1);
      }*/
      // Store the value in the matrix
      matrix(row, col) = value;
      col++;
    }
    row++;
    if (row%1000 == 0) {
        cout << row << "\n";
    }
  }
  /*for (int i=0; i<matrix.rows(); i++) {
    for (int j=0; j<matrix.cols(); j++) {
        if (abs(matrix(i,j)-int(matrix(i,j))) < 1e-6) {
            matrix(i,j) = int(matrix(i,j)+0.5);
        }
    }
  }*/

  // Close the file
  file.close();

  //Output the matrix
  return matrix;
}

void FileWrite(string filename, MatrixXd matrix) {
    ofstream file(filename);

    // Check if file opened successfully
    if (!file.is_open())
    {
        cerr << "Failed to open file" << endl;
    }

    // Write the matrix to the file
    for (int i = 0; i < matrix.rows(); i++)
    {
        for (int j = 0; j < matrix.cols(); j++)
        {
            file << matrix(i, j) << " ";
        }
        file << std::endl;
    }

    // Close the file
    file.close();
}

MatrixXd FreeState() {
    //Calculate Freestream State
    double g = 1.4;
    double Mf = 0.25;
    double alpha = 8*M_PI/180;
    MatrixXd n(1,2);
    n(0,0) = cos(alpha);
    n(0,1) = sin(alpha);
    double a0 = 1;
    double rho0 = 1;
    double p0 = pow(a0,2)*rho0/g;
    double pf = p0*pow((1+(g-1)/2*pow(Mf,2)),-g/(g-1));
    double rhof = rho0*pow((pf/p0),1/g);
    double af = sqrt(g*pf/rhof);
    MatrixXd vf = Mf*af*n;
    MatrixXd vft = vf.transpose();
    MatrixXd vf2 = vf*vft;
    double Ef = (pf/(g-1) + rhof*vf2(0,0)/2);
    MatrixXd U(1,4);
    U(0,0) = rhof; U(0,1) = rhof*vf(0,0); U(0,2) = rhof*vf(0,1); U(0,3) = Ef;
    return U;
}

struct Flux
{
    double F[4];
    double smag;
};

struct Flux RoeFlux(double* UL, double* UR, float gamma, double* n)
{
    // declaring variables:
    double gmi, rL, uL, vL, unL, qL, pL, rHL, HL, cL, FL[4], 
        rR, uR, vR, unR, qR, pR, rHR, HR, cR, FR[4],
        du[4], di, d1, ui, vi, Hi, af, ucp, c2, ci, ci1, l[3], l3, 
        epsilon, s1, s2, G1, G2, C1, C2, F[4], smag;

    gmi = gamma - 1.0;

    // process left state
    rL = UL[0];
    uL = UL[1]/rL;
    vL = UL[2]/rL;
    unL = uL*n[0] + vL*n[1];
    qL = sqrt(pow(UL[1], 2) + pow(UL[2], 2))/rL;
    pL = gmi*(UL[3] - 0.5*rL*pow(qL, 2));

    if ((pL<0) || (rL<0))
    {
        //cout << "Non-physical state!";
        //cout << UL;
    }

    if (pL<0)
    {
        pL = -pL;
    }

    if (rL<0)
    {
        rL = -rL;
    }

    rHL = UL[3] + pL;
    HL = rHL/rL;
    cL = sqrt(gamma*pL/rL);

    // left flux
    FL[0] = rL*unL;
    FL[1] = UL[1]*unL + pL*n[0];
    FL[2] = UL[2]*unL + pL*n[1];
    FL[3] = rHL*unL;

    // process right state
    rR = UR[0];
    uR = UR[1]/rR;
    vR = UR[2]/rR;
    unR = uR*n[0] + vR*n[1];
    qR = sqrt(pow(UR[1],2) + pow(UR[2], 2))/rR;
    pR = gmi*(UR[3] - 0.5*rR*pow(qR,2));

    if ((pR<0) || (rR<0))
    {
        //cout << "Non-physical state!";
        //cout << UR;
    }

    if (pR<0)
    {
        pR = -pR;
    }

    if (rR<0)
    {
        rR = -rR;
    }

    rHR = UR[3] + pR;
    HR = rHR/rR;
    cR = sqrt(gamma*pR/rR);

    // right flux
    FR[0] = rR*unR;
    FR[1] = UR[1]*unR + pR*n[0];
    FR[2] = UR[2]*unR + pR*n[1];
    FR[3] = rHR*unR;

    // difference in states
    
    du[0] = UR[0] - UL[0];
    du[1] = UR[1] - UL[1];
    du[2] = UR[2] - UL[2];
    du[3] = UR[3] - UL[3];

    // Roe average
    di = sqrt(rR/rL);
    d1 = 1.0/(1.0+di);

    ui = (di*uR + uL)*d1;
    vi = (di*vR + vL)*d1;
    Hi = (di*HR + HL)*d1;

    af = 0.5*(ui*ui + vi*vi);
    ucp = ui*n[0] + vi*n[1];
    c2 = gmi*(Hi - af);

    if (c2 < 0)
    {
        //cout << "Non-physical state!";
        c2 = -c2;
    }

    ci = sqrt(c2);
    ci1 = 1.0/ci;

    // eigenvalues
    l[0] = ucp+ci;
    l[1] = ucp-ci;
    l[2] = ucp;

    // entropy fix
    epsilon = ci*0.1;

    for (int i=0; i<3; i++)
    {
        if ((l[i] < epsilon) && (l[i] > -epsilon))
        {
            l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon);
        }
    }

    l[0] = abs(l[0]);
    l[1] = abs(l[1]);
    l[2] = abs(l[2]);
    l[3] = abs(l[3]);
    l3 = l[2];

    // average and half-difference of 1st and 2nd eigenvalues
    s1 = 0.5*(l[0] + l[1]);
    s2 = 0.5*(l[0] - l[1]);

    // left eigenvector product generators
    G1 = gmi*(af*du[0] - ui*du[1] - vi*du[2] + du[3]);
    G2 = -ucp*du[0] + du[1]*n[0] + du[2]*n[1];

    // required functions of G1 and G2
    C1 = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
    C2 = G1*s2*ci1 + G2*(s1-l3);

    // flux assembly
    F[0] = 0.5*(FL[0] + FR[0]) - 0.5*(l3*du[0] + C1);
    F[1] = 0.5*(FL[1] + FR[1]) - 0.5*(l3*du[1] + C1*ui + C2*n[0]);
    F[2] = 0.5*(FL[2] + FR[2]) - 0.5*(l3*du[2] + C1*vi + C2*n[1]);
    F[3] = 0.5*(FL[3] + FR[3]) - 0.5*(l3*du[3] + C1*Hi + C2*ucp);

    // max wave speed
    double smag0 = max(l[0], l[1]);
    double smag1 = max(smag0, l[2] );
    double smag2 = max(smag1,l[3]);
    smag = smag2;
    
    Flux R1;
    for(int i=0; i<4; i++) {
        R1.F[i] = F[i];
    }
    R1.smag = smag;
    return R1;
}

struct Flux RusanovFlux(double* UL, double* UR, float gamma, double* n)
{
    // declaring variables
    double gmi, rho_L, u_L, v_L, unL, v_mag_sq_L, pL, rho_R, u_R, v_R, 
        unR, v_mag_sq_R, pR, FL[4], FR[4], cL, cR, eig1, eig2, smax, F[4], smag;

    gmi = gamma-1.0;

    // left flux
    rho_L = UL[0];
    u_L = UL[1]/rho_L;
    v_L = UL[2]/rho_L;
    unL = u_L*n[0] + v_L*n[1];
    v_mag_sq_L = u_L*u_L + v_L*v_L;
    pL = gmi*(UL[3] - 0.5*rho_L*v_mag_sq_L);
    cL = sqrt(gamma*pL/rho_L);

    FL[0] = rho_L*unL;
    FL[1] = UL[1]*unL + pL*n[0];
    FL[2] = UL[2]*unL + pL*n[1];
    FL[3] = unL*(UL[3] + pL);

    // right flux
    rho_R = UR[0];
    u_R = UR[1]/rho_R;
    v_R = UR[2]/rho_R;
    unR = u_R*n[0] + v_R*n[1];
    v_mag_sq_R = u_R*u_R + v_R*v_R;
    pR = gmi*(UR[3] - 0.5*rho_R*v_mag_sq_R);
    cR = sqrt(gamma*pR/rho_R);

    FR[0] = rho_R*unR;
    FR[1] = UR[1]*unR + pR*n[0];
    FR[2] = UR[2]*unR + pR*n[1];
    FR[3] = unR*(UR[3] + pR);

    // calculating the maximum wave speed of system
    eig1 = unL + cL;
    eig2 = unR + cR;
    smag = max(eig1, eig2);

    // assigning the Flux
    Flux R1;
    for(int i=0; i<4; i++)
    {
        R1.F[i] = 0.5*(FL[i]+FR[i]) - 0.5*smag*(UR[i] - UL[i]);
    }

    R1.smag = smag;

    return R1;

}

struct Flux HLLEFlux(double* UL, double* UR, float gamma, double* n)
{
    // declaring variables
    double gmi, rho_L, u_L, v_L, unL, v_mag_sq_L, pL, rho_R, u_R, v_R, 
        unR, v_mag_sq_R, pR, FL[4], FR[4], cL, cR, eig1, eig2, smax, 
        smin, sLmax, sLmin, sRmax, sRmin, F[4], smag;

    gmi = gamma-1.0;

    // left flux
    rho_L = UL[0];
    u_L = UL[1]/rho_L;
    v_L = UL[2]/rho_L;
    unL = u_L*n[0] + v_L*n[1];
    v_mag_sq_L = u_L*u_L + v_L*v_L;
    pL = gmi*(UL[3] - 0.5*rho_L*v_mag_sq_L);
    cL = sqrt(gamma*pL/rho_L);

    FL[0] = rho_L*unL;
    FL[1] = UL[1]*unL + pL*n[0];
    FL[2] = UL[2]*unL + pL*n[1];
    FL[3] = unL*(UL[3] + pL);

    // right flux
    rho_R = UR[0];
    u_R = UR[1]/rho_R;
    v_R = UR[2]/rho_R;
    unR = u_R*n[0] + v_R*n[1];
    v_mag_sq_R = u_R*u_R + v_R*v_R;
    pR = gmi*(UR[3] - 0.5*rho_R*v_mag_sq_R);
    cR = sqrt(gamma*pR/rho_R);

    FR[0] = rho_R*unR;
    FR[1] = UR[1]*unR + pR*n[0];
    FR[2] = UR[2]*unR + pR*n[1];
    FR[3] = unR*(UR[3] + pR);

    // maximum and minimum wave speeds
    sLmin = min(double(0), (unL - cL));
    sRmin = min(double(0), (unR - cR));
    smin = min(sLmin, sRmin);

    sLmax = max(double(0), (unL + cL));
    sRmax = max(double(0), (unR + cR));
    smax = max(sLmax, sRmax);

    // calculating the flux
    Flux R1;
    for(int i=0; i<4; i++)
    {
        R1.F[i] = 0.5*(FL[i] + FR[i]) - 0.5*((smax + smin)/(smax-smin))*(FR[i] - FL[i]) + ((smax*smin)/(smax-smin))*(UR[i]-UL[i]);
    }
    R1.smag = max((abs(unL)+cL),(abs(unR)+cR)); 

    return R1;
}

struct Flux WallFlux(double* u_plus, double gamma, double* n)
{

    // declaring variables
    double gmi, F[4], v[2], vb[2], vb_mag_sq, pb, rho_b, cb, smag;

    gmi = gamma - 1.0;

    // velocity components next to boundary
    v[0] = u_plus[1]/u_plus[0];
    v[1] = u_plus[2]/u_plus[0];

    // velocity components at boundary
    vb[0] = v[0] - (v[0]*n[0] + v[1]*n[1])*n[0];
    vb[1] = v[1] - (v[0]*n[0] + v[1]*n[1])*n[1];
    vb_mag_sq = vb[0]*vb[0] + vb[1]*vb[1];

    // pressure at boundary
    pb = gmi*(u_plus[3] - 0.5*u_plus[0]*vb_mag_sq);

    // flux calculation normal to boundary
    Flux R1;
    R1.F[0] = 0;
    R1.F[1] = pb*n[0];
    R1.F[2] = pb*n[1];
    R1.F[3] = 0;

    // wave speed calculation
    rho_b = u_plus[0];
    cb = sqrt(gamma*pb/rho_b);
    R1.smag = cb + abs(vb[0]*n[0] + vb[1]*n[1]);

    return R1;

}

MatrixXd ResidCalc(MatrixXd uS, MatrixXd S) {
    double g = 1.4;
    MatrixXd Uinf = FreeState();
    MatrixXd R(uS.rows()/3,uS.cols());
    R.setZero();
    MatrixXd F(S.rows(),uS.cols());
    MatrixXd s(uS.rows()/3,1);
    s.setZero();
    MatrixXd W(S.rows(),1);
    Flux FW;
    for (int i=0; i<S.rows(); i++) {
        int k = floor(i/3);
        if (S(i,1) > 0) {
            double* uL = MatrixtoArray(uS(i,all));
            double* uR = MatrixtoArray(uS(int(S(i,6)+.5)-1,all));
            double* n = MatrixtoArray(S(i,seq(2,3)));
            if (int(S(i,5)+.5) == 0) {
                FW = RoeFlux(uL,uR,g,n);
                F(i,all) = ArraytoMatrix(FW.F);
                W(i) = (double) FW.smag;
            } else {
                F(i,all) = -F(int(S(i,6)+.5)-1,all);
                W(i) = W(int(S(i,6)+.5)-1);
            }
        } else if (S(i,1) == -1) {
            double* uL = MatrixtoArray(uS(i,all));
            double* uR = MatrixtoArray(Uinf);
            double* n = MatrixtoArray(S(i,seq(2,3)));
            FW = RoeFlux(uL,uR,g,n);
            F(i,all) = ArraytoMatrix(FW.F);
            W(i) = (double) FW.smag;
        } else {
            double* uL = MatrixtoArray(uS(i,all));
            double* n = MatrixtoArray(S(i,seq(2,3)));
            FW = WallFlux(uL,g,n);
            F(i,all) = ArraytoMatrix(FW.F);
            W(i) = (double) FW.smag;
        }
        R(k,all) = R(k,all) + F(i,all)*S(i,4);
        s(k) = s(k) + W(i)*S(i,4);
        //cout << k << "\n";
    }
    MatrixXd Rs(R.rows(),R.cols()+1);
    Rs << R,s;
    /*cout << "State \n"  << uS << "\n";
    cout << "Flux \n" << F << "\n";
    cout << "Residual \n" << R << "\n";*/
    return Rs;
}

MatrixXd EdgeState(MatrixXd u, MatrixXd grx, MatrixXd gry, MatrixXd Z, MatrixXd M) {
    MatrixXd uS(M.rows(),4);
    for (int i=0; i<uS.rows(); i++) {
        int k = floor(i/3);
        uS(i,all) = u(k,all) + grx(k,all)*(M(i,0)-Z(k,0)) + gry(k,all)*(M(i,1)-Z(k,1));
    }
    return uS;
}

double min(double arr[], int size)
{
    double min = arr[0];
    for (int i = 1; i < size; i++)
    {
        if (arr[i] < min) min = arr[i];
    }

    return min;
}

double max(double arr[], int size)
{
    double max = arr[0];
    for (int i = 1; i < size; i++)
    {
        if (arr[i] > max) max = arr[i];
    }

    return max;
}

void centroid(int dim, double p1[], double p2[], double p3[], double pc[])
{
    for (int i = 0; i < dim; i++)
    {
        pc[i] = (p1[i] + p2[i] + p3[i]) / 3;
    }
}

void vect(int dim, double p1[], double p2[], double r12[], double rc12[], double mag)
{
    double sum = 0;
    for (int i = 0; i < dim; i++)
    {
        r12[i] = p2[i] - p1[i];
        sum += r12[i] * r12[i];
    }
    mag = sqrt(sum);
    for (int i = 0; i < dim; i++)
    {
        rc12[i] = r12[i]/mag;
    }
}

void normal_vect(int dim, double p1[], double p2[], double r12[])
{
    if (dim == 2)
    {
        double sum = sqrt(pow(p2[0] - p1[0],2.)+ pow(p2[1] - p1[1], 2.));
        r12[0] = (p2[1] - p1[1])/sum;
        r12[1] = (p1[0] - p2[0])/sum;
    }
    else
    {
        for (int i = 0; i < dim; i++) r12[i] = 0;
    }
}

double dotprod(int dim, double p1[], double p2[])
{
    double sum = 0;
    for (int i = 0; i < dim; i++)
    {
        sum += p1[i]*p2[i];
    }
    return sum;
}

void crossprod(int dim, double p1[], double p2[], double crossprod[])
{
    if (dim == 3)
    {
        crossprod[0] = (p1[1] * p2[2]) - (p1[2] * p2[1]);
        crossprod[1] = (p1[2] * p2[0]) - (p1[0] * p2[2]);
        crossprod[2] = (p1[0] * p2[1]) - (p1[1] * p2[0]);
    }
    else
    {
        for (int i = 0; i < dim; i++) crossprod[i] = 0;
    }
}

void gradient_b(int state, int N_elem, int E2N[][3], int EnE[][3], double u[], double node_coord[][2], double grad_u[][2]) {
    MatrixXd uinf = FreeState();
    for (int i=0; i<N_elem; i++) {
        for (int j=0; j<2; j++) {
            grad_u[i][j] = 0;
        }
        for (int k=0; k<3; k++) {
            double p1[2] = {node_coord[E2N[i][(k+1)%3]-1][0],node_coord[E2N[i][(k+1)%3]-1][1]};
            double p2[2] = {node_coord[E2N[i][(k+2)%3]-1][0],node_coord[E2N[i][(k+2)%3]-1][1]};
            double n[2];
            double dl = sqrt(pow((p2[1]-p1[1]),2)+pow((p2[0]-p1[0]),2));
            normal_vect(2,p1,p2,n);
            double uh;
            if (EnE[i][k] > 0) {
                uh = (u[i]+u[EnE[i][k]-1])/2;
            } else if (EnE[i][k] == -1) {
                uh = (u[i]+uinf(state))/2;
            } else {
                uh = u[i];
            }
            grad_u[i][0] = grad_u[i][0] + uh*n[0]*dl;
            //cout << uh << " " << n[0] << " " << dl << "\n";
            grad_u[i][1] = grad_u[i][1] + uh*n[1]*dl;
        }
        double v1[3], v2[3], v3[3], A;
        v1[0] = node_coord[E2N[i][1]-1][0]-node_coord[E2N[i][0]-1][0];
        v1[1] = node_coord[E2N[i][1]-1][1]-node_coord[E2N[i][0]-1][1];
        v2[0] = node_coord[E2N[i][2]-1][0]-node_coord[E2N[i][0]-1][0];
        v2[1] = node_coord[E2N[i][2]-1][1]-node_coord[E2N[i][0]-1][1];
        crossprod(3,v1,v2,v3);
        A = (v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2])/2;
        //cout << A << "\n";
        grad_u[i][0] = grad_u[i][0]/A;
        grad_u[i][1] = grad_u[i][1]/A;
    }
}

void gradient_pn(int N_elem, int EnE[][3], double u[], double Ecent[][2], double grad_u[][2])
{
    // gradient calculated using plane normal method
    int x = 0;

    double V1[3], v1[3], V2[3], v2[3], N[3],p1[3],p2[3],p3[3];
    double mag_V1 = 0;
    double mag_V2 = 0;
    for (int i = 0; i < N_elem; i++)
    {
        if (EnE[i][0] < 0)
        {
            double p1[3] = { Ecent[i][0],Ecent[i][1],u[i] };
            double p2[3] = { Ecent[EnE[i][1] - 1][0],Ecent[EnE[i][1] - 1][1],u[EnE[i][1] - 1] };
            double p3[3] = { Ecent[EnE[i][2] - 1][0],Ecent[EnE[i][2] - 1][1],u[EnE[i][2] - 1] };
            vect(3, p1, p2, V1, v1, mag_V1);
            vect(3, p3, p1, V2, v2, mag_V2);
            crossprod(3, V1, V2, N);
            cout << V1[0] << " " << V1[1] << " " << V1[2] << "\n";
            cout << V2[0] << " " << V2[1] << " " << V2[2] << "\n";
            grad_u[i][0] = -1 * N[0] / N[2];
            grad_u[i][1] = -1 * N[1] / N[2];
        }
        else if (EnE[i][1] < 0)
        {
            double p1[3] = { Ecent[EnE[i][0] - 1][0],Ecent[EnE[i][0] - 1][1],u[EnE[i][0] - 1] };
            double p2[3] = { Ecent[i][0],Ecent[i][1],u[i] };
            double p3[3] = { Ecent[EnE[i][2] - 1][0],Ecent[EnE[i][2] - 1][1],u[EnE[i][2] - 1] };
            vect(3, p1, p2, V1, v1, mag_V1);
            vect(3, p3, p1, V2, v2, mag_V2);
            crossprod(3, V1, V2, N);

            grad_u[i][0] = -1 * N[0] / N[2];
            grad_u[i][1] = -1 * N[1] / N[2];

        }
        else if (EnE[i][2] < 0)
        {
            double p1[3] = { Ecent[EnE[i][0] - 1][0],Ecent[EnE[i][0] - 1][1],u[EnE[i][0] - 1] };
            double p2[3] = { Ecent[EnE[i][1] - 1][0],Ecent[EnE[i][1] - 1][1],u[EnE[i][1] - 1] };
            double p3[3] = { Ecent[i][0],Ecent[i][1],u[i] };
            vect(3, p1, p2, V1, v1, mag_V1);
            vect(3, p3, p1, V2, v2, mag_V2);
            crossprod(3, V1, V2, N);

            grad_u[i][0] = -1 * N[0] / N[2];
            grad_u[i][1] = -1 * N[1] / N[2];

        }
        else
        {
            double p1[3] = { Ecent[EnE[i][0]-1][0],Ecent[EnE[i][0]-1][1],u[EnE[i][0]-1] };
            double p2[3] = { Ecent[EnE[i][1]-1][0],Ecent[EnE[i][1]-1][1],u[EnE[i][1]-1] };
            double p3[3] = { Ecent[EnE[i][2]-1][0],Ecent[EnE[i][2]-1][1],u[EnE[i][2]-1] };      
            vect(3, p1, p2, V1, v1, mag_V1);
            vect(3, p3, p1, V2, v2, mag_V2);
            crossprod(3, V1, V2, N);

            grad_u[i][0] = -1 * N[0] / N[2];
            grad_u[i][1] = -1 * N[1] / N[2];
        }
        
    }
    
}

void BarthJesperson(int state, int dim, int n_elem, int E2N[][3], double node_coord[][2], int EnE[][3], double u[], double Ecent[][2], double grad_u[][2], double alpha[], double grad_u_lim[][2])
{
    double u_local[4];
    double alpha_iN[3], uiN[3];;
    double r1N[2], r2N[2], r3N[2];
    double ur1N[2], ur2N[2], ur3N[2];
    double l_r1N = 0;
    double l_r2N = 0;
    double l_r3N = 0;

    // Element centroid calculation

    for (int i = 0; i < n_elem; i++)
    {
        int p1 = E2N[i][0] - 1;
        int p2 = E2N[i][1] - 1;
        int p3 = E2N[i][2] - 1;
        double a[2] = { node_coord[p1][0],node_coord[p1][1] };
        double b[2] = { node_coord[p2][0],node_coord[p2][1] };
        double c[2] = { node_coord[p3][0],node_coord[p3][1] };
        double d[2];
        centroid(dim, a, b, c, d);
        Ecent[i][0] = d[0];
        Ecent[i][1] = d[1];
    }
    gradient_b(state, n_elem, E2N, EnE, u, node_coord, grad_u);

    for (int i = 0; i < n_elem; i++)
    {
        double cp0[2] = { Ecent[i][0] ,Ecent[i][1] };
        
        u_local[0] = u[i];
        for (int j = 0; j < 3; j++)
        {
            if (EnE[i][j] < 0) u_local[j + 1] = u[i];
            else u_local[j + 1] = u[EnE[i][j] - 1];
        }
        double u_max = max(u_local, 4);
        double u_min = min(u_local, 4);
        
        int p1 = E2N[i][0] - 1;
        int p2 = E2N[i][1] - 1;
        int p3 = E2N[i][2] - 1;

        double n1[2] = { node_coord[p1][0] , node_coord[p1][1] };
        double n2[2] = { node_coord[p2][0] , node_coord[p2][1] };
        double n3[2] = { node_coord[p3][0] , node_coord[p3][1] };

        vect(2, cp0, n1, r1N, ur1N, l_r1N);
        vect(2, cp0, n2, r2N, ur2N, l_r2N);
        vect(2, cp0, n3, r3N, ur3N, l_r3N);

        double grad_v[2] = { grad_u[i][0],grad_u[i][1]};

        uiN[0] = u[i] + dotprod(2, r1N, grad_v);
        uiN[1] = u[i] + dotprod(2, r2N, grad_v);
        uiN[2] = u[i] + dotprod(2, r3N, grad_v);

        for (int j = 0; j < 3; j++)
        {
            if (uiN[j] > u[i])
            {
                double comp[2] = {1,(u_max-u[i])/(uiN[j]-u[i])};
                alpha_iN[j]= min(comp, 2);

            }
            else if (uiN[j] < u[i])
            {
                double comp[2] = { 1,(u_min - u[i]) / (uiN[j] - u[i]) };
                alpha_iN[j] = min(comp, 2);
            }
            else alpha_iN[j] = 1;
        }
        alpha[i] = min(alpha_iN, 3);
        grad_u_lim[i][0] = grad_u[i][0] * alpha[i];
        grad_u_lim[i][1] = grad_u[i][1] * alpha[i];
    }
}

vector< vector<double> > mpscale(vector< vector<double> > nodes, vector< vector<int> > elem, vector< vector<int> > conn, vector< vector<double> > state, int stateFlag) {
    int numCells = state.size(); // DOUBLE CHECK THIS CORRECTLY INDEXS NUMBER OF CELLS, AND ISNT COUNTING COLUMNS (i.e. should NOT be 4)
    vector< vector<double> > Lout(numCells, vector<double>(2));
    for (int i = 0; i < numCells; i++) {
        // ---------------------- Gradient Calculations -------------------- //
        
        // - Geometry Set Up - 
        vector<int> side(3);
        side[0] = elem[i][0];
        side[1] = elem[i][1];
        side[2] = elem[i][2];
        
        vector<double> x(3), y(3);
        x[0] = nodes[side[0]-1][0]; y[0] = nodes[side[0]-1][1];
        x[1] = nodes[side[1]-1][0]; y[1] = nodes[side[1]-1][1];
        x[2] = nodes[side[2]-1][0]; y[2] = nodes[side[2]-1][1];
        
        vector<double> l(3);
        l[0] = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2));
        l[1] = sqrt(pow(x[2]-x[1], 2) + pow(y[2]-y[1], 2));
        l[2] = sqrt(pow(x[0]-x[2], 2) + pow(y[0]-y[2], 2));
        
        // - Normal Calculations - 
        vector< vector<double> > n(3, vector<double>(2));
        n[0][0] = y[1] - y[0]; n[0][1] = x[0] - x[1];
        n[1][0] = y[2] - y[1]; n[1][1] = x[1] - x[2];
        n[2][0] = y[0] - y[2]; n[2][1] = x[2] - x[0];
        
        // - Semiperimeter and Area Calculation from Heron's Formula - 
        double s = (l[0] + l[1] + l[2]) / 2;
        double A = sqrt(s * (s-l[0]) * (s-l[1]) * (s-l[2]));
        
        // State Vector Generation
        double u0 = state[i][stateFlag];
        std::vector<double> u(3);
        for (int j = 0; j < 3; j++) {
            int connindex = conn[i][j];
            if (conn[i][j] >= 1) {
                u[j] = state[connindex-1][stateFlag];
            } else {
                u[j] = u0;
            }
        }

        // Find Cell Index Across Side
        std::vector<int> cellAcross(3);
        for (int j = 0; j < 3; j++) {
        cellAcross[j] = conn[i][j];
        if (cellAcross[j] < 1) {
            std::vector<int> connData(3);
            connData = {conn[i][0], conn[i][1], conn[i][2]};
            cellAcross[j] = max({conn[i][0], conn[i][1], conn[i][2]}); // Disregards boundary cell and double counts freestream cell
            // THIS IS A SIMPLE APPROXIMATION. A MORE PRECISE SOLUTION WOULD BE TO DISREGARD THE CELL ENTIRELY AND ONLY USE 2 (OR 1) CELLS
            // IF WE SEE CONVERGENCE ISSUES THIS CAN BE MODIFIED AS NEEDED
        }
        }

        // Coordinates for Other Cell
        std::vector<std::vector<double>> xs(3, std::vector<double>(3)),
            ys(3, std::vector<double>(3));
        for (int j = 0; j < 3; j++) {
        xs[j][0] = nodes[elem[cellAcross[j]-1][0]-1][0];
        ys[j][0] = nodes[elem[cellAcross[j]-1][0]-1][1];
        xs[j][1] = nodes[elem[cellAcross[j]-1][1]-1][0];
        ys[j][1] = nodes[elem[cellAcross[j]-1][1]-1][1];
        xs[j][2] = nodes[elem[cellAcross[j]-1][2]-1][0];
        ys[j][2] = nodes[elem[cellAcross[j]-1][2]-1][1];
        }

        // Centroid Calculations
        std::vector<std::vector<double>> G(4, std::vector<double>(2));
        for (int j = 0; j < 3; j++) {
        G[j][0] = (xs[j][0] + xs[j][1] + xs[j][2]) / 3; // Gx
        G[j][1] = (ys[j][0] + ys[j][1] + ys[j][2]) / 3; // Gy
        }
        G[3][0] = (x[0] + x[1] + x[2]) / 3;
        G[3][1] = (y[0] + y[1] + y[2]) / 3;

        // r Vectors to Edge Midpoints
        std::vector<std::vector<double>> mid(3, std::vector<double>(2));
        mid[0][0] = (x[0] + x[1]) / 2;
        mid[0][1] = (y[0] + y[1]) / 2;
        mid[1][0] = (x[1] + x[2]) / 2;
        mid[1][1] = (y[2] + y[1]) / 2;
        mid[2][0] = (x[2] + x[0]) / 2;
        mid[2][1] = (y[2] + y[0]) / 2;

        vector<vector<double>> r0k(3, vector<double>(2));
        for (int j = 0; j < 3; j++) {
            r0k[j][0] = G[3][0] - mid[j][0];
            r0k[j][1] = G[3][1] - mid[j][1];
        }
        
        // Cell Gradient Calculation - Gradient Theorem
        vector<double> L(2);
        L[0] = 0; L[1] = 0;
        for (int j = 0; j < 3; j++) {
            L[0] += ((u0 + u[j]))*n[j][0]/2;
            L[1] += ((u0 + u[j]))*n[j][1]/2;
        }
        if (i == 0) {
            //cout << u0 << " " << u[0] << " " << u[1] << " " << u[2] << "\n";
            //cout << n[0][0] << " " << n[1][0] << " " << n[2][0] << " " << L[0] << "\n";
        }
        L[0] /= A;
        L[1] /= A;

        // Alpha Parameter Calculation
        vector<double> alphak(3);
        for (int j = 0; j < 3; j++) {
            if (r0k[j][0]*L[0] + r0k[j][1]*L[1] > max(u[j] - u0, double(0))) {
            alphak[j] = max(u[j] - u0, double(0))/(r0k[j][0]*L[0] + r0k[j][1]*L[1]);
            } else if (r0k[j][0]*L[0] + r0k[j][1]*L[1] < min(u[j] - u0, double(0))) {
            alphak[j] = min(u[j] - u0, double(0))/(r0k[j][0]*L[0] + r0k[j][1]*L[1]);
            } else {
            alphak[j] = 1;
            }
        }
        /*if (i==0) {
            cout << alphak[0] << " " << alphak[1] << " " << alphak[2] << "\n";
            for (int j = 0; j < 3; j++) {
                cout << max(u[j] - u0, double(0))/(r0k[j][0]*L[0] + r0k[j][1]*L[1]) << "\n";
                cout << min(u[j] - u0, double(0))/(r0k[j][0]*L[0] + r0k[j][1]*L[1]) << "\n";
                cout << (r0k[j][0]*L[0] + r0k[j][1]*L[1]) << "\n";
                cout << u0 << " " << u[j] << "\n";
                cout << "j = " << j << "\n";
            }
        } */
        
        double alpha = min({alphak[0], alphak[1], alphak[2]});
        
        // Gradient Scaling
        for (int k = 0; k < 2; k++) {
        Lout[i][k] = L[k]*alpha;
        }
    }
    return Lout;
}

MatrixXd GradCalc(MatrixXd u, MatrixXd E, MatrixXd N, MatrixXd C) {
    MatrixXd gr(u.rows(),2*u.cols());
    int n_elem = u.rows();
    int n_nodes = N.rows();
    //I know this is annoying, but there isn't an easier way to use limiters with different variable types
    //Uncomment this for BJ Limiter:
    /*
    int dim = 2;
    int E2N[n_elem][3];
    for (int i=0; i<n_elem; i++) {
        for (int j=0; j<3; j++) {
            E2N[i][j]= (int) E(i,j);
        }
    }
    double node_coord[N.rows()][2];
    for (int i=0; i<N.rows(); i++) {
        for (int j=0; j<2; j++) {
            node_coord[i][j] = N(i,j);
        }
    }
    int EnE[n_elem][3];
    for (int i=0; i<n_elem; i++) {
        for (int j=0; j<3; j++) {
            EnE[i][j] = (int) C(i,j);
        }
    }
    double Ecent[n_elem][2];
    for (int k=0; k<4; k++) {
        double u1[n_elem];
        double grad_u[n_elem][2];
        double alpha[n_elem];    
        double grad_u_lim[n_elem][2];
        for (int i=0; i<n_elem; i++) {
            u1[i] = u(i,k);
        }
        BarthJesperson(k,dim,n_elem,E2N,node_coord,EnE,u1,Ecent,grad_u,alpha,grad_u_lim);
        for (int i=0; i<n_elem; i++) {
            gr(i,k) = grad_u_lim[i][0];
            gr(i,k+4) = grad_u_lim[i][1];
        }

    }*/
    //Uncomment This for MP Limiter
    vector<vector<int>> elem(n_elem, vector<int>(3));
    for (int i=0; i<n_elem; i++) {
        for (int j=0; j<3; j++) {
            elem[i][j] = (int) E(i,j);
        }
    }
    vector<vector<double>> nodes(n_nodes, vector<double>(2));
    for (int i=0; i<n_nodes; i++) {
        for (int j=0; j<2; j++) {
            nodes[i][j] = N(i,j);
        }
    }
    vector<vector<int>> conn(n_elem, vector<int>(3));
    for (int i=0; i<n_elem; i++) {
        for (int j=0; j<3; j++) {
            conn[i][j] = (int) C(i,j);
        }
    }
    vector<vector<double> > state(n_elem, vector<double>(4));
    for (int i=0; i<n_elem; i++) {
        for (int j=0; j<4; j++) {
            state[i][j] = u(i,j);
        }
    }
    vector<vector<double> > grad_u_lim(n_elem, vector<double>(2));
    for (int k=0; k<4; k++) {
        grad_u_lim = mpscale(nodes,elem,conn,state,k);
        for (int i=0; i<n_elem; i++) {
            gr(i,k) = grad_u_lim[i][0];
            gr(i,k+4) = grad_u_lim[i][1];
        }
    }
    return gr;
}

int main() {
    //Load Mesh Files
    MatrixXd N = FileRead("CoarseMeshNodes.txt",1105,2);
    MatrixXd E = FileRead("CoarseMeshElems.txt",2054,3);
    MatrixXd C = FileRead("CoarseMeshConnect.txt",2054,3);
    MatrixXd U0 = FileRead("CoarseMeshSol1.txt",2054,4);
    //MatrixXd N = FileRead("FineMeshNodes.txt",4266,2);
    //MatrixXd E = FileRead("FineMeshElems.txt",8216,3);
    //MatrixXd C = FileRead("FineMeshConnect.txt",8216,3);
    //MatrixXd U0 = FileRead("FineMeshSol1.txt",8216,4);
    //MatrixXd N = FileRead("FinerMeshNodes.txt",16750,2);
    //MatrixXd E = FileRead("FinerMeshElems.txt",32864,3);
    //MatrixXd C = FileRead("FinerMeshConnect.txt",32864,3);
    //MatrixXd U0 = FileRead("FinerMeshSol1.txt",32864,4);
    //MatrixXd N = FileRead("FinerMeshNodes.txt",66366,2);
    //MatrixXd E = FileRead("FinerMeshElems.txt",131456,3);
    //MatrixXd C = FileRead("FinerMeshConnect.txt",131456,3);
    //MatrixXd U0 = FileRead("FineMeshSol1.txt",131456,4);
    //MatrixXd N = FileRead("UnitMeshN.txt",4,2);
    //MatrixXd E = FileRead("UnitMeshE.txt",2,3);
    //MatrixXd C = FileRead("UnitMeshC.txt",2,3);
    //MatrixXd U0 = FileRead("UnitMeshSol1.txt",2,4);
    cout << "Read Files \n";
    
    //Find Edges, Normals, and Midpoints
    MatrixXd S(3*C.rows(),7);
    MatrixXd M(S.rows(),2);
    S.isZero();
    int k = 0;
    for (int i = 0; i < C.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int j1 = (j+2)%3;
            int j2 = (j+1)%3;
            S(k,0) = i+1;
            S(k,1) = C(i,j);
            S(k,2) = N(int(E(i,j1)+.5)-1,1)-N(int(E(i,j2)+.5)-1,1);
            S(k,3) = N(int(E(i,j2)+.5)-1,0)-N(int(E(i,j1)+.5)-1,0);
            S(k,4) = sqrt(pow(S(k,2),2)+pow(S(k,3),2));
            S(k,2) = S(k,2)/S(k,4);
            S(k,3) = S(k,3)/S(k,4);
            S(k,5) = 0;
            M(k,1) = (N(int(E(i,j1)+.5)-1,1)+N(int(E(i,j2)+.5)-1,1))/2;
            M(k,0) = (N(int(E(i,j2)+.5)-1,0)+N(int(E(i,j1)+.5)-1,0))/2;
            k = k + 1;
        }
    }
    MatrixXd v(S.rows(),1);
    for (int i = 0; i < S.rows(); i++) {
        if (S(i,1) > 0) {
            int t = int(S(i,1)+0.5);
            if (S(3*(t-1)+2,1) == S(i,0)) {
                v(i) = 3*(t-1)+3;
                if (S(i,5) == 0) {
                    S(3*(t-1)+2,5) = 1;
                }
            } else if (S(3*(t-1)+1,1) == S(i,0)) {
                v(i) = 3*(t-1)+2;
                if (S(i,5) == 0) {
                    S(3*(t-1)+1,5) = 1;
                }
            } else {
                v(i) = 3*(t-1)+1;
                if (S(i,5) == 0) {
                    S(3*(t-1),5) = 1;
                }
            }
        }
    }
    S(all,6) = v;
    cout << "Built S \n";

    //Find Cell Centroids
    MatrixXd Z(E.rows(),2);
    for (int i=0; i<E.rows(); i++) {
        Z(i,0) = (N(int(E(i,0)+.5)-1,0)+N(int(E(i,1)+.5)-1,0)+N(int(E(i,2)+.5)-1,0))/3;
        Z(i,1) = (N(int(E(i,0)+.5)-1,1)+N(int(E(i,1)+.5)-1,1)+N(int(E(i,2)+.5)-1,1))/3;
    }

    //Calculate Freestream State
    MatrixXd Uinf = FreeState();

    //Initialize Solver
    MatrixXd u(E.rows(),4);
    //for (int i=0; i<u.rows(); i++) {
    //    u(i,all) = Uinf;
    //}
    u = U0;

    //Converge Solution
    double err = 100;
    int c = 0;
    double CFL = 1;
    cout << "Starting Loop \n";
    auto start = chrono::high_resolution_clock::now();
    while (err > 1e-5 && c < 1000) {
        MatrixXd grx(E.rows(),4);
        MatrixXd gry(E.rows(),4);
        MatrixXd gr = GradCalc(u,E,N,C);
        MatrixXd du0(E.rows(),4);
        MatrixXd du1(E.rows(),4);
        //cout << gr << "\n";
        grx = gr(all,seq(0,3));
        gry = gr(all,seq(4,7));
        MatrixXd uS = EdgeState(u,grx,gry,Z,M);
        MatrixXd Rs = ResidCalc(uS,S);
        MatrixXd R = Rs(all,seq(0,3));
        MatrixXd s = Rs(all,4);
        MatrixXd s4(s.rows(),4);
        s4 = s.replicate(1,4);
        MatrixXd uh(E.rows(),4);
        du0 = -2*CFL*R.cwiseQuotient(s4);
        uh = u + du0;
        gr = GradCalc(uh,E,N,C);
        grx = gr(all,seq(0,3));
        gry = gr(all,seq(4,7));
        uS = EdgeState(uh,grx,gry,Z,M);
        Rs = ResidCalc(uS,S);
        R = Rs(all,seq(0,3));
        du1 = -2*CFL*R.cwiseQuotient(s4);
        u = u + (du0 + du1)/2;
        R = -(du0 + du1)/(2*CFL);
        R = R.cwiseProduct(s4);
        c = c + 1;
        err = (R.cwiseAbs()).sum();
        if (c%10 == 0) {
        cout << "Iteration Number: " << c << "\n";
        cout << "Residual Norm: " << err << "\n";
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<chrono::milliseconds>(stop - start);
    double dt = time.count();
    dt = dt/1000;
    cout << "Time to Reach Solution: " << dt << " seconds \n";
    cout << "Iteration Number: " << c << "\n";
    cout << "Residual Norm: " << err << "\n";
    FileWrite("UnitMeshSol2.txt",u);
    return 0;
}

