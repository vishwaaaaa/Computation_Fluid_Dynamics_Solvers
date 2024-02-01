//AERO623 Project 2
//First-order Solver

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace Eigen;

MatrixXd ArraytoMatrix(double* A) {
    int n = 4; //For some reason A is claiming to be half as big as it actually is
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
        cout << "Non-physical state!";
        cout << UL;
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
        cout << "Non-physical state!";
        cout << UR;
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
        cout << "Non-physical state!";
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

MatrixXd ResidCalc(MatrixXd u, MatrixXd S) {
    double g = 1.4;
    MatrixXd Uinf = FreeState();
    MatrixXd R(u.rows(),u.cols());
    R.setZero();
    MatrixXd F(S.rows(),u.cols());
    MatrixXd s(u.rows(),1);
    s.setZero();
    MatrixXd W(S.rows(),1);
    Flux FW;
    for (int i=0; i<S.rows(); i++) {
        int k = floor(i/3);
        if (S(i,1) > 0) {
            double* uL = MatrixtoArray(u(int(S(i,0)+.5)-1,all));
            double* uR = MatrixtoArray(u(int(S(i,1)+.5)-1,all));
            double* n = MatrixtoArray(S(i,seq(2,3)));
            if (int(S(i,5)+.5) == 0) {
                FW = RoeFlux(uL,uR,g,n);
                F(i,all) = ArraytoMatrix(FW.F);
                W(i) = (double) FW.smag;
            } else {
                F(i,all) = -F(int(S(i,5)+.5)-1,all);
                W(i) = W(int(S(i,5)+.5)-1);
            }
        } else if (S(i,1) == -1) {
            double* uL = MatrixtoArray(u(int(S(i,0)+.5)-1,all));
            double* uR = MatrixtoArray(Uinf);
            double* n = MatrixtoArray(S(i,seq(2,3)));
            FW = RoeFlux(uL,uR,g,n);
            F(i,all) = ArraytoMatrix(FW.F);
            W(i) = (double) FW.smag;
        } else {
            double* uL = MatrixtoArray(u(int(S(i,0)+.5)-1,all));
            double* n = MatrixtoArray(S(i,seq(2,3)));
            FW = WallFlux(uL,g,n);
            F(i,all) = ArraytoMatrix(FW.F);
            W(i) = (double) FW.smag;
        }
        R(k,all) = R(k,all) + F(i,all)*S(i,4);
        s(k) = s(k) + W(i)*S(i,4);
        //cout << k << "\n";
    }
    MatrixXd Rs(u.rows(),u.cols()+1);
    Rs << R,s;
    return Rs;
}

int main() {
    //Load Mesh Files
    MatrixXd N = FileRead("CoarseMeshNodes.txt",1105,2);
    MatrixXd E = FileRead("CoarseMeshElems.txt",2054,3);
    MatrixXd C = FileRead("CoarseMeshConnect.txt",2054,3);
    //MatrixXd N = FileRead("FineMeshNodes.txt",4266,2);
    //MatrixXd E = FileRead("FineMeshElems.txt",8216,3);
    //MatrixXd C = FileRead("FineMeshConnect.txt",8216,3);
    //MatrixXd N = FileRead("FinerMeshNodes.txt",16750,2);
    //MatrixXd E = FileRead("FinerMeshElems.txt",32864,3);
    //MatrixXd C = FileRead("FinerMeshConnect.txt",32864,3);
    //MatrixXd N = FileRead("FinerMeshNodes.txt",66366,2);
    //MatrixXd E = FileRead("FinerMeshElems.txt",131456,3);
    //MatrixXd C = FileRead("FinerMeshConnect.txt",131456,3);
    //MatrixXd N = FileRead("UnitMeshN.txt",4,2);
    //MatrixXd E = FileRead("UnitMeshE.txt",2,3);
    //MatrixXd C = FileRead("UnitMeshC.txt",2,3);
    cout << "Read Files \n";
    
    //Find Edges and Normals
    MatrixXd S(3*C.rows(),6);
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
            k = k + 1;
        }
    }
    for (int i = 0; i < S.rows(); i++) {
        if (S(i,1) > 0) {
            if (S(i,5) == 0) {
                int t = int(S(i,1)+0.5);
                if (S(3*(t-1)+2,1) == S(i,0)) {
                    S(3*(t-1)+2,5) = i+1;
                } else if (S(3*(t-1)+1,1) == S(i,0)) {
                    S(3*(t-1)+1,5) = i+1;
                } else {
                    S(3*(t-1),5) = i+1;
                }
            }
        }
    }
    cout << "Built S \n";

    //Calculate Freestream State
    MatrixXd Uinf = FreeState();

    //Initialize Solver
    MatrixXd u(E.rows(),4);
    for (int i=0; i<u.rows(); i++) {
        u(i,all) = Uinf;
    }
 
    //Converge Solution
    double err = 100;
    int c = 0;
    double CFL = 1;
    cout << "Starting Loop \n";
    auto start = chrono::high_resolution_clock::now();
    while (err > 1e-5) {
        MatrixXd Rs = ResidCalc(u,S);
        MatrixXd R = Rs(all,seq(0,3));
        MatrixXd s = Rs(all,4);
        MatrixXd s4(s.rows(),4);
        s4 = s.replicate(1,4);
        u = u - 2*CFL*R.cwiseQuotient(s4);
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
    //FileWrite("CoarseMeshSol1.txt",u);
    return 0;
}

