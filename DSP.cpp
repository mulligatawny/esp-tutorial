/*---------------------------------------------------------------------------*/
/* Testbed for perturbation of subgrid stress anisotropy towards resolved    */
/* stress anisotropy                                                         */
/*---------------------------------------------------------------------------*/

#include<iostream>
#include<cmath>
#include<iomanip>

void print(double (&A)[3][3]);
void form_anisotensor(double (&A)[3][3], double (&B)[3][3]);

void diagonalize(
  const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]);

void sort(
  double (&Q)[3][3], double (&D)[3][3]);

void perturb(
  double (&Q)[3][3], double (&D)[3][3], double (&P)[3][3]);

void form_perturbed_stress(
  const double (&D)[3][3], const double (&Q)[3][3], double (&A)[3][3]);

void matrix_matrix_multiply(
  const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3]);

int main() {

  // define and declare sgs and resolved tensors
  double S_sgs[3][3] = { {2.0,-1.0,0.0}, 
                         {-1.0,2.0,-1.0}, 
                         {0.0,-1.0,2.0} }; 
  double S_res[3][3] = { {24.0, 14.0, 4.0}, 
                         {14.0, 13.0, 17.0}, 
                         {4.0, 17.0, 51.0} }; 
  double A_sgs[3][3];
  double A_res[3][3];

  // form sgs and resolved anisotropy tensors
  form_anisotensor(S_sgs, A_sgs);
  form_anisotensor(S_res, A_res);

  // print to check
  /*
  std::cout << "The resolved anisotropy tensor is:" << std::endl;
  print(A_res);
  */
  std::cout << "The sgs anisotropy tensor is:" << std::endl;
  print(A_sgs);

  // declare some matrices for the spectral decomposition
  double Q_sgs[3][3];
  double Q_res[3][3];
  double D_sgs[3][3];
  double D_res[3][3];

  // get resolved stress eigenvalues
  diagonalize(A_res, Q_res, D_res);
  sort(Q_res, D_res);
  /*
  std::cout << "The resolved eigenvalues are:" << std::endl;
  print(D_res);
  std::cout << "The resolved eigenvectors are:" << std::endl;
  print(Q_res);
  */
  // get sgs stress eigenvalues
  diagonalize(A_sgs, Q_sgs, D_sgs);
  sort(Q_sgs, D_sgs);
  /*
  std::cout << "The sgs eigenvalues are:" << std::endl;
  print(D_sgs);
  std::cout << "The sgs eigenvectors are:" << std::endl;
  print(Q_sgs);
  */
  perturb(Q_sgs, D_sgs, D_res);
  form_perturbed_stress(D_sgs, Q_sgs, A_sgs);

  std::cout << "The perturbed sgs anisotropy tensor is:" << std::endl;
  print(A_sgs);  

  return 0;
}

void print(double (&A)[3][3]) {
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) { 
      std::cout << std::setprecision(8) << A[i][j] << "\t"; 
    }
    std::cout << std::endl; 
  }
}


void form_anisotensor(double (&A)[3][3], double (&B)[3][3]) {
  double tr = A[0][0] + A[1][1] + A[2][2];
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      if (i==j) {
        B[i][j] = A[i][j]/tr -1.0/3.0;
      }
      else {
        B[i][j] = A[i][j]/tr;
      }
    }
  }
}

void diagonalize(
  const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3])
{
  /*
    obtained from:
    http://stackoverflow.com/questions/4372224/
    fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
    A must be a symmetric matrix.
    returns Q and D such that
    Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
  */

  const int maxsteps=24;
  int k0, k1, k2;
  double o[3], m[3];
  double q [4] = {0.0,0.0,0.0,1.0};
  double jr[4];
  double sqw, sqx, sqy, sqz;
  double tmp1, tmp2, mq;
  double AQ[3][3];
  double thet, sgn, t, c;
  for(int i=0;i < maxsteps;++i) {
    // quat to matrix
    sqx      = q[0]*q[0];
    sqy      = q[1]*q[1];
    sqz      = q[2]*q[2];
    sqw      = q[3]*q[3];
    Q[0][0]  = ( sqx - sqy - sqz + sqw);
    Q[1][1]  = (-sqx + sqy - sqz + sqw);
    Q[2][2]  = (-sqx - sqy + sqz + sqw);
    tmp1     = q[0]*q[1];
    tmp2     = q[2]*q[3];
    Q[1][0]  = 2.0 * (tmp1 + tmp2);
    Q[0][1]  = 2.0 * (tmp1 - tmp2);
    tmp1     = q[0]*q[2];
    tmp2     = q[1]*q[3];
    Q[2][0]  = 2.0 * (tmp1 - tmp2);
    Q[0][2]  = 2.0 * (tmp1 + tmp2);
    tmp1     = q[1]*q[2];
    tmp2     = q[0]*q[3];
    Q[2][1]  = 2.0 * (tmp1 + tmp2);
    Q[1][2]  = 2.0 * (tmp1 - tmp2);

    // AQ = A * Q
    AQ[0][0] = Q[0][0]*A[0][0]+Q[1][0]*A[0][1]+Q[2][0]*A[0][2];
    AQ[0][1] = Q[0][1]*A[0][0]+Q[1][1]*A[0][1]+Q[2][1]*A[0][2];
    AQ[0][2] = Q[0][2]*A[0][0]+Q[1][2]*A[0][1]+Q[2][2]*A[0][2];
    AQ[1][0] = Q[0][0]*A[0][1]+Q[1][0]*A[1][1]+Q[2][0]*A[1][2];
    AQ[1][1] = Q[0][1]*A[0][1]+Q[1][1]*A[1][1]+Q[2][1]*A[1][2];
    AQ[1][2] = Q[0][2]*A[0][1]+Q[1][2]*A[1][1]+Q[2][2]*A[1][2];
    AQ[2][0] = Q[0][0]*A[0][2]+Q[1][0]*A[1][2]+Q[2][0]*A[2][2];
    AQ[2][1] = Q[0][1]*A[0][2]+Q[1][1]*A[1][2]+Q[2][1]*A[2][2];
    AQ[2][2] = Q[0][2]*A[0][2]+Q[1][2]*A[1][2]+Q[2][2]*A[2][2];
    // D = Qt * AQ
    D[0][0] = AQ[0][0]*Q[0][0]+AQ[1][0]*Q[1][0]+AQ[2][0]*Q[2][0];
    D[0][1] = AQ[0][0]*Q[0][1]+AQ[1][0]*Q[1][1]+AQ[2][0]*Q[2][1];
    D[0][2] = AQ[0][0]*Q[0][2]+AQ[1][0]*Q[1][2]+AQ[2][0]*Q[2][2];
    D[1][0] = AQ[0][1]*Q[0][0]+AQ[1][1]*Q[1][0]+AQ[2][1]*Q[2][0];
    D[1][1] = AQ[0][1]*Q[0][1]+AQ[1][1]*Q[1][1]+AQ[2][1]*Q[2][1];
    D[1][2] = AQ[0][1]*Q[0][2]+AQ[1][1]*Q[1][2]+AQ[2][1]*Q[2][2];
    D[2][0] = AQ[0][2]*Q[0][0]+AQ[1][2]*Q[1][0]+AQ[2][2]*Q[2][0];
    D[2][1] = AQ[0][2]*Q[0][1]+AQ[1][2]*Q[1][1]+AQ[2][2]*Q[2][1];
    D[2][2] = AQ[0][2]*Q[0][2]+AQ[1][2]*Q[1][2]+AQ[2][2]*Q[2][2];
    o[0]    = D[1][2];
    o[1]    = D[0][2];
    o[2]    = D[0][1];
    m[0]    = std::abs(o[0]);
    m[1]    = std::abs(o[1]);
    m[2]    = std::abs(o[2]);
    // index of largest element of offdiag
    k0      = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; 
    k1      = (k0+1)%3;
    k2      = (k0+2)%3;
    if (o[k0]==0.0) {
      break;  // diagonal already
    }
    thet    = (D[k2][k2]-D[k1][k1])/(2.0*o[k0]);
    sgn     = (thet > 0.0)?1.0:-1.0;
    thet   *= sgn; // make it positive
    // sign(T)/(|T|+sqrt(T^2+1))
    t       = sgn /(thet +((thet < 1.E6)? std::sqrt(thet*thet+1.0):thet)) ; 
    c       = 1.0/std::sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c
    if(c==1.0) {
      break;  // no room for improvement - reached machine precision.
    }
    jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
    // using 1/2 angle identity sin(a/2) = std::sqrt((1-cos(a))/2)
    jr[k0]  = sgn*std::sqrt((1.0-c)/2.0);  
    // since our quat-to-matrix convention was for v*M instead of M*v
    jr[k0] *= -1.0;
    jr[3 ]  = std::sqrt(1.0f - jr[k0] * jr[k0]);
    if(jr[3]==1.0) {
      break; // reached limits of floating point precision
    }
    q[0]    = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
    q[1]    = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
    q[2]    = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
    q[3]    = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
    mq      = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0]   /= mq;
    q[1]   /= mq;
    q[2]   /= mq;
    q[3]   /= mq;
  }
}

void
sort(
  double (&Q)[3][3], double (&D)[3][3])
{
  // Goal: sort diagonalization eigenvalues from high to low;
  // save off row in D from high to low; reorder Q tensor accordingly

  for(int i = 0; i < 2; ++i) {

    if(D[0][0] < D[1][1]) {

      double tempEigenvalue = D[0][0];
      D[0][0] = D[1][1];
      D[1][1] = tempEigenvalue;

      double tempEigenvector0 = Q[0][0]; 
      double tempEigenvector1 = Q[1][0]; 
      double tempEigenvector2 = Q[2][0];
      Q[0][0] = Q[0][1]; Q[1][0] = Q[1][1]; Q[2][0] = Q[2][1];
      Q[0][1] = tempEigenvector0; 
      Q[1][1] = tempEigenvector1; 
      Q[2][1] = tempEigenvector2;

    }

    if(D[1][1] < D[2][2]) {

      double tempEigenvalue = D[1][1];
      D[1][1] = D[2][2];
      D[2][2] = tempEigenvalue;

      double tempEigenvector0 = Q[0][1]; 
      double tempEigenvector1 = Q[1][1]; 
      double tempEigenvector2 = Q[2][1];
      Q[0][1] = Q[0][2]; Q[1][1] = Q[1][2]; Q[2][1] = Q[2][2];
      Q[0][2] = tempEigenvector0; 
      Q[1][2] = tempEigenvector1; 
      Q[2][2] = tempEigenvector2;

    }
  }
}

void
perturb(
  double (&Q)[3][3], double (&D)[3][3], double(&P)[3][3])
{
  // sgs stress eigenvalue perturbation

  // extract sorted by size (L1 > L2 > L3)
  const double Lambda1 = D[0][0];
  const double Lambda2 = D[1][1];
  const double Lambda3 = D[2][2];

  // extract resolved stress eigenvalues
  const double L_res1 = P[0][0];
  const double L_res2 = P[1][1];
  const double L_res3 = P[2][2];

  const double deltaB_ = 0.789;

  // perturbed sgs eigenvalues
  const double pLambda1 = (1.0 - deltaB_)*Lambda1 + deltaB_*L_res1;
  const double pLambda2 = (1.0 - deltaB_)*Lambda2 + deltaB_*L_res2;
  const double pLambda3 = (1.0 - deltaB_)*Lambda3 + deltaB_*L_res3;

  D[0][0] = pLambda1;
  D[1][1] = pLambda2;
  D[2][2] = pLambda3;

}

void
form_perturbed_stress(
  const double (&D)[3][3], const double (&Q)[3][3], double (&A)[3][3])
{
  // A = Q*D*QT
  double QT[3][3];
  double B[3][3];

  // compute QT
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      QT[j][i] = Q[i][j];
    }
  }
  //mat-vec, B = Q*D
  matrix_matrix_multiply(Q,D,B);

  // mat-vec, A = (Q*D)*QT = B*QT
  matrix_matrix_multiply(B,QT,A);
}

void
matrix_matrix_multiply(
  const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3])
{
  //C = A*B
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0;
      for (int k = 0; k < 3; ++k) {
        sum = sum + A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}
