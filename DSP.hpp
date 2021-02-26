class DSP
{
public:
  DSP(){
    std::cout << "C" << std::endl;
  }

  // data members
  double A[3][3];
  double D_[3][3];
  double Q_[3][3];
  double BinvXt_[3];
  int biasTowards;
  double deltaB_;

  // function members
  void diagonalize(const double (&A)[3][3],double (&Q)[3][3],double (&D)[3][3]);
  void sort(double (&Q)[3][3], double (&D)[3][3]);
  void perturb(double (&Q)[3][3], double (&D)[3][3]);  
  void form_perturbed_stress( const double (&D)[3][3], const double (&Q)[3][3], double (&A)[3][3]);
  void matrix_matrix_multiply( const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3]);

  ~DSP(){
    std::cout << "D" << std::endl;
  }
};
