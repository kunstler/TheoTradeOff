#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

int Torus(int i, int ii, int i_max);

int Reflexion(int i, int ii, int i_max);

double Prob_Law(double x, double y, double K);

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_sp, NumericVector c_,
               double K);



class CellAuto {
public:
    CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc,
             int nrow, int ncol,
	     Rcpp::NumericVector c_e, Rcpp::NumericVector c_, Rcpp::NumericVector c_s,
	     Rcpp::NumericVector ss,
	     double prob_distur, double prob_suc, double K);
    Rcpp::IntegerMatrix returnSp() const;
    Rcpp::IntegerMatrix returnSuc() const;
    void update();
    void updateCell(int ii, int jj);
    void ColonizeNeigh(int ii, int jj);
    void KillCellDisturb(int ii, int jj);
    void KillCellStress(int ii, int jj);
    void SucCell(int ii, int jj);
    void iterate(unsigned int iterations);
private:
    Rcpp::IntegerMatrix m_mat_sp;
    Rcpp::IntegerMatrix m_mat_suc;
    Rcpp:: NumericVector m_c_e;
    Rcpp:: NumericVector m_c_l;
    Rcpp:: NumericVector m_c_s;
    Rcpp:: NumericVector m_ss;
    int m_nrow;
    int m_ncol;
    double m_prob_distur;
    double m_prob_suc;
    double m_K;
};



int Torus(int i, int ii, int i_max){
 int i_n = ( i_max + ii + i ) % i_max;

 return i_n;
}

int Reflexion(int i, int ii, int i_max){
 int i_n = ii + i;
 if(i_n < 0) i_n = 1;
 if(i_n > i_max-1) i_n = i_max -2;

 return i_n;
}


CellAuto::CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc,
                   int nrow, int ncol,
		   Rcpp::NumericVector c_e, Rcpp::NumericVector c_l, Rcpp::NumericVector c_s,
		   Rcpp::NumericVector ss,
		   double prob_distur, double prob_suc, double K):
                                            m_mat_sp(Rcpp::clone(mat_sp)),
					    m_mat_suc(Rcpp::clone(mat_suc)),
					    m_c_e(Rcpp::clone(c_e)),
					    m_c_l(Rcpp::clone(c_l)),
					    m_c_s(Rcpp::clone(c_s)),
					    m_ss(Rcpp::clone(ss)),
                                            m_nrow(nrow),
                                            m_ncol(ncol),
					    m_prob_distur(prob_distur),
					    m_prob_suc(prob_suc),
                                            m_K(K)
{
}

Rcpp::IntegerMatrix CellAuto::returnSp() const{
  return m_mat_sp;
}

Rcpp::IntegerMatrix CellAuto::returnSuc() const{
  return m_mat_suc;
}



void CellAuto::iterate( unsigned int iterations ) {
  for ( int i = 0; i < iterations; i++ ) {
      update();
  }
}


void CellAuto::update() {
  for ( int r = 0; r < (m_nrow*m_ncol); r++ ) {
int i_r = floor(R::runif(0, m_nrow));
int j_r = floor(R::runif(0, m_ncol));
      CellAuto::updateCell(i_r , j_r);
  }
}


void CellAuto::updateCell(int ii, int jj) {
  // Rcpp::Rcout << "start colo" << std::endl;
  ColonizeNeigh(ii, jj);
  // Rcpp::Rcout << "start kill" << std::endl;
  KillCellDisturb(ii, jj);
  KillCellStress(ii, jj);
  SucCell(ii, jj);
}


void CellAuto::ColonizeNeigh(int ii, int jj) {
  int i, j, i_n, j_n, res;
int neig_r = floor(R::runif(0, 8));
  switch ( neig_r )
     {
        case 0:
           i = -1;
           j = -1;
           break;
        case 1:
           i = -1;
           j = 0;
           break;
        case 2:
           i = -1;
           j = 1;
           break;
        case 3:
           i = 0;
           j = -1;
           break;
        case 4:
           i = 0;
           j = 1;
           break;
        case 5:
           i = 1;
           j = -1;
           break;
        case 6:
           i = 1;
           j = 0;
           break;
        case 7:
           i = 1;
           j = 1;
           break;
     }
//select neighborhood cell to colonize
  i_n = Torus(i, ii, m_nrow);
  j_n = Reflexion(j, jj, m_ncol);
//competition interaction depending on the successional status
  if(m_mat_suc(i_n, j_n) == 1){
  	m_mat_sp(i_n, j_n) = Sampel_Law(i_n, j_n, ii, jj, m_mat_sp, m_c_e, m_K);
  }else{
  	m_mat_sp(i_n, j_n) = Sampel_Law(i_n, j_n, ii, jj, m_mat_sp, m_c_l, m_K);
  }
}


void CellAuto::KillCellDisturb(int ii, int jj){
  if(R::runif(0,1) < m_prob_distur){
	m_mat_sp(ii , jj) = 0;
	m_mat_suc(ii , jj) = 1;
  }
}


void CellAuto::KillCellStress(int ii, int jj){
    bool test;
    int res;
    if(m_mat_sp(ii, jj) != 0){
	if(R::runif(0,1) < (m_ss[jj] * m_c_s[m_mat_sp(ii, jj)])){
	    m_mat_sp(ii , jj) = 0;
	    m_mat_suc(ii , jj) = 1;
	}
    }
}

void CellAuto::SucCell(int ii, int jj){
  if(m_mat_sp(ii, jj) != 0){
     if(R::runif(0,1) < m_prob_suc){
        m_mat_suc(ii , jj) = 2;
     }
  }
}

double Prob_Law(double x, double y, double K){

  return  1 / (1 + exp(-K*(x - y)));
}

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_sp, NumericVector c_,
               double K){
  bool test;
  // Rcpp::Rcout << "ni and nj " << ni << " , " << nj << std::endl;
  // Rcpp::Rcout << "val sp " << mat_sp(ni, nj) << std::endl;
  // Rcpp::Rcout << "val c " << c_[mat_sp(ni, nj)]  << " and " << c_[mat_sp(ii, jj)]<< std::endl;
  // Rcpp::Rcout << "p law " << Prob_Law(c_[mat_sp(ni, nj)],c_[mat_sp(ii, jj)], K) << std::endl;
  int res;
  double pval = Prob_Law(c_[mat_sp(ni, nj)],c_[mat_sp(ii, jj)],K);
  if(c_[mat_sp(ii, jj)] != -100){
      if(c_[mat_sp(ni, nj)] != -100 ){
          if(R::runif(0,1) < pval){
	      res = mat_sp(ni, nj);
	  } else {
              res = mat_sp(ii, jj);
	  }
        } else {
          res = mat_sp(ii, jj);
        }
  }else{
      if(c_[mat_sp(ni, nj)] != -100 ){
          res = mat_sp(ni, nj);
        }else{
          res = mat_sp(ii, jj);
        }
  }
  // Rcpp::Rcout << "val sp col " <<  res << std::endl;

  return res;
}


//' Run cellular automaton based on initila condition and species parameters n times
//' 
//' @param mat_sp landscape matrix of cells with species code if occupied.
//' @param mat_succ landscape matrix of cells successional status.
//' @param c_e vector of species early successional competitive ability.
//' @param c_l vector of species late successional competitive ability.
//' @param c_s vector of species abiotic stress tolerance.
//' @param ss vector of abiotic stress.
//' @param prob_distur probablity of cell disturbance (mortality).
//' @param prob_suc probablity of cell succession.
//' @param K strength of the competitive asymmetry parameter.
//' @param n number of repetition of landscape updating.
//' @export
// [[Rcpp::export]]
Rcpp::List UpdateIterR(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
		       NumericVector c_e, NumericVector c_l,NumericVector c_s,
                       NumericVector ss,
		       double prob_distur, double prob_suc, double K, int n){

    int nrow = mat_sp.nrow();
    int ncol = mat_sp.ncol();

    CellAuto cells(mat_sp, mat_suc,
                   nrow, ncol,
                   c_e, c_l, c_s, ss,
		   prob_distur, prob_suc, K);
    cells.iterate(n);
    return Rcpp::List::create(Rcpp::Named("sp") = cells.returnSp(),
                              Rcpp::Named("suc") = cells.returnSuc());
}

