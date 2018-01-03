#ifndef MATRIX_H
#define MATRIX_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

///Helper class for matrix operations
struct Matrix
{
  ///Calculate the determinant
  //Adapted from the code Maik Beckmann posted at
  //  http://boost.2283326.n4.nabble.com/How-to-compute-determinant-td2710896.html
  static double CalcDeterminant(boost::numeric::ublas::matrix<double> m);

  ///Chop returns a std::vector of sub-matrices
  ///[ A at [0]   B at [1] ]
  ///[ C at [2]   D at [4] ]
  static const std::vector<boost::numeric::ublas::matrix<double> > Chop(
    const boost::numeric::ublas::matrix<double>& m);

  ///Create a n_rows x n_cols sized matrix from a std::vector,
  ///used for easy initialization
  static const boost::numeric::ublas::matrix<double> CreateMatrix(
    const std::size_t n_rows,
    const std::size_t n_cols,
    const std::vector<double>& v);

  ///Create a random-filled matrix
  static const boost::numeric::ublas::matrix<double> CreateRandomMatrix(
    const std::size_t n_rows, const std::size_t n_cols);

  ///Create a uBLAS vector from a std::vector,
  ///used for easy initialization
  static const boost::numeric::ublas::vector<double> CreateVector(const std::vector<double>& v);

  ///Calculate the inverse of a matrix
  static const boost::numeric::ublas::matrix<double> Inverse(
    const boost::numeric::ublas::matrix<double>& m);

  ///Check if two doubles are about equal
  static bool IsAboutEqual(const double a, const double b);

  ///Test these functions
  static void Test();

  ///Unchop merges the 4 std::vector of sub-matrices produced by Chop
  static const boost::numeric::ublas::matrix<double> Unchop(
    const std::vector<boost::numeric::ublas::matrix<double> >& v);

};

#endif // MATRIX_H
