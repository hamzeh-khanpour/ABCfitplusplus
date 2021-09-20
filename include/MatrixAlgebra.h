//  ABCfit++ -- General constrained kinematic fit using ABC-parametrisation
//  Copyright (C) 2021
//    - Jørgen Beck Hansen, <beck@nbi.ku.dk>
//    - Julie Munch Torndal, <julie@torndal.com>
//
//  This file is part of ABCfit++.
//
//  ABCfit++ is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ABCfit++ is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef _MATRIXALGEBRA_
#define _MATRIXALGEBRA_
//#include "CoordRepr.h"
#include <iostream>
#include <vector>

namespace ABCFit {
typedef std::vector<double>
    Coordinates; // First 4 coordinates are allocated for internal
                 // representation (4-momentum), auxiliary coordinates should lie
                 // afterwards
typedef std::vector<Coordinates> Matrix;

namespace MatrixAlgebra {
template <typename T> void print(std::vector<T> const &v, bool raw = false) {
  if (!raw)
    std::cout << " ( ";
  for (auto i : v) {
    std::cout << i << ' ';
  }
  if (!raw)
    std::cout << " )" << std::endl;
}
template <typename T>
void print(std::vector<std::vector<T>> const &v, bool size = false) {
  std::cout << " ( ";
  if (v.size() == 0)
    std::cout << "empty";
  for (unsigned int i = 0; i < v.size(); i++) {
    print(v[i], true);
    if (i < v.size() - 1)
      std::cout << std::endl;
  }
  std::cout << " )";
  if (size && v.size() != 0)
    std::cout << "_[" << v.size() << "," << v.at(0).size() << "]_";
  std::cout << std::endl;
}

template <typename T>
std::vector<T> SubVector(std::vector<T> const &v, int m, int n) {
  auto first = v.cbegin() + m;
  auto last = v.cbegin() + n + 1;

  std::vector<T> vec(first, last);
  return vec;
}

Coordinates ScaleVector(Coordinates vector, double scale = 1.0); // scale*vector
Matrix ScaleMatrix(Matrix matrix, double scale = 1.0);           // scale*matrix
Coordinates VectorSum(Coordinates vector1, Coordinates vector2,
                      double scale = 1.0); // vector1+scale*vector2
Matrix MatrixSum(Matrix matrix1, Matrix matrix2,
                 double scale = 1.0); // Matrix1 + scale*Matrix2
double VectorMultiplication(Coordinates vector1, Coordinates vector2);
Coordinates MatrixVectorMultiplication(Matrix matrix, Coordinates vector);
Coordinates VectorMatrixMultiplication(Coordinates vector, Matrix matrix);
Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2);
Matrix Cofactor(Matrix *matrix, unsigned int row, unsigned int col);
double Determinant(Matrix *matrix);
Matrix Adjoint(Matrix *matrix);
Matrix InverseMatrix(Matrix *matrix);
Matrix TransposeMatrix(Matrix matrix);
bool CheckRow(Coordinates row);
bool CheckMatrix(Matrix matrix);
} // namespace MatrixAlgebra
} // namespace ABCFit
#endif
