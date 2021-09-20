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

#include "MatrixAlgebra.h"
#include <cmath>
#include <algorithm>
#include <iostream>
namespace ABCFit{

  Coordinates MatrixAlgebra::ScaleVector(Coordinates vector, double scale) {
    Coordinates result;
    result.reserve(vector.size());
    for (unsigned int i=0; i<vector.size(); i++) 
      result.push_back(vector.at(i)*scale);
    return result;
  }

  Matrix MatrixAlgebra::ScaleMatrix(Matrix matrix, double scale) {
    Matrix result;
    result.reserve(matrix.size());
    for (unsigned int i=0;i<matrix.size();i++)
      result.push_back(ScaleVector(matrix[i],scale));
    return result;
  }

  Coordinates MatrixAlgebra::VectorSum(Coordinates vector1, Coordinates vector2, double scale) {
    Coordinates result;
    if (vector1.size() != vector2.size()) return result;
    result.reserve(vector1.size());
    for (unsigned int i=0; i<vector1.size(); i++) 
      result.push_back(vector1.at(i)+scale*vector2.at(i));
    return result;
  }

  Matrix MatrixAlgebra::MatrixSum(Matrix matrix1, Matrix matrix2, double scale) {
    Matrix result;
    if (matrix1.size() != matrix2.size()) return result;
    if (matrix1.at(0).size() != matrix2.at(0).size()) return result;
    result.reserve(matrix1.size());
    for (unsigned int i=0; i<matrix1.size(); i++)
      result.push_back(VectorSum(matrix1.at(i),matrix2.at(i),scale));
    return result;
  }
  
  double MatrixAlgebra::VectorMultiplication(Coordinates vector1, Coordinates vector2){
    double result = 0.0;
    if (vector1.size() != vector2.size()) return std::nan("");
    for (unsigned int i=0; i<vector1.size(); i++)
      result+=vector1.at(i)*vector2.at(i);
    return result;
  }
  
  Coordinates MatrixAlgebra::MatrixVectorMultiplication(Matrix matrix, Coordinates vector){
    if (matrix.at(0).size() != vector.size()) return Coordinates();
    Coordinates result(matrix.size(),0.0);
    
    for (unsigned int i=0; i<matrix.size(); i++)
      result[i]=VectorMultiplication(matrix.at(i),vector);
    return result;
  }

  Coordinates MatrixAlgebra::VectorMatrixMultiplication(Coordinates vector, Matrix matrix){
    if (matrix.size() != vector.size()) return Coordinates();
    Coordinates result;
    result.reserve(matrix.at(0).size());

    for (unsigned int i=0; i<matrix.at(0).size(); i++) {
      Coordinates col;
      col.reserve(matrix.size());
      for (unsigned int j=0; j<matrix.size(); j++)
        col.push_back(matrix[j][i]);
      result.push_back(VectorMultiplication(vector,col));

    }
    return result;
  }
  
  Matrix MatrixAlgebra::MatrixMultiplication(Matrix matrix1, Matrix matrix2){
    if (matrix1.at(0).size() != matrix2.size()) return Matrix();
    Matrix result(matrix1.size(),Coordinates());
    
    for (unsigned int i=0; i<matrix2.at(0).size(); i++) {
      Coordinates col;
      for (unsigned int j=0; j<matrix2.size(); j++)
	col.push_back(matrix2[j][i]);
      col=MatrixVectorMultiplication(matrix1,col);
      for (unsigned int j=0; j<col.size(); j++) result[j].push_back(col[j]);
    }
    return result;
  }
  
  Matrix MatrixAlgebra::Cofactor(Matrix* matrix, unsigned int row, unsigned int col) {
    unsigned int size = matrix->size();
    Matrix result;
    result.reserve(size);
    if (size < row || matrix->at(0).size() < col) return result;
    for (unsigned int i = (row == 0 ? 1:0); i < size; (row == ++i ? i++:i)) {
      Coordinates coordinates;
      coordinates.reserve(size);
      for (unsigned int j = (col == 0 ? 1:0); j < matrix->at(0).size(); (col == ++j ? j++:j))
	coordinates.push_back((*matrix)[i][j]); 
      result.push_back(coordinates);
    }
    return result;
  }
  
  double MatrixAlgebra::Determinant(Matrix* matrix) {
    unsigned int size = matrix->size();
    if (size != matrix->at(0).size()) return std::nan("");
    if (size == 1) return (*matrix)[0][0];

    double result = 0.0;
    for (unsigned int i = 0; i < size; i++) {
      Matrix coMatrix = Cofactor(matrix, 0, i);
      result += (((i)%2==0)? 1.0: -1.0) * (*matrix)[0][i] * Determinant(&coMatrix);
    }
    return result;
  }
  
  Matrix MatrixAlgebra::Adjoint(Matrix* matrix) {
    unsigned int size = matrix->size();
    Matrix result;
    result.reserve(size);
    if (size != matrix->at(0).size()) return result;
    if (size == 1) return Matrix(1,Coordinates(1,1.0));
    for (unsigned int i = 0; i < size; i++) {
      Coordinates coordinates;
      coordinates.reserve(size);
      for (unsigned int j = 0; j < size; j++) {
        Matrix coMatrix = Cofactor(matrix, j, i);
	coordinates.push_back((((i+j)%2==0)? 1.0: -1.0)*Determinant(&coMatrix));
      }
      result.push_back(coordinates);
    }
    return result;
  }
  
  Matrix MatrixAlgebra::InverseMatrix(Matrix* matrix) {
    double det = Determinant(matrix);
    if (!std::isnormal(det)) return {{std::nan("")}};
    
    Matrix adj = Adjoint(matrix);
    for (unsigned int i = 0; i < matrix->size(); i++) {
      std::transform(adj[i].begin(), adj[i].end(), adj[i].begin(), [&det](double& c){return c/det;});
      //result.push_back(std::transform(row.begin(), row.end(), row.begin(), [&det](double& c){return c/det;}));
    }
    return adj;
  }

  Matrix MatrixAlgebra::TransposeMatrix(Matrix matrix) {
    unsigned int rowSize = matrix.size();
    unsigned int colSize = matrix[0].size();
    if (rowSize == 0 || colSize == 0) return matrix;
    Matrix result(colSize, Coordinates (rowSize, 0.0));
    for (unsigned int row=0; row<rowSize; row++)
      for (unsigned int col=0; col<colSize; col++)
	result[col][row]=matrix[row][col];
    return result;
  }
  
  bool MatrixAlgebra::CheckRow(Coordinates row) {
    for (auto x: row) 
      if (!std::isfinite(x)) return false;
    return row.size()>0;
  }

  bool MatrixAlgebra::CheckMatrix(Matrix matrix) {
    for (auto row: matrix)
      if (!CheckRow(row)) return false;
    return matrix.size()>0;
  }

}
