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

#ifndef _COORDREPR_
#define _COORDREPR_
#include "MatrixAlgebra.h"

#include <string>
#include <vector>

namespace ABCFit {
const int m_NumberOfInternalCoordinates =
    4; // Internal representation is currently 3-momentum + mass
const double MassTolerance =
    1e-6; // Masses smaller than mass tolerance are set to be zero

class CoordRepr {
public:
  typedef Coordinates (*GetParFn)(Coordinates p);
  typedef Matrix (*GetParCovFn)(Coordinates p);
  virtual Coordinates
  Transform(Coordinates p) = 0; // translate from one represenation to internal
                                // representation
  virtual Coordinates
  InvTransform(Coordinates p) = 0; // translate from internal coordinate
                                   // representation to this represenation
  //Translate from this representation to output representation
  Coordinates Transform(Coordinates p, CoordRepr *output_CoordType) {
    return output_CoordType->InvTransform(this->Transform(p));
  }
  virtual Matrix
  Derivative(Coordinates p) = 0; // Derivative of internal representation w.r.t.
                                 // this representation
  virtual Matrix
  InvDerivative(Coordinates p) = 0; // Derivative of this representation w.r.t.
                                    // internal representation
  //Derivate of output representation w.r.t. this representation
  //d(output_CoordType)/d(this) = d(output_CoordType)/d(internal) * d(internal)/d(this)
  Matrix Derivative(Coordinates p, CoordRepr *output_CoordType) {
    return MatrixAlgebra::MatrixMultiplication(output_CoordType->InvDerivative(p),this->Derivative(p));
  }
  virtual std::string name() const = 0;
  virtual unsigned int GetNumberOfParameters() const {
    return 4;
  }; // Number of free parameters
  virtual unsigned int
  GetNumberOfCoordinates() const = 0; // Total number of coordinates (including
                                      // auxiliary coordinates)
  void SetParametrisationFunction(GetParFn GetPar) { m_GetPar = GetPar; };
  void SetParametrisationCovMatrixFunction(GetParCovFn GetParCov) {
    m_GetParCov = GetParCov;
  };
  Coordinates GetExpectation(Coordinates p) {
    if (m_GetPar != NULL)
      return m_GetPar(p);
    return GetDefaultExpectation(p);
  };
  Matrix GetExpectationCovMatrix(Coordinates p) {
    if (m_GetParCov != NULL)
      return m_GetParCov(p);
    return GetDefaultExpectationCovMatrix(p);
  };

private:
  virtual Coordinates GetDefaultExpectation(Coordinates p) = 0;
  virtual Matrix GetDefaultExpectationCovMatrix(Coordinates p) = 0;
  GetParFn m_GetPar;       // Function pointer
  GetParCovFn m_GetParCov; // Function pointer allows to change Parametrisation
                           // without inheritance
};
} // namespace ABCFit
#endif
