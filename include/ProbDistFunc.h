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

#ifndef _PROBDISTFUNC_
#define _PROBDISTFUNC_

namespace ABCFit {
// ProbDistFunc has to be greater than zero (as implementation uses the function
// g(x)=-2log(PDF(x))
class ProbDistFunc {
public:
  ProbDistFunc(double stepsize = -1.0) : m_stepsize(stepsize){};
  virtual double Value(double x) = 0;
  virtual double ExpectationValue() = 0; // Most probable value used for chi2 contribution
  virtual double ConstrainValue(double x) { return x; } // to allow coordinate transformation
  virtual double ConstrainDerivative(double x) { return 1.0; } // to allow coordinate transformation  
  virtual double Derivative(double x);
  double DerivativeLog(double x); // Derivative of g(x)=-2log(PDF(x))
  virtual double SecondDerivative(double x);
  double SecondDerivativeLog(double x); // Second derivative of g(x)=-2log(PDF(x))
private:
  double m_stepsize;
};
} // namespace ABCFit
#endif
