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

#include "ProbDistFunc.h"
#include <cmath>

double ABCFit::ProbDistFunc::Derivative(double x) {
  if (m_stepsize<=0.0) return std::nan("");
  return (Value(x-2*m_stepsize)-8*Value(x-m_stepsize)+8*Value(x+m_stepsize)-Value(x+2*m_stepsize))/12.0/m_stepsize; //Error is of order m_stepsize^4
}

double ABCFit::ProbDistFunc::SecondDerivative(double x) {
  if (m_stepsize<=0.0) return std::nan("");
  return (-Value(x-2*m_stepsize)+16*Value(x-m_stepsize)-30*Value(x)+16*Value(x+m_stepsize)-Value(x+2*m_stepsize))/12.0/m_stepsize/m_stepsize;//Error is of order m_stepsize^4
}

double ABCFit::ProbDistFunc::DerivativeLog(double x) {
  return -2.0*Derivative(x)/Value(x);
}

double ABCFit::ProbDistFunc::SecondDerivativeLog(double x) {
  return DerivativeLog(x)*DerivativeLog(x)/2.0 -2.0*SecondDerivative(x)/Value(x);
}
