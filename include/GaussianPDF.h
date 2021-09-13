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

#ifndef _GAUSSIANPDF_
#define _GAUSSIANPDF_
#include "ProbDistFunc.h"
#include <cmath>

namespace ABCFit {

class GaussianPDF : virtual public ProbDistFunc {
public:
  GaussianPDF(double mean = 0.0, double sigma = 1.0);
  double Value(double x) {
    return (m_norm / m_sigma) *
           exp(-0.5 * (x - m_mean) * (x - m_mean) / (m_sigma * m_sigma));
  };
  double ExpectationValue() {
    return m_mean;
  }; // Most probable value used for chi2 contribution
  double Derivative(double x) {
    return -((x - m_mean) / (m_sigma * m_sigma)) * Value(x);
  };
  double SecondDerivative(double x) {
    return -Value(x) / (m_sigma * m_sigma) +
           ((x - m_mean) * (x - m_mean) /
            (m_sigma * m_sigma * m_sigma * m_sigma)) *
               Value(x);
  };
  void SetParameters(double mean = 0.0, double sigma = 1.0) {
    m_mean = mean;
    m_sigma = sigma;
  };

private:
  double m_mean;
  double m_sigma;
  const double m_norm = 0.3989422804; // 1/sqrt(2*pi)
};
} // namespace ABCFit
#endif
