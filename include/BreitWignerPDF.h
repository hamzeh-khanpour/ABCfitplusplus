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

#ifndef _BREITWIGNERPDF_
#define _BREITWIGNERPDF_
#include "ProbDistFunc.h"
#include <cmath>

// Probaility distribution for a relativistic Breit-Wigner with optional running width
//
//       prob = norm*Gam*mref/((m^2-mref^2)^2 + Gam^2*mref^2)
//with
//       gam = width * (m/mref)^widthpower
//
// To improve fit stability around peak-value a coordinate-transformation is used:
//          Mass = tan(ParameterValue) + Peak_Mass

namespace ABCFit {

class BreitWignerPDF : virtual public ProbDistFunc {
public:
  BreitWignerPDF(double mass  = 91.183, double width = 2.5, double widthpower = 0.0);
  double Value(double x);
  double ExpectationValue() {return 0.0; } // Most probable parameter value used for chi2 contribution
  double ConstraintValue(double x) { return tan(x) + m_mass; } // Coordinate transformation
  double ConstraintDerivative(double x) { return tan(x)*tan(x) + 1.0; } // Derivative of coordinate transformation  
  double Derivative(double x);
  double SecondDerivative(double x);
  void SetParameters(double mass = 91.183, double width = 2.5, double widthpower = 0.0) {
    m_mass = mass;
    m_width = width;
    m_widthpower = widthpower;
  };

private:
  double m_mass;
  double m_width;
  double m_widthpower;
  const double m_norm = 1.0/3.14159265359; // 1/pi -- should include widthpower, but not relevant....
};
} // namespace ABCFit
#endif
