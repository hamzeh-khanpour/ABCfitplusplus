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

#include "BreitWignerPDF.h"

namespace ABCFit {
  BreitWignerPDF::BreitWignerPDF(double mass, double width, double widthpower) : m_mass(mass), m_width(width), m_widthpower(widthpower) {}
  
  double BreitWignerPDF::Value(double x) {
    double temp = ConstrainValue(x);
    double gam = m_width;
    if (m_widthpower!=0.0) gam = gam*pow(temp/m_mass,m_widthpower);
    gam = gam*m_mass;
    temp = temp*temp - m_mass*m_mass;
    return m_norm*gam/(temp*temp+gam*gam); //relativistic BW)
  }

  double BreitWignerPDF::Derivative(double x) {
    double temp = ConstraintValue(x);
    double gam = m_width;
    if (m_widthpower!=0.0) gam = gam*pow(temp/m_mass,m_widthpower);
    double denom = temp*temp - m_mass*m_mass;
    denom = denom*denom+m_mass*m_mass*gam*gam;
    return ConstraintDerivative(x) * m_norm*gam*m_mass*(m_widthpower*(denom-2.0*m_mass*m_mass*gam*gam)/temp + 4.0*(m_mass*m_mass-temp*temp)*temp)/denom/denom; // first term is from chain-rule!
  }

  double BreitWignerPDF::SecondDerivative(double x) {
    double temp = ConstraintValue(x);
    double gam = m_width;
    if (m_widthpower!=0.0) gam = gam*pow(temp/m_mass,m_widthpower);
    double denom = temp*temp - m_mass*m_mass;
    denom = denom*denom+m_mass*m_mass*gam*gam;

    //bla / denom

    
    return 2.0*tan(x)*Derivative(x);
  }

  //
  
}
