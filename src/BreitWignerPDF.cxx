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
    return m_norm*gam/(temp*temp+gam*gam); //relativistic BW
  }

  double BreitWignerPDF::Derivative(double x) {
    double temp = ConstraintValue(x);
    double gam = m_width;
    if (m_widthpower!=0.0) gam = gam*pow(temp/m_mass,m_widthpower);
    gam = gam*m_mass;
    double denom = temp*temp - m_mass*m_mass;
    double factor = 2.0*(m_widthpower*gam*gam+2.0*temp*temp*denom);
    denom = denom*denom+gam*gam;
    return ConstraintDerivative(x) * m_norm*gam*(m_widthpower - factor/denom)/temp/denom; // first term is from chain-rule!
  }

  double BreitWignerPDF::SecondDerivative(double x) {
    double mass = ConstraintValue(x); //m
    double gam = m_width;
    if (m_widthpower!=0.0) gam = gam*pow(mass/m_mass,m_widthpower);
    gam = gam*m_mass;
    double temp = mass*mass; //m^2
    double denom = temp - m_mass*m_mass;
    double factor = 2.0*(m_widthpower*gam*gam+2.0*temp*denom);
    denom = denom*denom+gam*gam;
    /*
To test derivative using sympy@python:
----------------------------
from sympy import *
N, m, w, x, b = symbols('N m w x b')
init_printing()
bw = N*m*w*(x/m)**b/((x**2-m**2)**2+m**2*w**2*(x/m)**(2*b))
d2bwdx2 = diff(bw, x, 2)
temp = x**2
factor = 2*(b*m**2*w**2*(x/m)**(2*b)+2*x**2*(x**2-m**2))
denom = (x**2-m**2)**2+m**2*w**2*(x/m)**(2*b)
code = N*m*w*(x/m)**b*((factor*(2*factor/denom-4*b+1)+8*temp*((b-1)*(temp-m**2)-temp))/denom+b*(b-1))/temp/denom
testd2 = d2bwdx2*temp*denom**3/N/m/w/(x/m)**b
testcode = code*temp*denom**3/N/m/w/(x/m)**b
difference = testcode-testd2
difference.simplify()
-----------------------------
gives: 0
    */
    // using chain rule on BW(tan(param)+m_mass): d2BW/dparam2 = d(dBW/dM * dM/dparam) = d2BW/dM2 * (dM/dparam)^2 + dBW/dM * d2M/dparam2
    double dMdparam = ConstraintDerivative(x);
    return m_norm*gam*dMdparam*(
				dMdparam*((factor*(2.0*factor/denom - 4.0*m_widthpower + 1.0) + 8.0*temp*((m_widthpower-1.0)*(temp-m_mass*m_mass)-temp))/denom
					  + m_widthpower*(m_widthpower-1.0))/mass
				+ 2.0*tan(x)*(m_widthpower - factor/denom)    // d2M/dparam2 = 2tan*(tan^2 + 1) = 2tan*dM/dparam
				)/mass/denom;
  }
  
}
