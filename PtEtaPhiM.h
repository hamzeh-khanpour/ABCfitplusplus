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

#ifndef _PTETAPHIM_
#define _PTETAPHIM_
#include "CoordRepr.h"

namespace ABCFit {

class PtEtaPhiM : virtual public CoordRepr {
public:
  Coordinates Transform(Coordinates p); // translate from one represenation to
                                        // internal representation
  Coordinates
  InvTransform(Coordinates p); // translate from internal coordinate
                               // representation to this represenation
  std::vector<Coordinates>
  Derivative(Coordinates p); // Derivative of internal representation w.r.t.
                             // this representation
  std::vector<Coordinates>
  InvDerivative(Coordinates p); // Derivative of this representation w.r.t.
                                // internal representation
  std::string name() const { return "PtEtaPhiM"; };
  unsigned int GetNumberOfCoordinates() const {
    return 4;
  }; // Total number of coordinates (including auxiliary coordinates)
private:
  Coordinates GetDefaultExpectation(Coordinates p);
  Matrix GetDefaultExpectationCovMatrix(Coordinates p);
};
static CoordRepr *ptetaphim = new PtEtaPhiM();
} // namespace ABCFit
#endif
