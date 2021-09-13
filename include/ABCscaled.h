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

#ifndef _ABCSCALED_
#define _ABCSCALED_
#include "CoordRepr.h"

namespace ABCFit {

class ABCscaled : virtual public ABCD {
  public:
    std::string name() const { return "ABCscaled"; };
    unsigned int GetNumberOfParameters() const {
      return 3;
    }; // Number of free parameters
  private:
  };
  static CoordRepr *abcscaled = new ABCscaled();
} // namespace ABCFit
#endif
