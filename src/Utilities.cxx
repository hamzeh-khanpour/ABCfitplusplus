//Auxiliary functions useful for setting up constrained fit

#include "CoordRepr.h"
#include <cmath>

namespace ABCFit{

  //Rescale particles e.g. jets to have zero mass
  Coordinates RescaleParticle2ZeroMass(Coordinates particle){
    double scale = particle[0]/sqrt(particle[1]*particle[1]+particle[2]*particle[2]+particle[3]*particle[3]);
    return {particle[0], scale*particle[1], scale*particle[2], scale*particle[3]};
  }
  
  //Calculate parameters for a reconstructed particle in reference to the corresponding true particle (e.g. using matching, association etc.)
  //Distributions of the parameters are relevant for setting up the expectation values and covariance matrices
  //Step 1: Make distributions of calculated parameters. Can be split into particle/type. 
  //Step 2: Fit distributions by Gaussian (chi2 parameters must be Gaussian, non-Gaussian distribution must be suitably transformed)
  //Step 3: Set up expecation values/covariance matrix functions: 
  //        Expectation value template function: Coordinates TestParameterisation(Coordinates p)
  //        Covariance matrix template function: Matrix TestParameterisationCovMatrix(Coordinates p)
  //NOTE: paramtetrisation functions are unspecified user functions and can be made as detailed as required (e.g. energy/direction dependent)
  Coordinates CalculateParameters(Coordinates _reco, Coordinates _truth, CoordRepr* input_repr, CoordRepr* parametrisation_repr) {
    Coordinates result = {4,-10.0};
    Coordinates reco{_reco[1],_reco[2],_reco[3],_reco[0]};
    Coordinates truth{_truth[1],_truth[2],_truth[3],_truth[0]};
    Coordinates ref = parametrisation_repr->InvTransform(input_repr->Transform(reco)); //Setup reference 
    Coordinates MC  = input_repr->Transform(truth);
    for (unsigned int i=0; i<3; i++) {
      ref[i]=MC[i];
    }
    ref=parametrisation_repr->InvTransform(ref);
    for (unsigned int i=0; i<3; i++) {
      result[i]=ref[i];
    }
    result[3]=truth[3]-reco[3];
    return result;
  }
}
