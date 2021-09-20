//Auxiliary functions useful for setting up constrained fit

#include "CoordRepr.h"
#include "PxPyPzE.h"
#include "ABCD.h"
#include "ParticleObject.h"
#include <cmath>

namespace ABCFit{

  //Rescale particles e.g. jets to have zero mass where coodinates are given in PxPyPzE representation
  Coordinates RescaleParticle2ZeroMass(Coordinates particle){
    double scale = particle[3]/sqrt(particle[0]*particle[0]+particle[1]*particle[1]+particle[2]*particle[2]);
    return {scale*particle[0], scale*particle[1], scale*particle[2], particle[3]};
  }

  //Create particle object from px, py, pz and E parameters e.g. read from flat ntuple file produced in FCCAnalyses
  ParticleObject* CreateParticleObject(float px, float py, float pz, float e) {
    Coordinates particle = {px,py,pz,e};
    Coordinates temp = RescaleParticle2ZeroMass(particle);
    ParticleObject* result = new ParticleObject(temp, pxpypze, abcd);
    result->FixParameter(3); //zero mass should not be a parameter
    return result;
  }
  
  //Calculate parameters for a reconstructed particle in reference to the corresponding true particle (e.g. using matching, association etc.)
  //Here _reco and _truth are represented in input_repr
  //Distributions of the parameters are relevant for setting up the expectation values and covariance matrices
  //Step 1: Make distributions of calculated parameters. Can be split into particle/type. 
  //Step 2: Fit distributions by Gaussian (chi2 parameters must be Gaussian, non-Gaussian distribution must be suitably transformed)
  //Step 3: Set up expecation values/covariance matrix functions: 
  //        Expectation value template function: Coordinates TestParameterisation(Coordinates p)
  //        Covariance matrix template function: Matrix TestParameterisationCovMatrix(Coordinates p)
  //NOTE: paramtetrisation functions are unspecified user functions and can be made as detailed as required (e.g. energy/direction dependent)
  Coordinates CalculateParameters(Coordinates _reco, Coordinates _truth, CoordRepr* input_repr, CoordRepr* parametrisation_repr) {
    Coordinates result = {4,-10.0};
    Coordinates ref = parametrisation_repr->InvTransform(input_repr->Transform(_reco)); //Setup reference 
    Coordinates MC  = input_repr->Transform(_truth);
    for (unsigned int i=0; i<3; i++) {
      ref[i]=MC[i];
    }
    ref=parametrisation_repr->InvTransform(ref);
    for (unsigned int i=0; i<3; i++) {
      result[i]=ref[i];
    }
    result[3]=_truth[3]-_reco[3];
    return result;
  }
}
