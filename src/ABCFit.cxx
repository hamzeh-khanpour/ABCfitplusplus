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

#include "ABCFit.h"
#include "MatrixAlgebra.h"
#include "ABCD.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>

#define HACK
namespace ABCFit{
  FitResult::FitResult (std::vector<ParticleObject> _FittedParticles, double _chi2, int _Ndof, statuscode _status, unsigned int _Niter) : FittedParticles(_FittedParticles), chi2(_chi2), Ndof(_Ndof), status(_status), Niter(_Niter) {};

  void FitResult::print(bool PrintParticles) {
    std::cout << "Status = " << status << std::endl;
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << "Ndof = " << Ndof << std::endl;
    std::cout << "Niter = " << Niter << std::endl;
    if (PrintParticles) {
      for (unsigned int i=0; i<FittedParticles.size(); i++) {
	Coordinates coor = FittedParticles[i].GetCoordinates();
	std::cout << "Particle(" << i << "): " << coor[0] << "   " << coor[1]<<"   " << coor[2]<<"   " << coor[3]<< std::endl;
      }
    }
  }

  FitKernel::FitKernel(std::vector<CompositeConstraint> ListOfConstraints, unsigned int maxiter) : m_ListOfConstraints(ListOfConstraints), m_MaxIterations(maxiter) {
    if (!Initialize()) m_State = statuscode::InitiateFailure;
    else m_State = statuscode::Initiate;    
  };

  FitResult FitKernel::Fit(unsigned int maxiter) {
    if (maxiter>0) m_MaxIterations=maxiter;
    FitResult result(std::vector<ParticleObject> {}, -99.0, -99, statuscode::Initiate, 0);
    Coordinates ExpectationValues;
    Coordinates UnmeasuredExpectationValues;
    for (auto p : m_AllParticles) {
      if (!p->IsUnmeasured()) {
	Coordinates ExpValues = p->GetParametrisation()->GetExpectation(p->GetCoordinates());
	for (unsigned int i=0; i<ExpValues.size(); i++)
	  if (!p->IsParameterFixed(i)) 
	    ExpectationValues.push_back(ExpValues[i]);
      }
      else {
	Coordinates ExpValues = p->GetCoordinates(); //Unmeasured particles takes input internal representation as initial values
	for (unsigned int i=0; i<ExpValues.size(); i++)
	  if (!p->IsParameterFixed(i)) 
	    UnmeasuredExpectationValues.push_back(ExpValues[i]);
      }
    }
    if (ExpectationValues.size()!=m_NumberOfParameters) {
      result.status=statuscode::InconsistentExpectationValues;
      m_State=statuscode::InconsistentExpectationValues;
      return result;
    }
    if (m_NumberOfConstraintParameters>0) {
      for (auto c : m_ListOfConstraints) 
	if (c.HasProbDistFunc()) ExpectationValues.push_back(c.GetProbDistFunc()->ExpectationValue());
    }
    m_Covariance = Matrix(m_NumberOfParameters+m_NumberOfConstraintParameters, Coordinates (m_NumberOfParameters+m_NumberOfConstraintParameters, 0.0));
    bool nondiag=false;
    unsigned int Nrow=0;
    unsigned int Ncol=0;
    for (auto p : m_AllParticles) {
      if (p->IsUnmeasured()) continue;
      unsigned int Nparam = p->TotalNumberOfParameters();
      Matrix matrix = p->GetParametrisation()->GetExpectationCovMatrix(p->GetCoordinates());
      for (unsigned int j = 0; j<Nparam; j++) {
	if (p->IsParameterFixed(j)) continue;
	unsigned int Nofp=0;
	Ncol=Nrow;
	for (unsigned int k=j; k<Nparam; k++) {
	  if (p->IsParameterFixed(k)) continue;
	  if (j!=k && matrix[j][k]!=0.0) nondiag=true;
	  m_Covariance[Nrow][Ncol+Nofp]=matrix[j][k];
	  m_Covariance[Ncol+Nofp][Nrow]=matrix[j][k]; //Force covariance matrix to be symmetric
	  Nofp++;
	}
	Nrow++;
      }
    }
    for (unsigned int i=m_NumberOfParameters; i<m_NumberOfParameters+m_NumberOfConstraintParameters; i++)
      m_Covariance[i][i]=1.0;
    //std::cout <<"v = "; MatrixAlgebra::print(m_Covariance);
    Matrix InverseCovariance;
    if (nondiag) InverseCovariance = MatrixAlgebra::InverseMatrix(&m_Covariance);
    else {
      InverseCovariance = m_Covariance;
      for (unsigned int i=0; i<m_NumberOfParameters; i++) InverseCovariance[i][i] = 1.0/InverseCovariance[i][i];
    }
    if (!MatrixAlgebra::CheckMatrix(InverseCovariance)) {
      result.status=statuscode::NonInvertibleCovMatrix;
      m_State=statuscode::NonInvertibleCovMatrix;
      return result;
    }
    //Setting diagonal for constraint parameters to zero such that there is no contribution to chi2
    for (unsigned int i=m_NumberOfParameters; i<m_NumberOfParameters+m_NumberOfConstraintParameters; i++)
      InverseCovariance[i][i]=0.0; 
    //Setup of constraint pdf parameters

    //Remember to save input particle object's momenta
    //Init Fitting Procedure
    m_ParameterValues=ExpectationValues;
    m_UnmeasuredParameterValues=UnmeasuredExpectationValues;
    double chi2min = 99999999.0;
    double chi2 = 99999999.0;
    for (unsigned int iteration=0; iteration<m_MaxIterations; iteration++) {
      statuscode BAstate = CalculateDerivativeMatrix();
      if (BAstate!=statuscode::success) {
	result.status=BAstate;
	m_State=BAstate;
	return result;
      }
      if (m_NumberOfConstraintParameters>0) {
	unsigned int NProbDist = 0;
	for (unsigned int i=0; i<m_NumberOfConstraints; i++){
	  if (m_ListOfConstraints[i].HasProbDistFunc()) {
	    m_ConstraintValues[i]-=m_ParameterValues[m_NumberOfParameters+NProbDist]; //subtract ProbDist Parameter Value (initial value is most probable value)
	    NProbDist++;
	  }
	}
      }
      //compute global constrain fulfillment estimator
      double constraintSq = std::inner_product( m_ConstraintValues.begin(),  m_ConstraintValues.end(),  m_ConstraintValues.begin(), 0.0);
      if (m_NumberOfConstraintParameters>0) {
	unsigned int NProbDist=m_NumberOfParameters; 
	for (auto c : m_ListOfConstraints)
	  if (c.HasProbDistFunc()) {
#ifndef HACK
	    ExpectationValues[NProbDist] = m_ParameterValues[NProbDist]; 
#endif
	    //std::cout << "Gaus deriv = " << c.GetProbDistFunc()->DerivativeLog(m_ParameterValues[NProbDist]) << " Gaus 2nd deriv = " << c.GetProbDistFunc()->SecondDerivativeLog(m_ParameterValues[NProbDist]) << " param val = " << m_ParameterValues[NProbDist] << std::endl;
	    m_Covariance[NProbDist][NProbDist] = 2.0/(c.GetProbDistFunc()->SecondDerivativeLog(m_ParameterValues[NProbDist]));
	    NProbDist++;
	  }
      } 
      
      //now start the tedious matrix calculations
      
      //Full calculation work with non diagonal matrix
      
      //compute vbt = V*B{T} FOR SPECIAL CASE
      Matrix vbt = MatrixAlgebra::MatrixMultiplication(m_Covariance,MatrixAlgebra::TransposeMatrix(m_Bmatrix));
      //std::cout <<"b = "; MatrixAlgebra::print(m_Bmatrix);
      //std::cout <<"a = "; MatrixAlgebra::print(m_Amatrix);
      //std::cout <<"v = "; MatrixAlgebra::print(m_Covariance);
      //std::cout <<"vbt = "; MatrixAlgebra::print(vbt);
      //compute bvbt = B*VB{T}
      Matrix bvbt = MatrixAlgebra::MatrixMultiplication(m_Bmatrix,vbt);
      //std::cout <<"bvbt = "; MatrixAlgebra::print(bvbt);
      //invert bvb{T}
      Matrix Wb = MatrixAlgebra::InverseMatrix(&bvbt);
      //Chech if matrix could be inverted. In some cases, one row of b is zero (essentially at the first step, when b{i} and c{i} are equal to zero)
      if (!MatrixAlgebra::CheckMatrix(Wb)) {
	result.status=statuscode::NonInvertibleMatrix;
	m_State=statuscode::NonInvertibleMatrix;
	return result;
      }
      //std::cout <<"Wb = "; MatrixAlgebra::print(Wb);
      //std::cout <<"Unit = "; MatrixAlgebra::print(MatrixAlgebra::MatrixMultiplication(Wb,bvbt));
      //compute VB{T}(BVB{T})-1 
      m_vbtWb = MatrixAlgebra::MatrixMultiplication(vbt,Wb);
      //std::cout <<"vbtWb = "; MatrixAlgebra::print(m_vbtWb);
      //Recalculate VB{T}(BVB{T})-1 to include unmeasured parameters
      if (m_NumberOfFreeParameters>0) {
	//compute WbA where Wb=(BVB{T})-1
	Matrix WbA = MatrixAlgebra::MatrixMultiplication(Wb,m_Amatrix);
	//std::cout <<"WbA = "; MatrixAlgebra::print(WbA);
	//compute At = A{T} to be used several times
	Matrix At = MatrixAlgebra::TransposeMatrix(m_Amatrix);
	//std::cout <<"At = "; MatrixAlgebra::print(At);
	//std::cout <<"At dim = " << At.size() << " " << At.at(0).size() << " WbA dim = " << WbA.size() << " " << WbA.at(0).size() << std::endl;
	//compute Wa^-1 =(A{T}WbA)^-1 for general case
	Matrix mult = MatrixAlgebra::MatrixMultiplication(At,WbA);
	m_WaInv = MatrixAlgebra::InverseMatrix(&mult);
	//std::cout <<"WaInv = "; MatrixAlgebra::print(m_WaInv);
	//Chech if matrix could be inverted. In some cases, one row of b is zero (essentially at the first step, when b{i} and c{i} are equal to zero)
	if (!MatrixAlgebra::CheckMatrix(m_WaInv)) {
	  result.status=statuscode::NonInvertibleMatrix;
	  m_State=statuscode::NonInvertibleMatrix;
	  return result;
	}
	//COMPUTE WA^-1 A{T} WB (SAVED IN WAWB)
	m_WAWB = MatrixAlgebra::MatrixMultiplication(m_WaInv,MatrixAlgebra::MatrixMultiplication(At,Wb)); 
	//std::cout <<"WAWB = "; MatrixAlgebra::print(m_WAWB);
	//CALCULATE NEW vbtWb = vbtWb*(1-A*WAWB)
	m_vbtWb = MatrixAlgebra::MatrixSum(m_vbtWb,MatrixAlgebra::MatrixMultiplication(m_vbtWb,MatrixAlgebra::MatrixMultiplication(m_Amatrix,m_WAWB)),-1.0);
	//std::cout <<"vbtWb2 = "; MatrixAlgebra::print(m_vbtWb);
      }
      //update parameters values for constraint parameters
#ifndef HACK
      for (unsigned int i=0; i<m_NumberOfConstraintParameters; i++) {
	m_ParameterValues[m_NumberOfParameters+i] = -(m_ListOfConstraints[m_HasProbDistFunc[i]].GetProbDistFunc()->DerivativeLog(m_ParameterValues[m_NumberOfParameters+i]))/(m_ListOfConstraints[m_HasProbDistFunc[i]].GetProbDistFunc()->SecondDerivativeLog(m_ParameterValues[m_NumberOfParameters+i]))+ExpectationValues[m_NumberOfParameters+i]; //matches sign convention & adding expectationvalues to cancel later in c vector calculation so that the reference value is zero for x0
      }
#endif
      //compute c(y{l}) = A(a(l)-a(0))+B(y{l}-y{0})-f{y{l}}
      Coordinates cvector = 
	MatrixAlgebra::VectorSum(
				MatrixAlgebra::VectorSum(
							MatrixAlgebra::MatrixVectorMultiplication(m_Amatrix,MatrixAlgebra::VectorSum(m_UnmeasuredParameterValues,UnmeasuredExpectationValues,-1.0)), 
							MatrixAlgebra::MatrixVectorMultiplication(m_Bmatrix,MatrixAlgebra::VectorSum(m_ParameterValues,ExpectationValues,-1.0)),
							1.0),
				m_ConstraintValues,-1.0);		
      //std::cout << "b*diff = "; MatrixAlgebra::print(MatrixAlgebra::MatrixVectorMultiplication(m_Bmatrix,MatrixAlgebra::VectorSum(m_ParameterValues,ExpectationValues,-1.0)));
      //std::cout <<"cvector = "; MatrixAlgebra::print(cvector);
      //std::cout <<"c val = "; MatrixAlgebra::print(m_ConstraintValues);
      //update parameters values for constraint parameters 
#ifndef HACK 
      for (unsigned int i=0; i<m_NumberOfConstraintParameters; i++)
	ExpectationValues[m_NumberOfParameters+i] =  m_ParameterValues[m_NumberOfParameters+i];
#endif
      //compute y{l+1} = y{0}+vb{T}(bvb{T})-1*c(y{l})
      m_ParameterValues = MatrixAlgebra::MatrixVectorMultiplication(m_vbtWb, cvector); //Delta ParamaterValues
      //std::cout <<"param val = "; MatrixAlgebra::print(m_ParameterValues);
      if (m_NumberOfFreeParameters>0) {
	//std::cout << "wawb "; MatrixAlgebra::print(MatrixAlgebra::MatrixVectorMultiplication(m_WAWB, cvector),true);
	m_UnmeasuredParameterValues = MatrixAlgebra::VectorSum(UnmeasuredExpectationValues,MatrixAlgebra::MatrixVectorMultiplication(m_WAWB, cvector),1.0);
	//std::cout <<"Unmeasured param val = "; MatrixAlgebra::print(m_UnmeasuredParameterValues); 
      }
      //Compute chi2
      chi2 = MatrixAlgebra::VectorMultiplication(m_ParameterValues,MatrixAlgebra::MatrixVectorMultiplication(InverseCovariance, m_ParameterValues)); //inverse covariance has zeros for constraint parameters! Only measured parameters in classic chi2
      //std::cout << "chi2 = " << chi2 << std::endl;
      m_ParameterValues = MatrixAlgebra::VectorSum(ExpectationValues,m_ParameterValues,1.0);
      //Now add contribution from pdfs subtracting value at expectation point
      for (unsigned int i=0; i<m_HasProbDistFunc.size(); i++) {
	ProbDistFunc* pdf = m_ListOfConstraints[m_HasProbDistFunc[i]].GetProbDistFunc();
	chi2+=-2.0*log(pdf->Value(m_ParameterValues[m_NumberOfParameters+i]))+2.0*log(pdf->Value(pdf->ExpectationValue()));
      }
      //std::cout << "chi2total = " << chi2 << std::endl;
      //std::cout <<"param val = "; MatrixAlgebra::print(m_ParameterValues);
      //Update Particle Momenta
      for (unsigned int i=0;i<m_AllParticles.size();i++) {
	//std::cout << "particle "; MatrixAlgebra::print(m_AllParticles[i]->GetCoordinates(abcd));
	if (m_AllParticles[i]->NumberOfFreeParameters()>0) { //Only update particles that have free parameters
	  if (m_AllParticles[i]->IsUnmeasured()) 
	    m_AllParticles[i]->UpdateCoordinates(Coordinates (m_UnmeasuredParameterValues.begin()+m_Particle2Index[i],m_UnmeasuredParameterValues.begin()+m_Particle2Index[i]+m_AllParticles[i]->NumberOfFreeParameters()));
	  else 
	    m_AllParticles[i]->UpdateCoordinates(Coordinates (m_ParameterValues.begin()+m_Particle2Index[i],m_ParameterValues.begin()+m_Particle2Index[i]+m_AllParticles[i]->NumberOfFreeParameters()));
	}
      }
      //Chi2 convergrence 
      double deltachi2 = abs(chi2-chi2min)/(chi2!=0.0 ? chi2 : 1.0);
      if (abs(deltachi2)<m_DeltaChiConvergence && constraintSq<m_ConstraintConvergence*m_ConstraintConvergence*double(m_NumberOfConstraints*m_NumberOfConstraints)) {
	result.status=statuscode::success;
        m_State=statuscode::success;
	result.chi2=chi2;
	result.Ndof = m_NumberOfConstraints-m_NumberOfFreeParameters;
	result.Niter = iteration+1;
	for (auto p : m_AllParticles) 
	  result.FittedParticles.push_back(ParticleObject(p->GetCoordinates(p->GetParametrisation()),p->GetParametrisation(),p->GetParametrisation()));
        return result;
      }
      chi2min=abs(chi2); //protection against numerical features
    }
    result.status=statuscode::MaxNumIter;
    m_State=statuscode::MaxNumIter;
    result.chi2=chi2;
    result.Ndof = m_NumberOfConstraints-m_NumberOfFreeParameters;
    result.Niter = m_MaxIterations+1;
    for (auto p : m_AllParticles)
      result.FittedParticles.push_back(ParticleObject(p->GetCoordinates(p->GetParametrisation()),p->GetParametrisation(),p->GetParametrisation()));
    return result;
  }

  statuscode FitKernel::CalculateDerivativeMatrix() {
    //Constraint derivative (df/dParameters=df/dp*dp/dParameters
    m_dpdParameters = Matrix(internalRepresentation->GetNumberOfParameters()*m_AllParticles.size()+m_NumberOfConstraintParameters,Coordinates (m_TotalNumberOfParticleParameters+m_NumberOfConstraintParameters,0.0));  
    unsigned int row = 0;
    unsigned int col = 0;
    for (unsigned int i=0; i<m_AllParticles.size(); i++) {
      ParticleObject* p = m_AllParticles[i];
      //Unmeasured particles have zero derivative
      Matrix matrix = p->GetParametrisation()->Derivative(p->GetParametrisation()->InvTransform(p->GetCoordinates()));
      for (unsigned int j=0; j<matrix.size(); j++) std::copy(matrix[j].begin(),matrix[j].end(),m_dpdParameters[row++].begin()+col);
      col+=matrix[0].size();
    }
    //std::cout << "dpdParam = "; MatrixAlgebra::print(m_dpdParameters);
    //compute df/dp and constraint function distance
    Matrix dfdp(m_NumberOfConstraints, Coordinates (internalRepresentation->GetNumberOfParameters()*m_AllParticles.size()+m_NumberOfConstraintParameters,0.0));
    for (unsigned int idx_con=0; idx_con<m_NumberOfConstraints; idx_con++) {
      m_ConstraintValues[idx_con] = (m_ListOfConstraints[idx_con].ConstraintFunction()-m_ListOfConstraints[idx_con].GetConstraintValue());
      //Warning: Fix for ProbDist Constraints -- has been fixed
      Coordinates derivative=m_ListOfConstraints[idx_con].ConstraintInternalDerivative();
      //std::cout << "deriv = "; MatrixAlgebra::print(derivative);
      for (unsigned int i=0; i<m_Constraint2Particle[idx_con].size(); i++) {
	std::copy(derivative.begin()+i*internalRepresentation->GetNumberOfParameters(), derivative.begin()+(i+1)*internalRepresentation->GetNumberOfParameters(), dfdp[idx_con].begin()+m_Constraint2Particle[idx_con][i]*internalRepresentation->GetNumberOfParameters());
      }
    }
    //std::cout << "dfdp = "; MatrixAlgebra::print(dfdp);
    //To do : check matrix
    //B+A matrix:
    Matrix BAmatrix = MatrixAlgebra::MatrixMultiplication(dfdp,m_dpdParameters);
    //std::cout << "BAmatrix = "; MatrixAlgebra::print(BAmatrix);
    unsigned int Acol =0, Bcol=0;
    col=0;
    for (unsigned int i = 0; i<m_AllParticles.size(); i++) {
      if (col>BAmatrix[0].size()) return statuscode::InconsistentDerivativeMatrix;
      if (m_AllParticles[i]->NumberOfFreeParameters()==0) {
	col+=m_AllParticles[i]->TotalNumberOfParameters();
	continue;
      }
      for (unsigned int j=0; j<m_AllParticles[i]->TotalNumberOfParameters(); j++) {
	if (m_AllParticles[i]->IsParameterFixed(j)) {
	  col++;
	  continue;
	}
	if (m_AllParticles[i]->IsUnmeasured()) {
	  for (unsigned int row=0; row<BAmatrix.size(); row++) 
	    m_Amatrix[row][Acol] = BAmatrix[row][col];	    
	  Acol++;
	}
	else {
	  for (unsigned int row=0; row<BAmatrix.size(); row++) {
	    m_Bmatrix[row][Bcol] = BAmatrix[row][col];
	  }
	  Bcol++;
	}
	col++;
      }
    }
    if (Bcol!=m_NumberOfParameters || Acol!=m_NumberOfFreeParameters) return statuscode::InconsistentDerivativeMatrix;
    //Add constraint parameters
    /*if (m_NumberOfConstraintParameters>0) {
      for (unsigned int i=0; i<m_NumberOfConstraintParameters; i++) {
	for (unsigned int row=0; row<BAmatrix.size(); row++)
	  m_Bmatrix[row][Bcol] = BAmatrix[row][col];
	Bcol++;
	col++;
      }
      } */
    return statuscode::success;
  }
  
  Matrix FitKernel::GetParameterCovMatrix() { //covariance matrix in chi2-parameters
    if (m_State!=statuscode::success) return Matrix();
    Matrix BV = MatrixAlgebra::MatrixMultiplication(m_Bmatrix,m_Covariance);
    //TO DO: test that m_vbtWb is symmetric i.e. m_vbtWb=?=m_vbtWb^T
    Matrix C11 = MatrixAlgebra::MatrixSum(m_Covariance,MatrixAlgebra::MatrixMultiplication(m_vbtWb,BV),-1.0);
    Matrix C21;
    Matrix C22;
    if (m_NumberOfFreeParameters>0) {
      C21 = MatrixAlgebra::MatrixMultiplication(m_WAWB,BV);
      C22 = m_WaInv;
    }
    unsigned int Dim = m_NumberOfParameters+m_NumberOfFreeParameters+m_NumberOfConstraintParameters;
    Matrix result(Dim, Coordinates(Dim,0.0));
    for (unsigned int row=0; row<Dim; row++)
      for (unsigned int col=0; col<Dim; col++) {
	if (row<C11.size() && col<C11.size()) result[row][col]=C11[row][col];
	if (row>=C11.size() && col<C11.size()) result[row][col]=-C21[row-C11.size()][col];
	if (row<C11.size() && col>=C11.size()) result[row][col]=-C21[col-C11.size()][row];
	if (row>=C11.size() && col>=C11.size()) result[row][col]=C22[row-C11.size()][col-C11.size()];
      }
    return result;
  }

  Matrix FitKernel::GetParticleCovMatrix() { //covariance matrix for fitted momenta in internal representation (PxPyPxM)
    if (m_State!=statuscode::success) return Matrix();
    //Sort dp/dParam matrix
    unsigned int Dim1 = m_AllParticles.size()*internalRepresentation->GetNumberOfParameters();
    unsigned int Dim2 = m_NumberOfParameters+m_NumberOfFreeParameters+m_NumberOfConstraintParameters;
    Matrix sorted(Dim1, Coordinates(Dim2,0.0));
    // Dim for dpdParam: Matrix(internalRepresentation->GetNumberOfParameters()*m_AllParticles.size()+NumberOfConstraintParameters,Coordinates (NumberOfParameters+NumberOfConstraintParameters,0.0));  
    for (unsigned int row=0; row<Dim1; row++) {
      unsigned int Acol =m_NumberOfParameters+m_NumberOfConstraintParameters, Bcol=0, col=0;
      for (unsigned int i=0; i<m_AllParticles.size(); i++) {
	for (unsigned int j=0; j<m_AllParticles[i]->TotalNumberOfParameters(); j++) {
	  if (!m_AllParticles[i]->IsParameterFixed(j)) {
	    if (m_AllParticles[i]->IsUnmeasured())
	      sorted[row][Acol++]=m_dpdParameters[row][col];
	    else 
	      sorted[row][Bcol++]=m_dpdParameters[row][col];
	  }
	  col++;
	}
      }
    }
    //std::cout << "Sorted = " ; MatrixAlgebra::print(sorted,true);
    return MatrixAlgebra::MatrixMultiplication(sorted,MatrixAlgebra::MatrixMultiplication(GetParameterCovMatrix(),MatrixAlgebra::TransposeMatrix(sorted)));
  }

  bool FitKernel::Initialize() {
    std::vector<Constraint*> AllConstraints; //Only allow unique constraints
    m_AllParticles.clear();
    m_Particle2Index.clear();
    m_Constraint2Particle.clear();
    m_ParameterValues.clear();
    m_UnmeasuredParameterValues.clear();
    m_NumberOfConstraints = m_ListOfConstraints.size();
    m_NumberOfFreeParameters = 0;
    m_NumberOfParameters = 0;
    m_NumberOfConstraintParameters = 0;
    m_TotalNumberOfParticleParameters = 0;
    for (unsigned int idx_con =0; idx_con<m_ListOfConstraints.size(); idx_con++) {
      if (m_ListOfConstraints[idx_con].HasProbDistFunc()) m_HasProbDistFunc.push_back(idx_con);
      m_Constraint2Particle.push_back(std::vector<unsigned int> ());
      for (auto c : m_ListOfConstraints[idx_con].GetListOfConstraints()) { //c=constraint
	if (std::find(AllConstraints.begin(),AllConstraints.end(), c.first) == AllConstraints.end()) AllConstraints.push_back(c.first);
	else return false;
	for (auto p : c.first->GetListOfParticleObjects()) {
	  unsigned int j = std::find(m_AllParticles.begin(), m_AllParticles.end(), p)-m_AllParticles.begin();  
	  if (j == m_AllParticles.size()) m_AllParticles.push_back(p);
	  if (std::find(m_Constraint2Particle[idx_con].begin(), m_Constraint2Particle[idx_con].end(), j) == m_Constraint2Particle[idx_con].end()) m_Constraint2Particle[idx_con].push_back(j);
	}
      }
    }
    for (unsigned int i=0; i<m_AllParticles.size(); i++) {
      m_TotalNumberOfParticleParameters+=m_AllParticles[i]->TotalNumberOfParameters();
      m_Particle2Index.push_back(m_AllParticles[i]->IsUnmeasured() ? m_NumberOfFreeParameters : m_NumberOfParameters);
      unsigned int NParam = m_AllParticles[i]->NumberOfFreeParameters();
      if (m_AllParticles[i]->IsUnmeasured()) {
	m_NumberOfFreeParameters += NParam;
      }
      else {
        m_NumberOfParameters += NParam;
      }
    }
    m_NumberOfConstraintParameters = m_HasProbDistFunc.size();
    m_Bmatrix = Matrix(m_NumberOfConstraints, Coordinates (m_NumberOfParameters+m_NumberOfConstraintParameters,0.0));
    m_Amatrix = Matrix(m_NumberOfConstraints, Coordinates (m_NumberOfFreeParameters,0.0));
    m_ConstraintValues = Coordinates (m_NumberOfConstraints, 0.0);
    unsigned int NProbDist = m_NumberOfParameters;
    for (unsigned int idx_con=0; idx_con<m_NumberOfConstraints; idx_con++) {
      if (m_ListOfConstraints[idx_con].HasProbDistFunc()) {
	m_Bmatrix[idx_con][NProbDist++] = -1.0; //sign convention for Probability dist value
      }
    }
    /*std::cout << "Testing ABCFIT: Constraint2Particle = "<< m_Constraint2Particle.size() << " "; MatrixAlgebra::print(m_Constraint2Particle);
    std::cout << "Testing ABCFIT: Particle2Index = "; MatrixAlgebra::print(m_Particle2Index); 
    std::cout << "Number of free parameters = " << m_NumberOfFreeParameters  << std::endl;
    std::cout << "Number of measured parameters = " << m_NumberOfParameters << std::endl;
    std::cout << "Number of constraint parameters = " << m_NumberOfConstraintParameters << std::endl;*/
    return true;
  }



}
