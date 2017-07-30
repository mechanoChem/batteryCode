#ifndef model_h
#define model_h
#include <Sacado.hpp>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/table.h>
#include "supplementary/dataStruct.h"
#include "supplementary/parameters.h"
#include "supplementary/supplementaryFunctions.h"
#include "supplementary/functionEvaluations.h"

template <class T, int dim>
class model
{
public:
  model ();
  ~model();
	
	parametersClass* params;
	
	void reinit(parametersClass& _params);
	void refresh(double _currentTime, double _dt);
		
	//mechanics
	//constitutive model:Saint-Venant Kirchhoff as defaults
	virtual double SVK3D(unsigned int i, unsigned int j, unsigned int k, unsigned int l, double E);
	virtual double SVK2D(unsigned int i, unsigned int j, double E);
	virtual void evaluateStrain(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap);
	virtual void evaluateStress(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap);
	void residualForMechanics(const FEValues<dim>& _fe_values, unsigned int DOF, deformationMap<T, dim>& defMap, Table<1, T >& R);
	virtual void residualForNewmmanBC(Table<1, T >& R);
	
	void scalling(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& R, double scallingFactor);
	
	Table<3, T > P;
	Table<3, T > Fe;
	Table<3, T > E;
	double currentTime, dt;
	//unsigned int iteration;
	//material properity
	double nu, youngModule_neg, youngModule_sep, youngModule_pos;
	double l_neg, l_sep, l_pos;
	
};

#endif