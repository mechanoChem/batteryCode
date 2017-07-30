#ifndef electricChemo_h
#define electricChemo_h

#include "supplementary/parameters.h"
#include "supplementary/dataStruct.h"
#include "supplementary/functionEvaluations.h"
#include <Sacado.hpp>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>


template <class T, int dim>
class electricChemo
{
public:
  electricChemo ();
  ~electricChemo();
	
	parametersClass* params;
	void reinit(parametersClass& _params);
	virtual void update(T _c_li ,T _c_li_plus, Table<1,T> _c_li_plus_grad, T _Temp, T _phi_s, Table<1, T> _phi_s_grad, T _phi_e, Table<1, T> _phi_e_grad, int _domainflag);
	
	/*
	*expression, formula and equations
	*/
	
	virtual void formula_jn();
	virtual void formula_j0();
	virtual void formula_Usc(T UnitC, int Tdomainflag);
	
	virtual void c_li_surface_parabolic();
	
	virtual void formula_Q_ohm();
	virtual void formula_Q_rev();
	virtual void formula_Q_rxn();
	virtual void formula_dUdt();
	
	virtual void formula_Ke();
	virtual void formula_D_l();


	/*
	*1:positive electrode, 0:separator, -1 negative electrode
	*/
	int domainflag;
	
	bool afterUpdate;
	
	T c_li, c_li_plus, Temp, phi_s, phi_e;
	
	Table<1,T>* phi_s_grad;
	Table<1,T>* phi_e_grad;
	Table<1,T>* c_li_plus_grad;
	
	/*
	*variables
	*/
	T UnitC, UnitC_surface;
	T l_beta_s, l_beta, l_beta_s_T, l_beta_T;
	
	T l_a, l_R_e, l_eps_s, l_eps_l;
	
	T l_c_li_surface;
	
	T l_jn, l_j0, l_Usc;
	
	T l_Q_ohm, l_Q_rev, l_Q_rxn, l_dUdt;
	
	T l_Ke, l_D_l;
	
	/*
	*return value of formula
	*/
	T Beta_s(bool afterUpdateCheck=false);
	T Beta(bool afterUpdateCheck=false);
	T Beta_s_T(bool afterUpdateCheck=false);
	T Beta_T(bool afterUpdateCheck=false);
	
	T R_e(bool afterUpdateCheck=false);
	T eps_s(bool afterUpdateCheck=false);
	T eps_l(bool afterUpdateCheck=false);
	
	T jn(bool afterUpdateCheck=false);
	T Usc(bool afterUpdateCheck=false);
	
	
	T Q_ohm(bool afterUpdateCheck=false);
	T Q_rev(bool afterUpdateCheck=false);
	T Q_rxn(bool afterUpdateCheck=false);
	
	T dUdt(bool afterUpdateCheck=false);
	
	T Ke(bool afterUpdateCheck=false);
	T D_l(bool afterUpdateCheck=false);
	
	
};

#endif