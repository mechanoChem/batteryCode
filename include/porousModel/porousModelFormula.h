#ifndef porousModelFormula_h
#define porousModelFormula_h

#include "../electricChemo/electricChemo.h"

template <class T, int dim>
class porousModelFormula : public electricChemo<T, dim>
{
	 public:
	porousModelFormula();
	~porousModelFormula();
	
	void update(T _c_li ,T _c_li_plus, Table<1,T> _c_li_plus_grad, T _temp, T _phi_s, Table<1, T> _phi_s_grad, T _phi_e, Table<1, T> _phi_e_grad, T _defMap, int _domainflag, int _currentflag);
		
	void formula_Beta_s(T UnitC, int Tdomainflag, int Tcurrentflag);
	void formula_Beta(T UnitC,int Tdomainflag, int Tcurrentflag);
	void formula_Beta_s_T(T Temp, int Tdomainflag);
	void formula_Beta_T(T Temp, int Tdomainflag);
	
	void formula_SolidVolFracion(T UnitC, T Temp, T defVol, int Tdomainflag, int Tcurrentflag);
	void formula_PoreVolFraction(T defVol, int Tdomainflag);
	void formula_porosity(T UnitC, T Temp, T defVol, int Tdomainflag, int Tcurrentflag);
	
	int currentflag;
	
	T l_beta_s, l_beta, l_beta_s_T, l_beta_T, defVol;
	
	
	/*
	*return value of formula
	*/
	T Beta_s(bool afterUpdateCheck=false);
	T Beta(bool afterUpdateCheck=false);
	T Beta_s_T(bool afterUpdateCheck=false);
	T Beta_T(bool afterUpdateCheck=false);
};

#endif