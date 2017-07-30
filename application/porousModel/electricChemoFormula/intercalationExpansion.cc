#include"porousModel/porousModelFormula.h"

template <class T ,int dim>
void porousModelFormula<T,dim>::formula_Beta_s(T UnitC, int Tdomainflag, int Tcurrentflag)
{
	if(Tdomainflag==-1) {
		l_beta_s=1.496*std::pow(UnitC,3)-1.739*UnitC*UnitC+1.02*UnitC-0.03304*std::exp(2.972*UnitC)-0.04587*tanh((UnitC-0.1)/0.1)-0.003608*tanh((UnitC-0.3)/0.1)+0.0214*tanh((UnitC-0.65)/0.1);
		if(Tcurrentflag==1) l_beta_s=l_beta_s-0.0939;
	}
	else l_beta_s=0;
}

template <class T ,int dim>
void porousModelFormula<T,dim>::formula_Beta(T UnitC, int Tdomainflag, int Tcurrentflag)
{
	if(Tdomainflag==-1) {
		l_beta=0.1315*std::pow(UnitC,4)-0.1989*std::pow(UnitC,3)+0.06481*UnitC*UnitC+0.02793*UnitC-0.000655;
		if(Tcurrentflag==1) l_beta=l_beta-0.019;
	}
	else l_beta=0;
}

template <class T ,int dim>
void porousModelFormula<T,dim>::formula_Beta_s_T(T Temp, int Tdomainflag)
{
	double T_0=this->params->getDouble("T_0");
	double omega_s=this->params->getDouble("omega_s");
	l_beta_s_T=(Temp-T_0)*omega_s;
}

template <class T ,int dim>
void porousModelFormula<T,dim>::formula_Beta_T(T Temp, int Tdomainflag)
{
	double omega;
	if(Tdomainflag==1) omega=this->params->getDouble("omega_pos");
	else if(Tdomainflag==0) omega=this->params->getDouble("omega_sep");
	else if(Tdomainflag==-1) omega=this->params->getDouble("omega_neg");
	double T_0=this->params->getDouble("T_0");
	l_beta_T=(Temp-T_0)*omega; 

}

template class porousModelFormula<Sacado::Fad::DFad<double>, 1>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 2>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 3>;
