#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
void electricChemo<T,dim>::formula_Q_ohm()
{
	l_Q_ohm=0;
	double se, Rr, F, t_0;
	Rr=params->getDouble("Rr");
	F=params->getDouble("F");
	t_0=params->getDouble("t_0");
	
	if(domainflag==1) se=params->getDouble("se_pos");
	else if(domainflag==-1) se=params->getDouble("se_neg");
  for (unsigned int j = 0; j < dim; j++){
  	l_Q_ohm+=se*l_eps_s*(*phi_s_grad)[j]*(*phi_s_grad)[j]+(l_Ke*std::pow(l_eps_l, 1.5)*((*phi_e_grad)[j]-2.0*Rr*Temp/F*(1.0-t_0)/c_li_plus*((*c_li_plus_grad)[j])))*(*phi_e_grad)[j];
  }
}

template <class T,int dim>
void electricChemo<T,dim>::formula_Q_rev()
{
	double F;
	F=params->getDouble("F");
	formula_dUdt();
	l_Q_rev=F*l_a*l_jn*Temp*l_dUdt;
}

template <class T,int dim>
void electricChemo<T,dim>::formula_Q_rxn()
{
	double F;
	F=params->getDouble("F");
	T eta = phi_s-phi_e-l_Usc;
	l_Q_rxn=F*l_jn*eta;
}

template <class T,int dim>
void electricChemo<T,dim>::formula_dUdt()
{
	if(UnitC<=0.2)  l_dUdt=0.01442*UnitC*UnitC-0.00291*UnitC-0.000138;
	else if(UnitC<=0.4) l_dUdt=0.00634*UnitC*UnitC*UnitC-0.006625*UnitC*UnitC+0.002635*UnitC-0.0004554;
	else if(UnitC<=0.5) l_dUdt=0.001059*UnitC-0.0004793;
	else if(UnitC<=0.7) l_dUdt=0.00025*UnitC-7.5e-5;
	else if(UnitC<=0.8) l_dUdt=-0.001*UnitC+0.0008;
	else if(UnitC<=0.85) l_dUdt=0.0333*UnitC*UnitC-0.057*UnitC+0.02427;
	else if(UnitC<=0.95) l_dUdt=0.002*UnitC*UnitC-0.0039*UnitC+0.00177;
  else if (UnitC<=1) l_dUdt=-0.0014*UnitC+0.0012;
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;
