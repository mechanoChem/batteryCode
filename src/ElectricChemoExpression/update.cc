#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
void electricChemo<T,dim>::update(T _c_li ,T _c_li_plus, Table<1,T> _c_li_plus_grad, T _Temp, T _phi_s, Table<1, T> _phi_s_grad, T _phi_e, Table<1, T> _phi_e_grad, int _domainflag){
  
	afterUpdate=false;
	double c_li_max;
	if(domainflag==1) {
		l_eps_s=params->getDouble("eps_s_pos");
		l_eps_l=params->getDouble("eps_l_pos");
		l_R_e=params->getDouble("R_e_pos");
		c_li_max=params->getDouble("c_li_max_pos");
	}
	else if(domainflag==-1) {
		l_eps_s=params->getDouble("eps_s_neg");
		l_eps_l=params->getDouble("eps_l_neg");
		l_R_e=params->getDouble("R_e_neg");
		c_li_max=params->getDouble("c_li_max_neg");
	}
	else {std::cout<<"wrong domain domainflag"<<std::endl; exit(-1);}
	
	c_li=_c_li;
	c_li_plus=_c_li_plus;
	c_li_plus_grad=&_c_li_plus_grad;
	Temp=_Temp;
	phi_s=_phi_s;
	phi_s_grad=&_phi_s_grad;
	phi_e=_phi_e;
	phi_e_grad=&_phi_e_grad;
	
	UnitC=c_li/c_li_max;
	c_li_surface_parabolic();
	UnitC_surface=l_c_li_surface/c_li_max;
	domainflag=_domainflag;
	
	l_a=3.0/l_R_e;
	
	//formula_Usc();
	//formula_j0();
	formula_jn();
	
	formula_Q_ohm();
	formula_Q_rev();
	formula_Q_rxn();
	
	afterUpdate=true;
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;