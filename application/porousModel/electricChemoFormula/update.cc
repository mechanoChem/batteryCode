#include"porousModel/porousModelFormula.h"

template <class T,int dim>
void porousModelFormula<T,dim>::update(T _c_li ,T _c_li_plus, Table<1,T> _c_li_plus_grad, T _temp, T _phi_s, Table<1, T> _phi_s_grad, T _phi_e, Table<1, T> _phi_e_grad, T _defVol, int _domainflag, int _currentflag){
  
	this->afterUpdate=false;
	
	this->c_li=_c_li;
	this->c_li_plus=_c_li_plus;
	this->c_li_plus_grad=&_c_li_plus_grad;
	this->Temp=_temp;
	this->phi_s=_phi_s;
	
	this->phi_s_grad=&_phi_s_grad;
	this->phi_e=_phi_e;
	this->phi_e_grad=&_phi_e_grad;
	
	defVol=_defVol;

	this->domainflag=_domainflag;
	currentflag=_currentflag;
	
	double c_li_max;
	if(this->domainflag==1) {
		c_li_max=this->params->getDouble("c_li_max_pos");
	}
	else if(this->domainflag==-1) {
		c_li_max=this->params->getDouble("c_li_max_neg");
	}
	else if(this->domainflag==0) {
		c_li_max=1;
	}
	else {std::cout<<"wrong domain domainflag"<<std::endl; exit(-1);}
	
	
	this->UnitC=this->c_li/c_li_max;
	
	//formula_SolidVolFracion();
	//formula_PoreVolFraction();
	this->formula_porosity(this->UnitC, this->Temp, defVol, this->domainflag, currentflag);
	this->l_a=3.0/this->l_R_e;
	this->formula_Q_ohm();
	this->formula_Ke();
	this->formula_D_l();
	
	if(this->domainflag!=0){
		this->c_li_surface_parabolic();
		this->UnitC_surface=this->l_c_li_surface/c_li_max;
		this->l_a=3.0/this->l_R_e;
		this->formula_jn();
		this->formula_Q_rev();
		this->formula_Q_rxn();
	}	
	this->afterUpdate=true;
}

template class porousModelFormula<Sacado::Fad::DFad<double>, 1>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 2>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 3>;