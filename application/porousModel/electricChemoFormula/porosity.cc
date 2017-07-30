#include"porousModel/porousModelFormula.h"

template <class T,int dim>
void porousModelFormula<T,dim>::formula_SolidVolFracion(T UnitC, T Temp, T defVol, int Tdomainflag, int Tcurrentflag){
	double kappa, kappa_s, eps_s_0, Pl, Pb, R_s_0;
	if(Tdomainflag==1) {
		kappa=this->params->getDouble("kappa_pos");
		eps_s_0=this->params->getDouble("eps_s_0_pos");
		R_s_0=this->params->getDouble("R_s_0_pos");
	}
	else if(Tdomainflag==-1) {
		kappa=this->params->getDouble("kappa_neg");
		eps_s_0=this->params->getDouble("eps_s_0_neg");
		R_s_0=this->params->getDouble("R_s_0_neg");
	}
	else if(Tdomainflag==0){
		kappa=this->params->getDouble("kappa_sep");
		eps_s_0=this->params->getDouble("eps_s_0_sep");		
		R_s_0=this->params->getDouble("R_s_0_sep");	
	}
	kappa_s=this->params->getDouble("kappa_s");
	Pl=this->params->getDouble("pl");
	Pb=this->params->getDouble("pb");
	formula_Beta_s(UnitC, Tdomainflag, Tcurrentflag);
	formula_Beta(UnitC, Tdomainflag, Tcurrentflag);
	formula_Beta_s_T(Temp,Tdomainflag);
	formula_Beta_T(Temp,Tdomainflag);
	
	T factor=((kappa*(defVol/(1+l_beta)/(1+l_beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+l_beta_s)*(1+l_beta_s_T);
	this->l_eps_s=factor/defVol*eps_s_0;
	this->l_R_e=std::pow(factor,1/3)*R_s_0;
}

template <class T,int dim>
void porousModelFormula<T,dim>::formula_PoreVolFraction(T defVol, int Tdomainflag){
	double eps_b_0;
	if(Tdomainflag==1) {
		eps_b_0=this->params->getDouble("eps_b_0_pos");
	}
	else if(Tdomainflag==-1) {
		eps_b_0=this->params->getDouble("eps_b_0_neg");
	}
	else if(Tdomainflag==0){
		eps_b_0=this->params->getDouble("eps_b_0_sep");		
	}
	this->l_eps_l=1-this->l_eps_s-eps_b_0/defVol;
}

template <class T,int dim>
void porousModelFormula<T,dim>::formula_porosity(T UnitC, T Temp, T defVol, int Tdomainflag, int Tcurrentflag){
	formula_SolidVolFracion(UnitC,Temp, defVol, Tdomainflag, Tcurrentflag);
	formula_PoreVolFraction(defVol, Tdomainflag);
}

template class porousModelFormula<Sacado::Fad::DFad<double>, 1>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 2>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 3>;
