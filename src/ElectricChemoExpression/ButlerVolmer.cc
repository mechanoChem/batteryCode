#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
void electricChemo<T,dim>::formula_j0()
{
	double k, alpha, c1max;
	if(domainflag==1) {k=params->getDouble("k_pos");alpha=params->getDouble("alpha_pos");c1max=params->getDouble("c_li_max_pos");}
	else if(domainflag==-1) {k=params->getDouble("k_neg");alpha=params->getDouble("alpha_neg");c1max=params->getDouble("c_li_max_neg");}
	
	l_j0=k*std::pow(1000,0.5)*c1max*std::pow(c_li_plus*1.0e3,alpha)*std::pow((1-UnitC_surface),alpha)*std::pow(UnitC_surface,(1-alpha));
}

template <class T,int dim>
void electricChemo<T,dim>::formula_jn()
{
	if(domainflag==0)l_jn=0;
	else{
		double alpha, F, Rr;
		F=params->getDouble("F");
		Rr=params->getDouble("Rr");
		if(domainflag==1) alpha=params->getDouble("alpha_pos");
		else if(domainflag==-1) alpha=params->getDouble("alpha_neg");
	
		formula_j0();
		formula_Usc(UnitC, domainflag);
		T eta = phi_s-phi_e-l_Usc;
		l_jn=l_a*l_j0*(exp(alpha*F/Rr/Temp*(eta))-exp(-alpha*F/Rr/Temp*(eta)));
	}
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;