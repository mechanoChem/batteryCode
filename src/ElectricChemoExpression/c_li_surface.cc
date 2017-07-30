#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
void electricChemo<T,dim>::c_li_surface_parabolic()
{
	double D_s, R_s_0;

	if(domainflag==1) {D_s=params->getDouble("D_s_pos");R_s_0=params->getDouble("R_s_0_pos");}
	else if(domainflag==-1) {D_s=params->getDouble("D_s_neg");R_s_0=params->getDouble("R_s_0_neg");}
	else {std::cout<<"wrong domain domainflag for D_s"<<std::endl; exit(-1);}
	l_c_li_surface=c_li-R_s_0*l_jn/5/D_s;
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;