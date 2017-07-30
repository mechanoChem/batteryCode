#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
void electricChemo<T,dim>::formula_Ke()
{
	l_Ke=((34.5*exp(-798/Temp)*std::pow((c_li_plus*1.0e3),3)-485*exp(-1080/Temp)*std::pow((c_li_plus*1.0e3),2)+2440*exp(-1440/Temp)*(c_li_plus*1.0e3))/10)*1.0e6;
}

template <class T,int dim>
void electricChemo<T,dim>::formula_D_l()
{
	l_D_l=std::pow(10,(-4.43-54/(Temp-229-5e3*c_li_plus)-2.2e2*c_li_plus))*1.0e8;
}

template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;