#include"porousModel/porousModelFormula.h"

template <class T,int dim>
T porousModelFormula<T,dim>::Beta_s(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(this->afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_beta_s;
}

template <class T,int dim>
T porousModelFormula<T,dim>::Beta(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(this->afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_beta;
}

template <class T,int dim>
T porousModelFormula<T,dim>::Beta_s_T(bool afterUpdateCheck){
  if(afterUpdateCheck==true){
		if(this->afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_beta_s_T;
}

template <class T,int dim>
T porousModelFormula<T,dim>::Beta_T(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(this->afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_beta_T;
}
template class porousModelFormula<Sacado::Fad::DFad<double>, 1>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 2>;
template class porousModelFormula<Sacado::Fad::DFad<double>, 3>;