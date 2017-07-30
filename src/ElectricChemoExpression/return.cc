#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
T electricChemo<T,dim>::R_e(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_R_e;
}

template <class T,int dim>
T electricChemo<T,dim>::eps_s(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_eps_s;
}

template <class T,int dim>
T electricChemo<T,dim>::eps_l(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_eps_l;
}

template <class T,int dim>
T electricChemo<T,dim>::jn(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_jn;
}

template <class T,int dim>
T electricChemo<T,dim>::Usc(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_Usc;
}

template <class T,int dim>
T electricChemo<T,dim>::Q_ohm(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_Q_ohm;
}

template <class T,int dim>
T electricChemo<T,dim>::Q_rev(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_Q_rev;
}

template <class T,int dim>
T electricChemo<T,dim>::Q_rxn(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_Q_rxn;
}

template <class T,int dim>
T electricChemo<T,dim>::D_l(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_D_l;
}

template <class T,int dim>
T electricChemo<T,dim>::Ke(bool afterUpdateCheck){
	if(afterUpdateCheck==true){
		if(afterUpdate==false){std::cout<<"electric chemo return value before update"; exit(-1);}
	}
	return l_Ke;
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;
