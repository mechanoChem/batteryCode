#include"../../include/model.h"


template <class T, int dim>
model<T,dim>::model(){std::cout<<"MODEL GENERATED"<<std::endl; }

template <class T, int dim>
model<T, dim>::~model (){}


template <class T, int dim>
void model<T, dim>::reinit(parametersClass& _params)
{
	params=&_params;
	youngModule_neg=params->getDouble("youngModule_neg");
	youngModule_sep=params->getDouble("youngModule_sep");
	youngModule_pos=params->getDouble("youngModule_pos");
	l_neg=params->getDouble("l_neg");
	l_sep=params->getDouble("l_sep");
	l_pos=params->getDouble("l_pos");
	
	
}

template <class T, int dim>
void model<T, dim>::refresh(double _currentTime, double _dt)
{
	currentTime=_currentTime;
	dt=_dt;
	
}

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;