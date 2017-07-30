#include"ElectricChemo/electricChemo.h"

template <class T,int dim>
electricChemo<T,dim>::electricChemo(){}

template <class T,int dim>
electricChemo<T,dim>::~electricChemo(){}

template <class T,int dim>
void electricChemo<T,dim>::reinit(parametersClass& _params)
{
	params=&_params;
}


template class electricChemo<Sacado::Fad::DFad<double>,1>;
template class electricChemo<Sacado::Fad::DFad<double>,2>;
template class electricChemo<Sacado::Fad::DFad<double>,3>;