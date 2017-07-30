#include"../../../include/porousModel/porousModelResidual.h"

template <class T, int dim>
porousModelResidual<T, dim>::porousModelResidual():model<T,dim>()
{
	std::cout<<"battery porous model generated"<<std::endl;
}


template <class T, int dim>
porousModelResidual<T, dim>::~porousModelResidual (){}


template class porousModelResidual<Sacado::Fad::DFad<double>, 1>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 2>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 3>;

