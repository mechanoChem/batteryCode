#include"porousModel/porousModelFormula.h"

template <class T,int dim>
porousModelFormula<T,dim>::porousModelFormula(){}

template <class T,int dim>
porousModelFormula<T,dim>::~porousModelFormula(){}


template class porousModelFormula<Sacado::Fad::DFad<double>,1>;
template class porousModelFormula<Sacado::Fad::DFad<double>,2>;
template class porousModelFormula<Sacado::Fad::DFad<double>,3>;