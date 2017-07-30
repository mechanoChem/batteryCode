#include"../../include/model.h"

template <class T, int dim>
model<T, dim>::model(){};

template <class T, int dim>
model<T, dim>::~model(){};

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;