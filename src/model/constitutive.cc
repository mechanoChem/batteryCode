#include"../../include/model.h"

template <class T, int dim>
double model<T,dim>::SVK3D(unsigned int i, unsigned int j, unsigned int k, unsigned int l,double E)
{
	double nu=params->getDouble("nu");
  double lambda=(E*nu)/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
}

template <class T, int dim>
double model<T,dim>::SVK2D(unsigned int i, unsigned int j,double E)
{
  double lambda=(E*nu)/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
  if (i==j && i<2) return lambda + 2*mu;
  else if (i==2 && j==2) return mu;
  else if ((i+j)==1) return lambda;
  else return 0.0;
}

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;