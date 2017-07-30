#include"../../include/model.h"

//default evaluateStress No thermal chemical coupling

template <class T, int dim>
void model<T,dim>::scalling(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& R, double scallingFactor)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			R[i] = R[i]*scallingFactor;
    }
  }
	
}

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;