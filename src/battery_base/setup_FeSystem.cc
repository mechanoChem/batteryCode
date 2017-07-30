#include"../../include/battery.h"

template <int dim>
void battery<dim>::setup_FeSystem()
{
  electrode_fe.reset(new FESystem<dim>(FE_Q<dim>(1),dim,
	             												 FE_Q<dim>(1),1,
	       			 												 FE_Q<dim>(1),2,
	       			 												 FE_Q<dim>(1),1,
	       			 												 FE_Q<dim>(1),1));
  electrolyte_fe.reset(new FESystem<dim>(FE_Q<dim>(1),dim,
		 						 												 FE_Nothing<dim>(),1,
		 						 											   FE_Q<dim>(1),2,
		 						 												 FE_Nothing<dim>(),1,
		 						 												 FE_Q<dim>(1),1));
																				 
	fe_collection.push_back (*electrode_fe);
	fe_collection.push_back (*electrolyte_fe);
  q_collection.push_back (electrode_quadrature);
  q_collection.push_back (electrolyte_quadrature);
								 	
}

template class battery<1>;
template class battery<2>;
template class battery<3>;