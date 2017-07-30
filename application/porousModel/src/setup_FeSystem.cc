#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::setup_FeSystem()
{
  battery<dim>::electrode_fe.reset(new FESystem<dim>(FE_Q<dim>(1),dim,
	             												 FE_Q<dim>(1),1,
	       			 												 FE_Q<dim>(1),2,
	       			 												 FE_Q<dim>(1),1,
	       			 												 FE_Q<dim>(1),1));
  battery<dim>::electrolyte_fe.reset(new FESystem<dim>(FE_Q<dim>(1),dim,
		 						 												 FE_Nothing<dim>(),1,
		 						 											   FE_Q<dim>(1),2,
		 						 												 FE_Nothing<dim>(),1,
		 						 												 FE_Q<dim>(1),1));
																				 
	battery<dim>::fe_collection.push_back (*battery<dim>::electrode_fe);
	battery<dim>::fe_collection.push_back (*battery<dim>::electrolyte_fe);
  battery<dim>::q_collection.push_back (battery<dim>::electrode_quadrature);
  battery<dim>::q_collection.push_back (battery<dim>::electrolyte_quadrature);
								 	
}


template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;