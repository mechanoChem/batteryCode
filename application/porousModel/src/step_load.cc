#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::step_load()
{
  double dis_top0=battery<dim>::params->getDouble("dis_top0");
  typename hp::DoFHandler<dim>::active_cell_iterator cell = battery<dim>::dof_handler.begin_active(), endc=battery<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
    FEFaceValues<dim> fe_face_values (*battery<dim>::electrode_fe, battery<dim>::common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    cell->get_dof_indices (local_dof_indices);
    const Point<dim> face_center_upper = cell->face(3)->center();
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      //const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(face_center_upper[1]==battery<dim>::bY){
        fe_face_values.reinit (cell, 3);
        const unsigned int ckf = fe_face_values.get_fe().system_to_component_index(i).first;
        if(ckf==1) {
	  			if(battery<dim>::currentTime>=0.1) battery<dim>::U(local_dof_indices[i])=-0.5;
	  			if(battery<dim>::currentTime>=0.2) battery<dim>::U(local_dof_indices[i])=-0.8;
	  			if(battery<dim>::currentTime>=0.3) battery<dim>::U(local_dof_indices[i])=-1;
	  			if(battery<dim>::currentTime>=0.4) battery<dim>::U(local_dof_indices[i])=-1.2;
	  			if(battery<dim>::currentTime>=0.5) battery<dim>::U(local_dof_indices[i])=-1.4;
          if(battery<dim>::currentTime>=0.6) battery<dim>::U(local_dof_indices[i])=-1.6;
          if(battery<dim>::currentTime>=0.7) battery<dim>::U(local_dof_indices[i])=-1.8;
          if(battery<dim>::currentTime>=0.8) battery<dim>::U(local_dof_indices[i])=-2;
	  			if(battery<dim>::currentTime>=0.9) battery<dim>::U(local_dof_indices[i])=-2.2;
	  			// if(currentTime>=0.8) U(local_dof_indices[i])=-2.4;
	  			//if(currentTime>=0.3) IpA=44;
				}
      }
    }
  }
}

template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;