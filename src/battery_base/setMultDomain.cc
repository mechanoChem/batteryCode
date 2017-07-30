#include"../../include/battery.h"

template <int dim>
void battery<dim>::setMultDomain()
{
  typename Triangulation<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (;cell!=endc; ++cell){
    const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
		const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;
    if(face_center_upper[1]<=electrode_Y1 or face_center_low[1]>=electrode_Y2){
      cell->set_material_id (electrode_domain_id);
    }
    else cell->set_material_id (electrolyte_domain_id);
  }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;