#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::setup_constraints()
{
  battery<dim>::dof_handler.distribute_dofs (battery<dim>::fe_collection);
	double totalDOF=battery<dim>::totalDOF;

  std::vector<bool> x_component (totalDOF, false); x_component[0]=true; 
  std::vector<bool> y_component (totalDOF, false); y_component[1]=true;  
  std::vector<bool> z_component (totalDOF, false); z_component[2]=true;
  std::vector<bool> c_li_component (totalDOF, false); c_li_component[3]=true;
  std::vector<bool> c_li_plus_component (totalDOF, false); c_li_plus_component[4]=true;
	std::vector<bool> T_component (totalDOF, false); T_component[5]=true;
  std::vector<bool> phi_s_component (totalDOF, false); phi_s_component[6]=true;
  std::vector<bool> phi_e_component (totalDOF, false); phi_e_component[7]=true;
  
  battery<dim>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (battery<dim>::dof_handler, battery<dim>::constraints);

  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, y_component);
  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, z_component);
	
  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 5, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 5, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, y_component);
  VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 5, ZeroFunction<dim> (totalDOF),battery<dim>::constraints, z_component);
	
	VectorTools:: interpolate_boundary_values (battery<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),battery<dim>::constraints,phi_s_component);

	//VectorTools:: interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim> (totalDOF),constraints,T_component);
	//VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints,T_component);

  battery<dim>::constraints.close ();
  
  std::cout << "   Number of active cells:       " << battery<dim>::triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << battery<dim>::dof_handler.n_dofs() << std::endl; 
}


template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;