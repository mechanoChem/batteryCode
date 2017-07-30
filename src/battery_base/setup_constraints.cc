#include"../../include/battery.h"

template <int dim>
void battery<dim>::setup_constraints()
{
  dof_handler.distribute_dofs (fe_collection);

  std::vector<bool> x_component (totalDOF, false); x_component[0]=true; 
  std::vector<bool> y_component (totalDOF, false); y_component[1]=true;  
  std::vector<bool> z_component (totalDOF, false); z_component[2]=true;
  std::vector<bool> c_li_component (totalDOF, false); c_li_component[3]=true;
  std::vector<bool> c_li_plus_component (totalDOF, false); c_li_plus_component[4]=true;
	std::vector<bool> T_component (totalDOF, false); T_component[5]=true;
  std::vector<bool> phi_s_component (totalDOF, false); phi_s_component[6]=true;
  std::vector<bool> phi_e_component (totalDOF, false); phi_e_component[7]=true;
  
  constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);

  VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints, x_component);
  VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints, y_component);
  VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints, z_component);
	
   VectorTools:: interpolate_boundary_values (dof_handler, 5, ZeroFunction<dim> (totalDOF),constraints, x_component);
  VectorTools:: interpolate_boundary_values (dof_handler, 5, ZeroFunction<dim> (totalDOF),constraints, y_component);
  VectorTools:: interpolate_boundary_values (dof_handler, 5, ZeroFunction<dim> (totalDOF),constraints, z_component);
	
	VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints,phi_s_component);
	VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints,phi_e_component);
  
  constraints.close ();
  
  std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl; 
}

template class battery<1>;
template class battery<2>;
template class battery<3>;