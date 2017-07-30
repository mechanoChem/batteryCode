#include"../../include/battery.h"

template <int dim>
void battery<dim>::output_results (const unsigned int cycle) const
{
  std::cout<<"outPut"<<std::endl;
	std::vector<std::string> nodal_solution_names;
	std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
	
	//Nodal Solution names
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
	}
  nodal_solution_names.push_back("C_li"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("C_li_plus"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("T"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("phi_s"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
  nodal_solution_names.push_back("phi_e"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);  	
  //Write results to VTK file

  char filename1 [200]; sprintf (filename1, "output/output-%u.vtk", cycle); std::ofstream output1 (filename1);
  DataOut<dim,hp::DoFHandler<dim> >data_out; 
  data_out.attach_dof_handler (dof_handler);
  //Add nodal DOF data
  data_out.add_data_vector (Un, nodal_solution_names, DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}

template class battery<1>;
template class battery<2>;
template class battery<3>;