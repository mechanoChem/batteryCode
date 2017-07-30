#include"../../include/battery.h"

template <int dim>
void battery<dim>::apply_initial_condition()
{ 

  double c1max_neg=params->getDouble("csmax_neg");
	double c1max_pos=params->getDouble("csmax_pos");
  
  double cs_0_neg=params->getDouble("cs_0_neg");
	double cs_0_pos=params->getDouble("cs_0_pos");
	double cs_100_neg=params->getDouble("cs_100_neg");
        double cs_100_pos=params->getDouble("cs_100_pos");
	double c2_ini=params->getDouble("c2_ini");
	
	double T_0=params->getDouble("T_0");
        double eps_l_sep0=params->getDouble("eps_l_sep0");
        double eps_l_neg0=params->getDouble("eps_l_neg0");
        double eps_l_pos0=params->getDouble("eps_l_pos0");

	double dis_top0=params->getDouble("dis_top0");
	
  Un=0;
  //Un.compress(VectorOperation::add);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    FEFaceValues<dim> fe_face_values (*electrode_fe, common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
    //std::vector<types::global_dof_index> local_dof_indices;
    //local_dof_indices.resize (cell->get_fe().dofs_per_cell);
		
		std::vector<unsigned int> local_dof_indices (dofs_per_cell);
		unsigned int n_q_points= fe_values.n_quadrature_points;


    const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
    const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;
    
    cell->get_dof_indices (local_dof_indices);
   
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(face_center_upper[1]==bY){
				fe_face_values.reinit (cell, 3);
				const unsigned int ckf = fe_face_values.get_fe().system_to_component_index(i).first;
				if(ck==1) Un(local_dof_indices[i])=dis_top0;
			}

			//std::cout<<"A"<<std::endl;  
			if(ck==4){	
				Un(local_dof_indices[i])=c2_ini;//C_li_plus
    	}
			if(ck==5){	
				Un(local_dof_indices[i])=T_0;//C_li_plus
    	}
 	 	 if(ck==7){
			Un(local_dof_indices[i])=0;//Phi_e
			//Un(local_dof_indices[i])=-GetOpenCirculatePotential(cs_0_neg,-1).val();
		 }
			
      if(cell_is_in_electrode_domain(cell) and face_center_upper[1]<=electrode_Y1){
		  	if(ck==3){
					Un(local_dof_indices[i])=cs_100_neg*c1max_neg;//C_li
		  	}
	  	}
			
    	if(cell_is_in_electrode_domain(cell) and face_center_low[1]>=electrode_Y2){
		  	if(ck==3){
					Un(local_dof_indices[i])=cs_100_pos*c1max_pos;//C_li
		 	 }
				if(ck==6){
					Un(local_dof_indices[i])=0;//Phi_s
				}
      }
    } 
  }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;