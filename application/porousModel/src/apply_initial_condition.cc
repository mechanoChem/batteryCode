#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::apply_initial_condition()
{ 

  double c_li_max_neg=battery<dim>::params->getDouble("c_li_max_neg");
	double c_li_max_pos=battery<dim>::params->getDouble("c_li_max_pos");
  
	double c_li_100_neg=battery<dim>::params->getDouble("c_li_100_neg");
  double c_li_100_pos=battery<dim>::params->getDouble("c_li_100_pos");
	double c_li_plus_ini=battery<dim>::params->getDouble("c_li_plus_ini");
	
	double T_0=battery<dim>::params->getDouble("T_0");
  double eps_l_0_sep=battery<dim>::params->getDouble("eps_l_0_sep");
  double eps_l_0_neg=battery<dim>::params->getDouble("eps_l_0_neg");
  double eps_l_0_pos=battery<dim>::params->getDouble("eps_l_0_pos");

	double dis_top0=battery<dim>::params->getDouble("dis_top0");
	
  battery<dim>::Un=0;
  //Un.compress(VectorOperation::add);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = battery<dim>::dof_handler.begin_active(), endc=battery<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
    hp::FEValues<dim> hp_fe_values (battery<dim>::fe_collection, battery<dim>::q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    FEFaceValues<dim> fe_face_values (*battery<dim>::electrode_fe, battery<dim>::common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
		std::vector<unsigned int> local_dof_indices (dofs_per_cell);
		unsigned int n_q_points= fe_values.n_quadrature_points;


    const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
    const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;
    
    cell->get_dof_indices (local_dof_indices);

      
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(face_center_upper[1]==battery<dim>::bY){
				fe_face_values.reinit (cell, 3);
				const unsigned int ckf = fe_face_values.get_fe().system_to_component_index(i).first;
				if(ckf==1) battery<dim>::Un(local_dof_indices[i])=dis_top0;
      }
 
			if(ck==4){	
				battery<dim>::Un(local_dof_indices[i])=c_li_plus_ini;//C_li_plus
    	}
			if(ck==5){	
				battery<dim>::Un(local_dof_indices[i])=T_0;//C_li_plus
    	}
 	 	 	if(ck==7){
 				electricChemoFormula->formula_Usc(c_li_100_neg,-1);
 				battery<dim>::Un(local_dof_indices[i])=-electricChemoFormula->Usc().val();//Phi_e
		 	}
			
      if(battery<dim>::cell_is_in_electrode_domain(cell) and face_center_upper[1]<=battery<dim>::electrode_Y1){
		  	if(ck==3){
					battery<dim>::Un(local_dof_indices[i])=c_li_100_neg*c_li_max_neg;//C_li
		  	}
	  	}
			
    	if(battery<dim>::cell_is_in_electrode_domain(cell) and face_center_low[1]>=battery<dim>::electrode_Y2){
		  	if(ck==3){
					battery<dim>::Un(local_dof_indices[i])=c_li_100_pos*c_li_max_pos;//C_li
		 	 }
				if(ck==6){
					electricChemoFormula->formula_Usc(c_li_100_pos,1);
					double Usc_pos=electricChemoFormula->Usc().val();
					electricChemoFormula->formula_Usc(c_li_100_neg,-1);
					double Usc_neg=electricChemoFormula->Usc().val();
					battery<dim>::Un(local_dof_indices[i])=Usc_pos-Usc_neg;//Phi_s
				}
      }
    } 
  }
}


template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;