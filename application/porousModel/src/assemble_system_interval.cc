#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::assemble_system_interval(const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end)
{	
	double c_li_max_neg=battery<dim>::params->getDouble("c_li_max_neg");
	double c_li_max_pos=battery<dim>::params->getDouble("c_li_max_pos");
	
  hp::FEValues<dim> hp_fe_values (battery<dim>::fe_collection, battery<dim>::q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  FEFaceValues<dim> electrode_fe_face_values (*battery<dim>::electrode_fe, battery<dim>::common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*battery<dim>::electrolyte_fe, battery<dim>::common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 
  const unsigned int electrode_dofs_per_cell = battery<dim>::electrode_fe->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = battery<dim>::electrolyte_fe->dofs_per_cell;
  
  FullMatrix<double> local_matrix;
  Vector<double>            local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
		
  //loop over cells
  typename hp::DoFHandler<dim>::active_cell_iterator cell = battery<dim>::dof_handler.begin_active(), endc=battery<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    local_matrix.reinit (cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
    local_rhs.reinit (cell->get_fe().dofs_per_cell);
    local_dof_indices.resize (cell->get_fe().dofs_per_cell);
    cell->get_dof_indices (local_dof_indices);
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;

    //AD variables
    Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
		Table<1, double > ULocalConv(dofs_per_cell);
	  Table<1, double > U0Local(dofs_per_cell);
		
    for (unsigned int i=0; i<dofs_per_cell; ++i){
			if (std::abs(battery<dim>::U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
			else{ULocal[i]=battery<dim>::U(local_dof_indices[i]);}
			ULocal[i].diff (i, dofs_per_cell);
			ULocalConv[i]= battery<dim>::Un(local_dof_indices[i]);
			U0Local[i]= battery<dim>::U0(local_dof_indices[i]);
    }

    Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
		for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
		
 	 	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
 	 	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, 0, ULocal, defMap, battery<dim>::currentIteration);
 	 	deformationMap<double, dim> defMapConv(n_q_points); 
 	 	getDeformationMap<double, dim>(fe_values, 0, ULocalConv, defMapConv, battery<dim>::currentIteration);
		
		/*
		*evaluate primary fields
		*/
 	 	dealii::Table<1,double> c_li_conv(n_q_points), c_li_plus_conv(n_q_points),T_conv(n_q_points), phi_s_conv(n_q_points), phi_e_conv(n_q_points);
  	dealii::Table<1,Sacado::Fad::DFad<double> > c_li(n_q_points), c_li_plus(n_q_points), T(n_q_points), phi_s(n_q_points), phi_e(n_q_points), c_li_0(n_q_points);
  	dealii::Table<2,Sacado::Fad::DFad<double> > c_li_plus_grad(n_q_points, dim), T_grad(n_q_points, dim),phi_s_grad(n_q_points, dim),phi_e_grad(n_q_points, dim);
		
		evaluateScalarFunction<double,dim>(fe_values, 3, ULocalConv, c_li_conv);
		evaluateScalarFunction<double,dim>(fe_values, 4, ULocalConv, c_li_plus_conv);
		evaluateScalarFunction<double,dim>(fe_values, 5, ULocalConv, T_conv);
		evaluateScalarFunction<double,dim>(fe_values, 6, ULocalConv, phi_s_conv);
		evaluateScalarFunction<double,dim>(fe_values, 7, ULocalConv, phi_e_conv);
		
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 3, ULocal, c_li);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 4, ULocal, c_li_plus);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 5, ULocal, T);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 6, ULocal, phi_s);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 7, ULocal, phi_e);
		
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 4, ULocal, c_li_plus_grad, true, defMap);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 5, ULocal, T_grad, true, defMap);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 6, ULocal, phi_s_grad, true, defMap);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 7, ULocal, phi_e_grad, true, defMap);
			
		/*
		*update electrichemo formula
		*/
  	dealii::Table<1,Sacado::Fad::DFad<double> > jn(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > Q(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > eps_l(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > eps_s(n_q_points);	
    dealii::Table<1,Sacado::Fad::DFad<double> > D_l(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > Ke(n_q_points);
		
    dealii::Table<1,Sacado::Fad::DFad<double> > fac_chemo(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > fac_T(n_q_points);
		dealii::Table<1,double> eps_l_conv(n_q_points),eps_s_conv(n_q_points);
				
		
  	const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
  	const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;
		
		Sacado::Fad::DFad<double> c_li_max;
		double eps_s_0, eps_l_0;
		int domainflag;		
		//domain identification
 		if (battery<dim>::cell_is_in_electrode_domain(cell) and face_center_upper[1]<=battery<dim>::electrode_Y1){
			c_li_max=c_li_max_neg;
			domainflag=-1;
			eps_s_0=this->params->getDouble("eps_s_0_neg");
			eps_l_0=this->params->getDouble("eps_l_0_neg");
		}
		else if(battery<dim>::cell_is_in_electrode_domain(cell) and face_center_upper[1]>=battery<dim>::electrode_Y2){
			c_li_max=c_li_max_pos;
			domainflag=1;
			eps_s_0=this->params->getDouble("eps_s_0_pos");
			eps_l_0=this->params->getDouble("eps_l_0_pos");
		}
		else if(battery<dim>::cell_is_in_electrolyte_domain(cell)){
			c_li_max=c_li_max_pos;
			domainflag=0;
			eps_s_0=this->params->getDouble("eps_s_0_sep");
			eps_l_0=this->params->getDouble("eps_l_0_sep");
		}
		
		for (unsigned int q=0; q<n_q_points; ++q){
			Table<1,Sacado::Fad::DFad<double>>phi_s_grad_point(dim);
			Table<1,Sacado::Fad::DFad<double>>phi_e_grad_point(dim); 
			Table<1,Sacado::Fad::DFad<double>>c_li_plus_point(dim); 
			
			for (unsigned int i=0;i<dim;i++){
				phi_s_grad_point[i]=phi_s_grad[q][i];
				phi_e_grad_point[i]=phi_e_grad[q][i];
				c_li_plus_point[i]=c_li_plus_grad[q][i];
			}
			
			/*
			*direct call formula of porosity to get porosity at previous time step
			*/
			Sacado::Fad::DFad<double> UnitC_conv=c_li_conv[q]/c_li_max;
			Sacado::Fad::DFad<double> Temp_conv=T_conv[q];
			Sacado::Fad::DFad<double> defVol=defMapConv.detF[q];	
			electricChemoFormula->formula_porosity(UnitC_conv, Temp_conv, defVol,domainflag, currentflag);
			eps_s_conv[q]=electricChemoFormula->eps_s().val();
			eps_l_conv[q]=electricChemoFormula->eps_l().val();
			if(this->currentTime==0){eps_s_conv[q]=eps_s_0;eps_l_conv[q]=eps_l_0;}
			/*
			*update electricChemoFormula
			*/
			electricChemoFormula->update(c_li[q] , c_li_plus[q],c_li_plus_point,  T[q], phi_s[q], phi_s_grad_point, phi_e[q], phi_e_grad_point, defMap.detF[q], domainflag, currentflag);
			eps_s[q]=electricChemoFormula->eps_s();
			eps_l[q]=electricChemoFormula->eps_l();
			Ke[q]=electricChemoFormula->Ke();
			D_l[q]=electricChemoFormula->D_l();
			jn[q]=electricChemoFormula->jn()*eps_s[q];
			Q[q]=electricChemoFormula->Q_rxn()+electricChemoFormula->Q_rev()+electricChemoFormula->Q_ohm();
				
			fac_chemo[q]=electricChemoFormula->Beta()+1;
			fac_T[q]=std::pow(electricChemoFormula->Beta_T()+1,1.0/3.0);
		}	
		
		/*
		call residual
		*/	
		porousResidual->refresh(battery<dim>::currentIteration, battery<dim>::dt);
		porousResidual->IncludeStrain_chemo(fac_chemo);
		porousResidual->IncludeStrain_thermal(fac_T);
 	 	if(domainflag==1 or domainflag==-1){
   	 	if(dofs_per_cell!=electrode_dofs_per_cell){printf("dofs_per_cell!=electrode_dofs_per_cell"); exit(-1);}//check domain
   	 	porousResidual->residualForChemoCli(fe_values, dim, R, defMap, defMapConv, c_li, c_li_conv, jn, eps_s, eps_s_conv);
		 	porousResidual->residualForChemoCli_plus(fe_values, dim+1, R, defMap, defMapConv, c_li_plus, c_li_plus_conv, c_li_plus_grad, jn, eps_l, eps_l_conv, D_l);
			porousResidual->residualForChemoThermal(fe_values, dim+2, R, defMap, T, T_conv, T_grad, Q);
			porousResidual->residualForChemoPhis(fe_values, dim+3, R, defMap, phi_s, phi_s_grad, jn, eps_s);
			porousResidual->residualForChemoPhie(fe_values, dim+4, R, defMap, c_li_plus, c_li_plus_grad, T, phi_e, phi_e_grad, jn, eps_l, Ke);
			porousResidual->residualForMechanics(fe_values, 0,defMap, R);

		}   			
  	else if(domainflag==0){
    	if(dofs_per_cell!=electrolyte_dofs_per_cell){printf("dofs_per_cell!=electrolyte_dofs_per_cell"); exit(-1);}//check domain	
   	 	porousResidual->residualForChemoCli_plus(fe_values, dim+1, R, defMap, defMapConv, c_li_plus, c_li_plus_conv, c_li_plus_grad, jn, eps_l, eps_l_conv, D_l);
			porousResidual->residualForChemoThermal(fe_values, dim+2, R, defMap, T, T_conv, T_grad, Q);
   		porousResidual->residualForChemoPhie(fe_values, dim+4, R, defMap, c_li_plus, c_li_plus_grad, T, phi_e, phi_e_grad, jn, eps_l, Ke);
			porousResidual->residualForMechanics(fe_values, 0, defMap, R);
		}
		/*
		*boundary term for current input and heat convection
		*/
		for (unsigned int faceID=0; faceID<2*dim; faceID++){
	    if (cell->face(faceID)->at_boundary()){
			  FEFaceValues<dim>* fe_face_values;
				if(domainflag==1 or domainflag==-1) fe_face_values=&electrode_fe_face_values;
				else if(domainflag==0) fe_face_values=&electrolyte_fe_face_values;
				fe_face_values->reinit (cell, faceID);
				const unsigned int n_face_quadrature_points = fe_face_values->n_quadrature_points;
				dealii::Table<1,Sacado::Fad::DFad<double> > T_face(n_face_quadrature_points);
				for(unsigned int q=0;q<n_face_quadrature_points;q++){
					T_face[q]=0;
				  for (unsigned int i=0; i<dofs_per_cell; ++i) {
				    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - (dim+2);
						if(ck==0){
							T_face[q]+=fe_face_values->shape_value(i, q)*ULocal[i];
						}
					}
				}
				porousResidual->residualForHeatConv(fe_values, *fe_face_values, dim+2, R,T_face);
				/*
				*phi_e boundary condition for current input
				*/
	      const Point<dim> face_center = cell->face(faceID)->center();
				if (std::abs(face_center[1]-(battery<dim>::bY))<1.0e-3){
					porousResidual->residualForPhisBC(fe_values, *fe_face_values, dim+3, R, defMap, period, periodTime);
				}
			}
		}
		porousResidual->scalling(fe_values, 0, R, 1.0e-10);
		porousResidual->scalling(fe_values, 1, R, 1.0e-10);
		porousResidual->scalling(fe_values, 2, R, 1.0e-10);
		porousResidual->scalling(fe_values, dim, R, 1.0e1);
		porousResidual->scalling(fe_values, dim+1, R, 1.0e1);
		porousResidual->scalling(fe_values, dim+2, R, 1.0e-7);
		porousResidual->scalling(fe_values, dim+3, R, 1.0e-8);
		porousResidual->scalling(fe_values, dim+4, R, 1.0e-5);
		
    //Residual(R) and Jacobian(R')
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      for (unsigned int j=0; j<dofs_per_cell; ++j){
				// R' by AD
				local_matrix(i,j)= R[i].dx(j);
      }
      //R
      local_rhs(i) = -R[i].val(); 
    }
    //Global assemble with apply boundary condition and interface constraints
    battery<dim>::constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, battery<dim>::system_matrix, battery<dim>::system_rhs);
  }
}

template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;
