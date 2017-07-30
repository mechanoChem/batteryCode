#include"../../include/battery.h"

template <int dim>
void battery<dim>::assemble_system_interval(const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end)
{	
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);

  const unsigned int electrode_dofs_per_cell = electrode_fe->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = electrolyte_fe->dofs_per_cell;
  
  FullMatrix<double> local_matrix;
  
  Vector<double>            local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
		
  //loop over cells
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
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
			if (std::abs(U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
			else{ULocal[i]=U(local_dof_indices[i]);}
			ULocal[i].diff (i, dofs_per_cell);
			ULocalConv[i]= Un(local_dof_indices[i]);
			U0Local[i]= U0(local_dof_indices[i]);
    }

    Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
		for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
		
		/*
		*evaluate primary fields
		*/
 	 	dealii::Table<1,double> c_li_conv(n_q_points), c_li_plus_conv(n_q_points),T_conv(n_q_points), phi_s_conv(n_q_points), phi_e_conv(n_q_points);
  	dealii::Table<1,Sacado::Fad::DFad<double> > c_li(n_q_points), c_li_plus(n_q_points), T(n_q_points), phi_s(n_q_points), phi_e(n_q_points), c_li_0(n_q_points);
  	dealii::Table<2,Sacado::Fad::DFad<double> > c_li_plus_grad(n_q_points, dim), T_grad(n_q_points, dim),phi_s_grad(n_q_points, dim),phi_e_grad(n_q_points, dim);
		
		evaluateScalarFunction<double, dim>(fe_values, 3, ULocalConv, c_li_conv);
		evaluateScalarFunction<double, dim>(fe_values, 4, ULocalConv, c_li_plus_conv);
		evaluateScalarFunction<double, dim>(fe_values, 5, ULocalConv, T_conv);
		evaluateScalarFunction<double, dim>(fe_values, 6, ULocalConv, phi_s_conv);
		evaluateScalarFunction<double, dim>(fe_values, 7, ULocalConv, phi_e_conv);
		
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 3, ULocal, c_li);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 4, ULocal, c_li_plus);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 5, ULocal, T);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 6, ULocal, phi_s);
		evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 7, ULocal, phi_e);
		
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 4, ULocal, c_li_plus_grad);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 5, ULocal, T_grad);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 6, ULocal, phi_s_grad);
		evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 7, ULocal, phi_e_grad);
		
		
  	Sacado::Fad::DFad<double> UnitC;
		Sacado::Fad::DFad<double> dudt;
		Sacado::Fad::DFad<double> e;
		Sacado::Fad::DFad<double> a;
  	dealii::Table<1,Sacado::Fad::DFad<double> > jn(n_q_points), j0(n_q_points);
  	dealii::Table<1,Sacado::Fad::DFad<double> > eta(n_q_points), Usc(n_q_points);//,c_r(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > Q(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > R_e(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > eps_l(n_q_points);
    dealii::Table<1,Sacado::Fad::DFad<double> > eps_s(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > Q_ohm(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > Q_rev(n_q_points);
		dealii::Table<1,Sacado::Fad::DFad<double> > Q_rxn(n_q_points);
		
		dealii::Table<1,double> eps_l_conv(n_q_points),eps_s_conv(n_q_points);
		
 	 	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
 	 	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, 0, ULocal, defMap, currentIteration);
 	 	deformationMap<double, dim> defMapConv(n_q_points); 
 	 	getDeformationMap<double, dim>(fe_values, 0, ULocalConv, defMapConv, currentIteration);
		
		mechanicalResidual->residualForMechanics(fe_values, 0,defMap, R);

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
    constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
  }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;