#include"../../include/battery.h"

template <int dim>
bool battery<dim>::cell_is_in_electrode_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
    return (cell->material_id() == electrode_domain_id);
}
  
template <int dim>
bool battery<dim>::cell_is_in_electrolyte_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
    return (cell->material_id() == electrolyte_domain_id);
}

template <int dim>
void  battery<dim>::set_active_fe_indices ()
{
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end(); 
  for (;cell!=endc; ++cell){
    if (cell_is_in_electrode_domain (cell)) cell->set_active_fe_index (0);
    else if (cell_is_in_electrolyte_domain(cell)) cell->set_active_fe_index (1);
    else Assert (false, ExcNotImplemented());
  }
}

template <int dim>
void battery<dim>::mark_boundary()
{
  typename  Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  for (;cell!=endc; ++cell){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (cell->face(f)->at_boundary()){
				cell->face(f)->set_boundary_indicator (0);
				const Point<dim> face_center = cell->face(f)->center();
				if (face_center[0] == 0.0){
	  			cell->face(f)->set_boundary_indicator (1); //X-
				}
				if (face_center[1] == 0.0){
					cell->face(f)->set_boundary_indicator (2); //Y-
			 	}
			 	if (face_center[2] == 0.0){
					cell->face(f)->set_boundary_indicator (3); //Z-
			 	}
			 	if (face_center[0] == bX){
					cell->face(f)->set_boundary_indicator (4); //X+
				}
			 	if (face_center[1] == bY){
				 cell->face(f)->set_boundary_indicator (5); //Y+
			 	}
			 	if (face_center[2] == bZ){
					cell->face(f)->set_boundary_indicator (6); //Z+
			 	}
  		}
    }
  }
}

template <int dim>
void battery<dim>::setup_system()
{
  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  constraints.condense (csp);
  sparsity_pattern.copy_from (csp);
	
  system_matrix.reinit (sparsity_pattern);
  U.reinit (dof_handler.n_dofs()); 
  Un.reinit (dof_handler.n_dofs());  
  dU.reinit (dof_handler.n_dofs()); 
  system_rhs.reinit (dof_handler.n_dofs()); 
  U0.reinit (dof_handler.n_dofs());
}

template <int dim>
void battery<dim>::assemble_system()
{
  std::cout<<"begin assemble_system"<<std::endl;
  system_matrix=0; system_rhs=0; //boundary_values.clear();
	
  const unsigned int n_threads=dealii::MultithreadInfo::n_threads();
  Threads::ThreadGroup<> threads;
  typedef typename hp::DoFHandler<dim>::active_cell_iterator active_cell_iterator;
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (), dof_handler.end (), n_threads);
  for (unsigned int thread=0; thread<n_threads; ++thread){
    threads += Threads::new_thread (&battery<dim>::assemble_system_interval, *this, thread_ranges[thread].first, thread_ranges[thread].second);
    threads.join_all ();
  }
  std::cout<<"finish assemble_system"<<std::endl;
}

template class battery<1>;
template class battery<2>;
template class battery<3>;