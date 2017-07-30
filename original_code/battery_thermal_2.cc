#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <Sacado.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>
#include <tbb/task_scheduler_init.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include "chemo_thermal.h"
#include "mechanics.h"
#include "supplementaryFunctions.h"
#include "quadraturePointData.h"
#include "parameters.h"
using namespace dealii;

#define DIMS 3
//#define problemWidth
#define totalDOF 8

template <int dim>
class battery
{
  public:
    battery (const unsigned int li_degree, const unsigned int phi_degree,const unsigned int elasticity_degree, parametersClass* _params);
    ~battery();
    void run ();
		parametersClass* params;

    enum
    {
      electrode_domain_id,
      electrolyte_domain_id,
    };

  private:

    static bool
    cell_is_in_electrode_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);
    static bool
    cell_is_in_electrolyte_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell); 
    
    void make_grid();
    void mark_boundary();
    void set_active_fe_indices ();
    void setup_dofs();
    void setup_system();
    void apply_initial_condition();
    //void apply_phi_e_BC();
    void assemble_system ();
    void assemble_system_interval (const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end);
    void solve();
    void output_results (const unsigned int cycle) const;
    void apply_Dirichlet_BC();
    void Global_index_map();
    void step_load();
	//void reset_potential();(code is deleted see back up)
		
		//quadratic projection
		void setup_system_projection_electrode();
		void projection_quadrature_electrode();
		
		void setup_system_projection_whole();
		void projection_quadrature_whole();
		
		void projection_output_results_1 (const unsigned int cycle) const;
		void projection_output_results_2 (const unsigned int cycle) const;
		void projection_output_results_3 (const unsigned int cycle) const;
                void projection_output_results_4 (const unsigned int cycle) const;
                void projection_output_results_5 (const unsigned int cycle) const;
		
		dealii::Table<2,int >   Global_index;

    const unsigned int li_degree;
    const unsigned int phi_degree;
    const unsigned int elasticity_degree;
    

    
    Triangulation<dim>    triangulation;
    //FESystem<dim>*         electrode_fe;
    //FESystem<dim>*         electrolyte_fe;
std::shared_ptr<FESystem<dim>> electrode_fe;
std::shared_ptr<FESystem<dim>> electrolyte_fe;

    hp::FECollection<dim> fe_collection;
    hp::DoFHandler<dim>   dof_handler;
    hp::QCollection<dim>  q_collection;
    ConstraintMatrix      constraints;
    
		
    SparsityPattern       sparsity_pattern;
		
    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs, U, Un, dU, U0;
		
		
    const QGauss<dim> electrode_quadrature;
    const QGauss<dim> electrolyte_quadrature;
    const QGauss<dim-1> common_face_quadrature;
    
    
    //geometry information
    double bX,bY,bZ;
  double electrode_Y1,electrode_Y2;

    //solution variables
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
    unsigned int currentIncrement, currentIteration,period;
    double totalTime, currentTime, periodTime, dt;
    //Threads::ThreadMutex assembler_lock;
	
	
	double Potential_bot;
	double Potential_top;

	double U_bot;
	double U_top;
		
		//projection
		//projection variable for variable only lives in electrode
    FESystem<dim>         projection_electrode_fe_electrode;
    FESystem<dim>         projection_electrolyte_fe_electrode;
		
	  hp::FECollection<dim> projection_fe_collection_electrode;
	  hp::DoFHandler<dim>   projection_dof_handler_electrode;
	  ConstraintMatrix      projection_constraints_electrode;   
	  SparsityPattern       projection_sparsity_pattern_electrode;
	  SparseMatrix<double>  projection_mass_matrix_electrode;

          std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > C_surface;
		 
		Vector<double>        projection_system_rhs_1;
		Vector<double> 				projection_nodalValues_1;		
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_1;
		
		Vector<double>        projection_system_rhs_2;
		Vector<double> 				projection_nodalValues_2;		
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_2;
		
		Vector<double>        projection_system_rhs_3;
		Vector<double> 				projection_nodalValues_3;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_3;
		
		//projection variable for whole
    FESystem<dim>         projection_electrode_fe_whole;
    FESystem<dim>         projection_electrolyte_fe_whole;
		
	  hp::FECollection<dim> projection_fe_collection_whole;
	  hp::DoFHandler<dim>   projection_dof_handler_whole;
	  ConstraintMatrix      projection_constraints_whole;   
	  SparsityPattern       projection_sparsity_pattern_whole;
	  SparseMatrix<double>  projection_mass_matrix_whole; 
		
		Vector<double>        projection_system_rhs_4;
		Vector<double> 				projection_nodalValues_4;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_4;
		
  	Vector<double>        projection_system_rhs_5;
		Vector<double> 				projection_nodalValues_5;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_5;
		
  	Vector<double>        projection_system_rhs_6;
		Vector<double> 				projection_nodalValues_6;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, std::vector<double> > projection_quadratureData_6;
		
	//int cycle=0;
};

template <int dim>
battery<dim>::battery (const unsigned int li_degree, const unsigned int phi_degree,const unsigned int elasticity_degree, parametersClass* _params)
  :
  li_degree(li_degree), phi_degree(phi_degree), elasticity_degree(elasticity_degree),params(_params),
	projection_electrode_fe_electrode(FE_Q<dim>(li_degree),1),
	projection_electrolyte_fe_electrode(FE_Nothing<dim>(),1),
	projection_electrode_fe_whole(FE_Q<dim>(li_degree),1),
	projection_electrolyte_fe_whole(FE_Q<dim>(li_degree),1),
  electrode_quadrature(li_degree+2), electrolyte_quadrature(li_degree+2), common_face_quadrature(phi_degree+2),
  dof_handler (triangulation),
	projection_dof_handler_electrode(triangulation),projection_dof_handler_whole(triangulation)
{

  electrode_fe.reset(new FESystem<dim>(FE_Q<dim>(elasticity_degree),dim,
	             FE_Q<dim>(li_degree),1,
	       			 FE_Q<dim>(li_degree),2,
	       			 FE_Q<dim>(phi_degree),1,
	       			 FE_Q<dim>(phi_degree),1));
  electrolyte_fe.reset(new FESystem<dim>(FE_Q<dim>(elasticity_degree),dim,
		 						 FE_Nothing<dim>(),1,
		 						 FE_Q<dim>(li_degree),2,
		 						 FE_Nothing<dim>(),1,
		 						 FE_Q<dim>(phi_degree),1));
	
  fe_collection.push_back (*electrode_fe);
  fe_collection.push_back (*electrolyte_fe);
	projection_fe_collection_electrode.push_back (projection_electrode_fe_electrode);
	projection_fe_collection_electrode.push_back (projection_electrolyte_fe_electrode);
	projection_fe_collection_whole.push_back (projection_electrode_fe_whole);
	projection_fe_collection_whole.push_back (projection_electrolyte_fe_whole);
  q_collection.push_back (electrode_quadrature);
  q_collection.push_back (electrolyte_quadrature);
  //solution variables
	totalTime=params->getDouble("totalTime");
	dt=params->getDouble("dt");
  currentIncrement=0; currentTime=0;
  period=0;periodTime=0;
  //Nodal Solution names
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
	}
  nodal_solution_names.push_back("C_li"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("C_li_plus"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("T"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  nodal_solution_names.push_back("phi_s"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
  nodal_solution_names.push_back("phi_e"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);  	
}

template <int dim>
battery<dim>::~battery (){dof_handler.clear ();}

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
void battery<dim>::make_grid()
{
	//geometry
	double l_neg=params->getDouble("l_neg");//120.0 negtive electrode
	double l_s=params->getDouble("l_s");//23.0
	double l_pos=params->getDouble("l_pos");//92.0 positive electrode
	double w_b1=params->getDouble("w_b1");//1e3.0
	double w_b2=params->getDouble("w_b2");//1e3.0
	
  //subdivided_hyper
  bool colorize = false;
  std::vector< std::vector< double > > step_sizes;
  step_sizes.resize(3);
  step_sizes[0].push_back(w_b1);
  step_sizes[2].push_back(w_b2);
  //for(unsigned int cell=1;cell<=10;cell++) step_sizes[0].push_back(2);
  //for(unsigned int cell=1;cell<=10;cell++) step_sizes[2].push_back(2);
	
  for(unsigned int cell=1;cell<=1;cell++) step_sizes[1].push_back(60);
  for(unsigned int cell=1;cell<=1;cell++) step_sizes[1].push_back(23);
  for(unsigned int cell=1;cell<=1;cell++) step_sizes[1].push_back(45);
	
	/*
	step_sizes[1].push_back(tab_neg);
	step_sizes[1].push_back(l_neg);
	step_sizes[1].push_back(l_s);
	step_sizes[1].push_back(l_pos);
	step_sizes[1].push_back(tab_pos);
	*/
	
  bX=w_b1; bY=l_neg+l_s+l_pos; bZ=w_b2;
  electrode_Y1=l_neg;
	electrode_Y2=l_neg+l_s;

  GridGenerator::subdivided_hyper_rectangle (triangulation, step_sizes,
                                             Point<3>(0.0,0.0,0.0),
                                             Point<3>(bX,bY,bZ),
					     colorize);
  //assign material_id to each cell
  typename Triangulation<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (;cell!=endc; ++cell){
    const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
		const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;
    if(face_center_upper[1]<=l_neg or face_center_low[1]>=l_neg+l_s){
			cell->set_material_id (electrode_domain_id);
		}
		else {
			cell->set_material_id (electrolyte_domain_id);
		}
  }
	
  typename Triangulation<dim>::active_cell_iterator projection_cell = projection_dof_handler_electrode.begin_active(), projection_endc = projection_dof_handler_electrode.end();
  for (;projection_cell!=projection_endc; ++projection_cell){
    const Point<dim> face_center_low = projection_cell->face(2)->center();//face n=y-;
		const Point<dim> face_center_upper = projection_cell->face(3)->center();//face n=y+;
    if(face_center_upper[1]<=electrode_Y1 or face_center_low[1]>=electrode_Y2){
      projection_cell->set_material_id (electrode_domain_id);
    }
    else projection_cell->set_material_id (electrolyte_domain_id);
  }
	
  typename Triangulation<dim>::active_cell_iterator projection_cell1 = projection_dof_handler_whole.begin_active(), projection_endc1 = projection_dof_handler_whole.end();
  for (;projection_cell1!=projection_endc1; ++projection_cell1){
    const Point<dim> face_center_low = projection_cell1->face(2)->center();//face n=y-;
		const Point<dim> face_center_upper = projection_cell1->face(3)->center();//face n=y+;
    if(face_center_upper[1]<=electrode_Y1 or face_center_low[1]>=electrode_Y2){
      projection_cell1->set_material_id (electrode_domain_id);
    }
    else projection_cell1->set_material_id (electrolyte_domain_id);
  }
	
}

template <int dim>
void battery<dim>::mark_boundary()
{
	double l_neg=params->getDouble("l_neg");//120.0 negtive electrode
	double l_s=params->getDouble("l_s");//23.0
	double l_pos=params->getDouble("l_pos");//92.0 positive electrode
	double w_b1=params->getDouble("w_b1");
	double w_b2=params->getDouble("w_b2");
  typename  Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  for (;cell!=endc; ++cell){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (cell->face(f)->at_boundary()){
				cell->face(f)->set_boundary_indicator (0);
				const Point<dim> face_center = cell->face(f)->center();
				if (face_center[0] == 0.0){
	  			cell->face(f)->set_boundary_indicator (1); //X
				}
				if (face_center[1] == 0.0){
	  			cell->face(f)->set_boundary_indicator (2); //Y
				}
				if (face_center[2] == 0.0){
	  			cell->face(f)->set_boundary_indicator (3); //Z
				}
				if (std::abs(face_center[0] - bX)<1e-3){
	  			cell->face(f)->set_boundary_indicator (4); //X
				}
				if (std::abs(face_center[1] - bY)<1e-3){
	  			cell->face(f)->set_boundary_indicator (5); //Y
				}
				if (std::abs(face_center[2] - bZ)<1e-3){
	  			cell->face(f)->set_boundary_indicator (6); //Z
				}
      }
    }
  }
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
	
  typename hp::DoFHandler<dim>::active_cell_iterator projection_cell = projection_dof_handler_electrode.begin_active(), projection_endc=projection_dof_handler_electrode.end(); 
  for (;projection_cell!=projection_endc; ++projection_cell){
    if (cell_is_in_electrode_domain (projection_cell)) projection_cell->set_active_fe_index (0);
    else if (cell_is_in_electrolyte_domain(projection_cell)) projection_cell->set_active_fe_index (1);
    else Assert (false, ExcNotImplemented());
  }
	
  typename hp::DoFHandler<dim>::active_cell_iterator projection_cell1 = projection_dof_handler_whole.begin_active(), projection_endc1=projection_dof_handler_whole.end(); 
  for (;projection_cell1!=projection_endc1; ++projection_cell1){
    if (cell_is_in_electrode_domain (projection_cell1)) projection_cell1->set_active_fe_index (0);
    else if (cell_is_in_electrolyte_domain(projection_cell1)) projection_cell1->set_active_fe_index (1);
    else Assert (false, ExcNotImplemented());
  }
	
}

template <int dim>
void battery<dim>::setup_dofs()
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
  
	//VectorTools:: interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim> (totalDOF),constraints,T_component);
	//VectorTools:: interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim> (totalDOF),constraints,T_component);
	
  //interface constraints
  std::vector<types::global_dof_index> local_face_dof_indices (electrode_fe->dofs_per_face);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    if(cell_is_in_electrode_domain(cell)){
      for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
		  if(!cell->at_boundary(f)){
			  bool face_is_on_interface = false;
			  if((cell->neighbor(f)->has_children() == false) and (cell_is_in_electrolyte_domain (cell->neighbor(f)))){
				  face_is_on_interface = true;
			  }
			  else if (cell->neighbor(f)->has_children() == true){printf("Did mesh refiments"); exit(-1);}
			  if (face_is_on_interface){
				  cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
				  for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
					  //if (electrode_fe.face_system_to_component_index(i).first ==3) constraints.add_line (local_face_dof_indices[i]);//c_li
					  //if (electrode_fe.face_system_to_component_index(i).first ==5) constraints.add_line (local_face_dof_indices[i]);//phi_s
				  }    
			  }
		  }
     }
    }
		
    if(cell_is_in_electrolyte_domain(cell)){
      for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
		  if(!cell->at_boundary(f)){
			  bool face_is_on_interface = false;
			  if((cell->neighbor(f)->has_children() == false) and (cell_is_in_electrode_domain (cell->neighbor(f)))){
				  face_is_on_interface = true;
			  }
			  else if (cell->neighbor(f)->has_children() == true){printf("Did mesh refiments"); exit(-1);}
			  if (face_is_on_interface){
				  cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
				  for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
					  //if (electrode_fe.face_system_to_component_index(i).first ==3) constraints.add_line (local_face_dof_indices[i]);//c_li
					  //if (electrode_fe.face_system_to_component_index(i).first ==6) constraints.add_line (local_face_dof_indices[i]);//phi_e
				  }    
			  }
		  }
      }
    }
  }
  constraints.close ();
  
  std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl; 
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
void battery<dim>::setup_system_projection_electrode()
{
	projection_dof_handler_electrode.distribute_dofs (projection_fe_collection_electrode);	
	projection_constraints_electrode.clear ();
	DoFTools::make_hanging_node_constraints (projection_dof_handler_electrode, projection_constraints_electrode);
	projection_constraints_electrode.close ();
	CompressedSimpleSparsityPattern projection_csp_electrode (projection_dof_handler_electrode.n_dofs(), projection_dof_handler_electrode.n_dofs());

	DoFTools::make_sparsity_pattern (projection_dof_handler_electrode, projection_csp_electrode, projection_constraints_electrode, false);
  projection_constraints_electrode.condense (projection_csp_electrode);
  projection_sparsity_pattern_electrode.copy_from (projection_csp_electrode);
  projection_mass_matrix_electrode.reinit(projection_sparsity_pattern_electrode);
  MatrixCreator::create_mass_matrix(projection_dof_handler_electrode, q_collection, projection_mass_matrix_electrode);
  projection_constraints_electrode.condense (projection_mass_matrix_electrode);
	
	//projection_system_rhs_JN.reinit (projection_dof_handler_electrode.n_dofs()); 
	//projection_nodalValues_JN.reinit (projection_dof_handler_electrode.n_dofs());
	
	//projection_system_rhs_Re.reinit (projection_dof_handler_electrode.n_dofs()); 
	//projection_nodalValues_Re.reinit (projection_dof_handler_electrode.n_dofs());
	
	std::cout << "   Number of degrees of freedom of JN projection: " << projection_dof_handler_electrode.n_dofs() << std::endl;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    unsigned int n_q_points= fe_values.n_quadrature_points;
		
		//projection quadratureData resize.
    std::vector<double> temp1(n_q_points, 0.0); 
    //projection_quadratureData_JN.insert(std::make_pair(cell, temp1));
		//projection_quadratureData_Re.insert(std::make_pair(cell, temp1));
	}
}

template <int dim>
void battery<dim>::setup_system_projection_whole()
{
	projection_dof_handler_whole.distribute_dofs (projection_fe_collection_whole);	
	projection_constraints_whole.clear ();
	DoFTools::make_hanging_node_constraints (projection_dof_handler_whole, projection_constraints_whole);
	projection_constraints_whole.close ();	
	CompressedSimpleSparsityPattern projection_csp_whole (projection_dof_handler_whole.n_dofs(), projection_dof_handler_whole.n_dofs());
	DoFTools::make_sparsity_pattern (projection_dof_handler_whole, projection_csp_whole, projection_constraints_whole, false);
  projection_constraints_whole.condense (projection_csp_whole);
  projection_sparsity_pattern_whole.copy_from (projection_csp_whole);
	
  projection_mass_matrix_whole.reinit(projection_sparsity_pattern_whole);
  MatrixCreator::create_mass_matrix(projection_dof_handler_whole, q_collection, projection_mass_matrix_whole);
  projection_constraints_whole.condense (projection_mass_matrix_whole);
	
	projection_system_rhs_1.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_1.reinit (projection_dof_handler_whole.n_dofs());
	projection_system_rhs_2.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_2.reinit (projection_dof_handler_whole.n_dofs());
	projection_system_rhs_3.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_3.reinit (projection_dof_handler_whole.n_dofs());
	projection_system_rhs_4.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_4.reinit (projection_dof_handler_whole.n_dofs());
	projection_system_rhs_5.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_5.reinit (projection_dof_handler_whole.n_dofs());
	projection_system_rhs_6.reinit (projection_dof_handler_whole.n_dofs());
	projection_nodalValues_6.reinit (projection_dof_handler_whole.n_dofs());

	std::cout << "   Number of degrees of freedom of whole projection: " << projection_dof_handler_whole.n_dofs() << std::endl;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    unsigned int n_q_points= fe_values.n_quadrature_points;
		
		//projection quadratureData resize.
    std::vector<double> temp1(n_q_points, 0.0); 
		projection_quadratureData_1.insert(std::make_pair(cell, temp1));
		projection_quadratureData_2.insert(std::make_pair(cell, temp1));	
		projection_quadratureData_3.insert(std::make_pair(cell, temp1));
		projection_quadratureData_4.insert(std::make_pair(cell, temp1));	
		projection_quadratureData_5.insert(std::make_pair(cell, temp1));
		projection_quadratureData_6.insert(std::make_pair(cell, temp1));
		C_surface.insert(std::make_pair(cell, temp1));
	}
}

template <int dim>
void battery<dim>::assemble_system()
{	
	
	double Rr=params->getDouble("Rr");
	double F=params->getDouble("F");
	double T_0=params->getDouble("T_0");
	double IpA=params->getDouble("IpA");

	double alpha_neg=params->getDouble("alpha_neg");
	double alpha_pos=params->getDouble("alpha_pos");
	double k_neg=params->getDouble("k_neg");
	double k_pos=params->getDouble("k_pos");

	double c1max_neg=params->getDouble("csmax_neg");
	double c1max_pos=params->getDouble("csmax_pos");
	
  double cs_100_neg=params->getDouble("cs_100_neg");
	double cs_0_neg=params->getDouble("cs_0_neg");
	double cs_100_pos=params->getDouble("cs_100_pos");
	double cs_0_pos=params->getDouble("cs_0_pos");
	
	//R_s Radius of solid particals(m);esp_s volume fraction active material
	double R_s_neg=params->getDouble("R_s_neg");
	double R_s_pos=params->getDouble("R_s_pos");
	
	double eps_s_neg0=params->getDouble("eps_s_neg0");
	double eps_s_pos0=params->getDouble("eps_s_pos0");
	double eps_s_sep0=params->getDouble("eps_s_sep0");
	double eps_l_sep0=params->getDouble("eps_l_sep0");
	double eps_l_neg0=params->getDouble("eps_l_neg0");
	double eps_l_pos0=params->getDouble("eps_l_pos0");
	
	double kappa_neg=params->getDouble("kappa_neg");
	double kappa_pos=params->getDouble("kappa_pos");
	double kappa_l=params->getDouble("kappa_l");
	double kappa_s=params->getDouble("kappa_s");
	double D_s_neg=params->getDouble("D_s_neg");
	double D_s_pos=params->getDouble("D_s_pos");
	
	double se_neg=params->getDouble("se_neg");
	double se_pos=params->getDouble("se_pos");
	double t_0=params->getDouble("t_0");
	
	
	system_matrix=0; system_rhs=0;
	
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);

  FEFaceValues<dim> electrode_fe_face_values (*electrode_fe, common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*electrolyte_fe, common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 

  const unsigned int electrode_dofs_per_cell = electrode_fe->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = electrolyte_fe->dofs_per_cell;
  
  FullMatrix<double> local_matrix;
  //FullMatrix<double> local_interface_matrix (electrode_dofs_per_cell, electrolyte_dofs_per_cell);
  Vector<double>            local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  //std::vector<types::global_dof_index> neighbor_dof_indices (electrolyte_dofs_per_cell);
	
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
    	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell); Table<1, double > U0Local(dofs_per_cell);
			
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
			if (std::abs(U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
			else{ULocal[i]=U(local_dof_indices[i]);}
			ULocal[i].diff (i, dofs_per_cell);
			ULocalConv[i]= Un(local_dof_indices[i]);
			U0Local[i]= U0(local_dof_indices[i]);
	    }
			
   	 	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
   	 	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, 0, ULocal, defMap, currentIteration);
   	 	deformationMap<double, dim> defMapConv(n_q_points); 
   	 	getDeformationMap<double, dim>(fe_values, 0, ULocalConv, defMapConv, currentIteration);
			
    	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
			for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
   
   	 	dealii::Table<1,double> c_li_conv(n_q_points), c_li_plus_conv(n_q_points),T_conv(n_q_points), phi_s_conv(n_q_points), phi_e_conv(n_q_points);
    	dealii::Table<1,Sacado::Fad::DFad<double> > c_li(n_q_points), c_li_plus(n_q_points), T(n_q_points), phi_s(n_q_points), phi_e(n_q_points), c_li_0(n_q_points);
    	dealii::Table<2,Sacado::Fad::DFad<double> > c_li_plus_grad(n_q_points, dim), T_grad(n_q_points, dim),phi_s_grad(n_q_points, dim),phi_e_grad(n_q_points, dim);
    	dealii::Table<1,Sacado::Fad::DFad<double> > jn(n_q_points), j0(n_q_points);
    	dealii::Table<1,Sacado::Fad::DFad<double> > eta(n_q_points), Usc(n_q_points);//,c_r(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > Q(n_q_points);
								
    	for (unsigned int q=0; q<n_q_points; ++q){
      	c_li_conv[q]=0.0; c_li_plus_conv[q]=0.0; T_conv[q]=0; phi_s_conv[q]=0.0; phi_e_conv[q]=0.0;
      	c_li[q]=0.0; c_li_plus[q]=0.0; T[q]=0; phi_s[q]=0.0; phi_e[q]=0.0; c_li_0[q]=0.0;
			
				//reset projection value
				projection_quadratureData_1[cell][q]=0.0;
				projection_quadratureData_2[cell][q]=0.0;
				projection_quadratureData_3[cell][q]=0.0;
				//projection_quadratureData_4[cell][q]=0.0;
				//projection_quadratureData_5[cell][q]=0.0;	
				//projection_quadratureData_6[cell][q]=0.0;

      	for (unsigned int j=0; j<dim; j++) {
					c_li_plus_grad[q][j]=0.0; T_grad[q][j]=0.0; phi_s_grad[q][j]=0;phi_e_grad[q][j]=0;
				}
				for (unsigned int i=0; i<dofs_per_cell; ++i) {
					const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - dim;
					if (ck==0) {
						c_li[q]+=fe_values.shape_value(i, q)*ULocal[i];
						c_li_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; 
         		c_li_0[q]+=fe_values.shape_value(i, q)*U0Local[i];
     			}
					if (ck==1) {
						c_li_plus[q]+=fe_values.shape_value(i, q)*ULocal[i]; c_li_plus_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
		  			for (unsigned int j=0; j<dim; j++) {
							c_li_plus_grad[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
		  			}
      		}
					if (ck==2) {
						T[q]+=fe_values.shape_value(i, q)*ULocal[i]; T_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
		  			for (unsigned int j=0; j<dim; j++) {
							T_grad[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
		  			}
      		}
     			if (ck==3) {
						phi_s[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_s_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
		  			for (unsigned int j=0; j<dim; j++) {
							phi_s_grad[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
		  			}
     			}
     			if (ck==4) {
						phi_e[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_e_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
		  			for (unsigned int j=0; j<dim; j++) {
							phi_e_grad[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
		  			}
     			}
   			}    		 
    	}

    	Sacado::Fad::DFad<double> UnitC;
			Sacado::Fad::DFad<double> e;
			Sacado::Fad::DFad<double> a;
			dealii::Table<1,Sacado::Fad::DFad<double> > R_e(n_q_points);
	    dealii::Table<1,Sacado::Fad::DFad<double> > eps_l(n_q_points);
	    dealii::Table<1,Sacado::Fad::DFad<double> > eps_s(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > Q_ohm(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > Q_rev(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > Q_rxn(n_q_points);
			
			dealii::Table<1,double> eps_l_conv(n_q_points),eps_s_conv(n_q_points);
			
			Sacado::Fad::DFad<double> dudt;
			
    	const Point<dim> face_center_low = cell->face(2)->center();//face n=y-;
    	const Point<dim> face_center_upper = cell->face(3)->center();//face n=y+;

    	for (unsigned int q=0; q<n_q_points; ++q){
      	// -
      	Usc[q]=0.0; j0[q]=0.0; jn[q]=0.0; eta[q] = 0-phi_e[q];
				eps_s[q]=1-eps_l_sep0;
				eps_l[q]=eps_l_sep0;
				R_e[q]=0.0;
				Q_rxn[q]=0;
				Q_rev[q]=0;
				Q_ohm[q]=0;
				Q[q]=0;

				
				//negtive electrode
     		if (cell_is_in_electrode_domain(cell) and face_center_upper[1]<=electrode_Y1){	
				  double Pl=0;
					double Pb=0;	
					double Rs=R_s_neg;
					
					double UnitC_conv;
					double beta_conv, beta_s_conv,beta_T_conv,beta_s_T_conv;
					Sacado::Fad::DFad<double> UnitC;
					Sacado::Fad::DFad<double> UnitC_surface;
					Sacado::Fad::DFad<double> beta;
					Sacado::Fad::DFad<double> beta_s;
					Sacado::Fad::DFad<double> beta_T;
					Sacado::Fad::DFad<double> beta_s_T;

					Sacado::Fad::DFad<double> factor;
					beta_s=0;
					beta=0;
					beta_T=0;
					beta_s_T=0;
					beta_conv=0;
					beta_s_conv=0;
					beta_T_conv=0;
					beta_s_T_conv=0;
					

					beta_T_conv=(T_conv[q]-298)*9.615e-6;
					beta_s_T_conv=(T_conv[q]-298)*6e-6;
					beta_T=(T[q]-298)*9.615e-6;
					beta_s_T=(T[q]-298)*6e-6;
					
				  UnitC_conv=c_li_conv[q]/c1max_neg;  
				   beta_s_conv=GetBeta_s(UnitC_conv,1).val();//1 discharging				
				  beta_conv=GetBeta(UnitC_conv,1).val();
				  if(currentTime==0){beta_s_conv=0;beta_conv=0;}
					//eps_s_conv[q] =(1+(kappa_neg*(defMapConv.detF[q]-1-beta_conv)-Pl-Pb)/kappa_s+beta_s_conv)/defMapConv.detF[q]*eps_s_neg0;
					//eps_l_conv[q]=1-eps_s_conv[q]-(1-eps_s_neg0-eps_l_neg0)/defMapConv.detF[q];
				  eps_s_conv[q]=((kappa_neg*(defMapConv.detF[q]/(1+beta_conv)/(1+beta_T_conv)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s_conv)*(1+beta_s_T_conv)/defMapConv.detF[q]*eps_s_neg0;
		eps_l_conv[q]=1-eps_s_conv[q]-(1-eps_s_neg0-eps_l_neg0)/defMapConv.detF[q];
	       
					UnitC=c_li[q]/c1max_neg;
					  beta_s=GetBeta_s(UnitC,1);//1 discharging 
					  beta=GetBeta(UnitC,1);
					// eps_s[q] =(1+(kappa_neg*(defMap.detF[q]-1-beta-beta_T)-Pl-Pb)/kappa_s+beta_s+beta_s_T)/defMap.detF[q]*eps_s_neg0; //1.0-0.3-0.044;
					eps_s[q]=((kappa_neg*(defMap.detF[q]/(1+beta)/(1+beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_s_T)/defMap.detF[q]*eps_s_neg0;
					eps_l[q]=1-eps_s[q]-(1-eps_s_neg0-eps_l_neg0)/defMap.detF[q];
					
					//a_p
					factor=((kappa_neg*(defMap.detF[q]/(1+beta)/(1+beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_s_T);
					R_e[q]=std::pow(factor,1/3)*R_s_neg;
                                        a=3.0/R_e[q];

					UnitC_surface=C_surface[cell][q]/c1max_neg;
					j0[q] = eps_s[q]*a*k_neg*std::pow(1000,0.5)*c1max_neg*std::pow(c_li_plus[q]*1.0e3,alpha_neg)*std::pow((1-UnitC_surface),alpha_neg)*std::pow(UnitC_surface,(1-alpha_neg));

					Usc[q]=GetOpenCirculatePotential(UnitC,-1);
					eta[q] = phi_s[q]-phi_e[q]-Usc[q];
					jn[q] = j0[q]*(exp(alpha_neg*F/Rr/T[q]*(eta[q]))-exp(-alpha_neg*F/Rr/T[q]*(eta[q])));
				       
					C_surface[cell][q]=c_li[q].val()-R_e[q].val()*jn[q].val()/5/D_s_neg;
					//thermal terms Q
					Q_rxn[q]=F*jn[q]*eta[q];
					dudt=GetQrev(UnitC);//discharging
					Q_rev[q]=F*a*jn[q]*T[q]*dudt;
					
					Sacado::Fad::DFad<double> Ke;
					Ke=((34.5*exp(-798/T[q])*std::pow((c_li_plus[q]*1.0e3),3)-485*exp(-1080/T[q])*std::pow((c_li_plus[q]*1.0e3),2)+2440*exp(-1440/T[q])*(c_li_plus[q]*1.0e3))/10)*1.0e6;
			    for (unsigned int j = 0; j < dim; j++){
			    	Q_ohm[q]+=se_neg*eps_s[q]*phi_s_grad[q][j]*phi_s_grad[q][j]+(Ke*std::pow(eps_l[q], 1.5)*(phi_e_grad[q][j]-2.0*Rr*T[q]/F*(1.0-t_0)/c_li_plus[q]*(c_li_plus_grad[q][j])))*phi_e_grad[q][j];
			    }
			    //Q[q]=0;
			    Q[q]=Q_rxn[q]+Q_rev[q]+Q_ohm[q];
					
				}
				//positive electrode
      	else if (cell_is_in_electrode_domain(cell) and face_center_low[1]>=electrode_Y2){			
				  double Pl=0;
					double Pb=0;	
					double Rs=R_s_pos;
					Sacado::Fad::DFad<double> UnitC_conv;
					double beta_conv, beta_s_conv,beta_T_conv,beta_s_T_conv;
					Sacado::Fad::DFad<double> UnitC;
					Sacado::Fad::DFad<double> UnitC_surface;
					Sacado::Fad::DFad<double> beta;
					Sacado::Fad::DFad<double> beta_s;
					Sacado::Fad::DFad<double> beta_T;
					Sacado::Fad::DFad<double> beta_s_T;
					Sacado::Fad::DFad<double> factor;
					beta_s=0;
					beta=0;
					beta_T=0;
					beta_s_T=0;
					beta_conv=0;
					beta_s_conv=0;
					beta_T_conv=0;
					beta_s_T_conv=0;
					
				      
					beta_T_conv=(T_conv[q]-298)*6.025e-6;
                                        beta_s_T_conv=(T_conv[q]-298)*6e-6;
                                        beta_T=(T[q]-298)*6.025e-6; 
                                        beta_s_T=(T[q]-298)*6e-6; 
					
					//eps_s_conv[q] =(1+(kappa_pos*(defMapConv.detF[q]-1-beta_conv)-Pl-Pb)/kappa_s+beta_s_conv)/defMapConv.detF[q]*eps_s_pos0;
		eps_s_conv[q]=((kappa_pos*(defMapConv.detF[q]/(1+beta_conv)/(1+beta_T_conv)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s_conv)*(1+beta_s_T_conv)/defMapConv.detF[q]*eps_s_pos0;
					eps_l_conv[q]=1-eps_s_conv[q]-(1-eps_s_pos0-eps_l_pos0)/defMapConv.detF[q];
					
					UnitC=c_li[q]/c1max_pos;
					//eps_s[q] =(1+(kappa_pos*(defMap.detF[q]-1-beta-beta_T)-Pl-Pb)/kappa_s+beta_s+beta_s_T)/defMap.detF[q]*eps_s_pos0; //1.0-0.3-0.044;
					eps_s[q]=((kappa_pos*(defMap.detF[q]/(1+beta)/(1+beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_s_T)/defMap.detF[q]*eps_s_pos0;
					eps_l[q]=1-eps_s[q]-(1-eps_s_pos0-eps_l_pos0)/defMap.detF[q];

                                        factor=((kappa_pos*(defMap.detF[q]/(1+beta)/(1+beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_s_T);
                                        R_e[q]=std::pow(factor,1/3)*R_s_pos;
                                        a=3.0/R_e[q];

					UnitC_surface=C_surface[cell][q]/c1max_pos;
					j0[q] = eps_s[q]*a*k_pos*std::pow(1000,0.5)*c1max_pos*std::pow(c_li_plus[q]*1.0e3,alpha_pos)*std::pow((1-UnitC_surface),alpha_pos)*std::pow(UnitC_surface,(1-alpha_pos));

   
					Usc[q]=GetOpenCirculatePotential(UnitC,1);
		  		eta[q] = phi_s[q]-phi_e[q]-Usc[q];
          //Pore-wall flux across interface (mol m-2 s-1)
          jn[q] = j0[q]*(exp(alpha_pos*F/Rr/T[q]*(eta[q]))-exp(-alpha_pos*F/Rr/T[q]*(eta[q])));
	  C_surface[cell][q]=c_li[q].val()-R_e[q].val()*jn[q].val()/5/D_s_pos;
					//thermal terms Q
					Q_rxn[q]=F*jn[q]*eta[q];
					dudt=GetQrev(UnitC);
					Q_rev[q]=F*a*jn[q]*T[q]*dudt;
					
					Sacado::Fad::DFad<double> Ke;
					Ke=((34.5*exp(-798/T[q])*std::pow((c_li_plus[q]*1.0e3),3)-485*exp(-1080/T[q])*std::pow((c_li_plus[q]*1.0e3),2)+2440*exp(-1440/T[q])*(c_li_plus[q]*1.0e3))/10)*1.0e6;
			    for (unsigned int j = 0; j < dim; j++){
			    	Q_ohm[q]+=se_pos*eps_s[q]*phi_s_grad[q][j]*phi_s_grad[q][j]+(Ke*std::pow(eps_l[q], 1.5)*(phi_e_grad[q][j]-2.0*Rr*T[q]/F*(1.0-t_0)/c_li_plus[q]*(c_li_plus_grad[q][j])))*phi_e_grad[q][j];
			    }
			    //Q[q]=0;
			    Q[q]=Q_rxn[q]+Q_rev[q]+Q_ohm[q];
     	 	}
				//sepeartor
				else if(cell_is_in_electrolyte_domain(cell)){
				  double Pl=0;
					double Pb=0;	
					Sacado::Fad::DFad<double> UnitC_conv;
					double beta_conv, beta_s_conv, beta_T_conv, beta_s_T_conv;
					Sacado::Fad::DFad<double> UnitC;
					Sacado::Fad::DFad<double> beta;
					Sacado::Fad::DFad<double> beta_s;
					Sacado::Fad::DFad<double> beta_T;
					Sacado::Fad::DFad<double> beta_s_T;
					
					beta_s=0;
					beta=0;
					beta_T=0;
					beta_s_T=0;
					beta_conv=0;
					beta_s_conv=0;
					beta_T_conv=0;
					beta_s_T_conv=0;
					
					beta_T_conv=(T_conv[q]-298)*82.46e-6;
                                        beta_s_T_conv=(T_conv[q]-298)*6e-6;
                                        beta_T=(T[q]-298)*82.46e-6;
                                        beta_s_T=(T[q]-298)*6e-6;

					//eps_s_conv[q] =(1+(kappa_l*(defMapConv.detF[q]-1-beta_conv)-Pl-Pb)/kappa_s+beta_s_conv)/defMapConv.detF[q]*eps_s_sep0;
		 eps_s_conv[q]=((kappa_l*(defMapConv.detF[q]/(1+beta_conv)/(1+beta_T_conv)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s_conv)*(1+beta_s_T_conv)/defMapConv.detF[q]*eps_s_sep0;
					eps_l_conv[q]=1-eps_s_conv[q]-(1-eps_s_sep0-eps_l_sep0)/defMapConv.detF[q];
					
					//eps_s[q] =(1+(kappa_l*(defMap.detF[q]-1-beta-beta_T)-Pl-Pb)/kappa_s+beta_s+beta_s_T)/defMap.detF[q]*eps_s_sep0;
					eps_s[q]=((kappa_l*(defMap.detF[q]/(1+beta)/(1+beta_T)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_s_T)/defMap.detF[q]*eps_s_sep0;
					eps_l[q]=1-eps_s[q]-(1-eps_s_sep0-eps_l_sep0)/defMap.detF[q];
					
					Sacado::Fad::DFad<double> Ke;
					Ke=((34.5*exp(-798/T[q])*std::pow((c_li_plus[q]*1.0e3),3)-485*exp(-1080/T[q])*std::pow((c_li_plus[q]*1.0e3),2)+2440*exp(-1440/T[q])*(c_li_plus[q]*1.0e3))/10)*1.0e6;
			    for (unsigned int j = 0; j < dim; j++){
			    	Q_ohm[q]+=(Ke*std::pow(eps_l[q], 1.5)*(phi_e_grad[q][j]-2.0*Rr*T[q]/F*(1.0-t_0)/c_li_plus[q]*(c_li_plus_grad[q][j])))*phi_e_grad[q][j];
			    }
			    //Q[q]=0;
			    Q[q]=Q_ohm[q];
     	 	}
		//Q[q]=0;
				projection_quadratureData_1[cell][q]=eps_s[q].val();
				projection_quadratureData_2[cell][q]=eps_l[q].val();
				projection_quadratureData_3[cell][q]=jn[q].val();
				//projection_quadratureData_Eps_l[cell][q]=jn[q].val(); 
    	}
	//std::cout<<"jn=="<<jn[0].val()<<std::endl;  
			//Q_rxn one order bigger than Q_ohm
			//std::cout<<"Q_rxn"<<Q_rxn[0].val()<<std::endl;
			//std::cout<<"Q_ohm"<<Q_ohm[0].val()<<std::endl;
	//std::cout<<"eps_s_conv[q]="<<eps_s_conv[0]<<"eps_l_conv[q]="<<eps_l_conv[0]<<std::endl;
	//std::cout<<"eps_s_conv[q]="<<eps_s[0].val()<<"eps_l_conv[q]="<<eps_l[0].val()<<std::endl;			
    	//chemo
			
			//std::cout<<"jn="<<jn[0].val()<<std::endl;
			//std::cout<<"eps_s="<<eps_s[0].val()<<std::endl;
			//std::cout<<"eps_l="<<eps_l[0].val()<<std::endl;
			
			Table<1, double > Pyy(n_q_points);
   	 	if(cell_is_in_electrode_domain(cell)){
     	 	if(dofs_per_cell!=electrode_dofs_per_cell){printf("dofs_per_cell!=electrode_dofs_per_cell"); exit(-1);}//check domain
				std::cout<<"jn="<<jn[0].val()<<std::endl;
				std::cout<<"eps_s="<<eps_s[0].val()<<std::endl;
				std::cout<<"eps_s_conv="<<eps_l[0].val()<<std::endl;
				std::cout<<"c_li="<<c_li[0].val()<<std::endl;
				std::cout<<"c_li_conv="<<c_li_conv[0]<<std::endl;
     	 	residualForChemo1(fe_values, dim+0, dt, currentTime, totalTime, R, defMap, defMapConv, c_li, c_li_conv, jn, eps_s, eps_s_conv, params);
				residualForChemo2(fe_values, dim+1, dt, currentTime, totalTime, R, defMap, defMapConv,  c_li_plus, c_li_plus_conv, c_li_plus_grad, T,jn, eps_l, eps_l_conv, params);
				residualForChemoThermal(fe_values, dim+2, electrode_fe_face_values,ULocal, cell, dt, currentTime, totalTime, R, defMap, T_conv, T, T_grad, Q, params);
     		residualForChemo3(fe_values, dim+3, electrode_fe_face_values, cell, dt, currentTime, totalTime,periodTime, period, R, defMap, phi_s, phi_s_grad, jn, eps_s, params);
     		residualForChemo4(fe_values, dim+4, electrode_fe_face_values, dt, currentTime, totalTime, R, defMap, c_li_plus, c_li_plus_grad, T, phi_e, phi_e_grad, jn, eps_l, params);
		double youngsModulus = 0.5e-3; // Pa/um^2
     		residualForMechanics(Pyy, fe_values, 0, electrode_fe_face_values, ULocal, ULocalConv, R, defMap, youngsModulus, cell,c_li_conv, c_li, c_li_0, T, currentTime, totalTime, dt);

			}   			
    	if(cell_is_in_electrolyte_domain(cell)){
      	if(dofs_per_cell!=electrolyte_dofs_per_cell){printf("dofs_per_cell!=electrolyte_dofs_per_cell"); exit(-1);}//check domain	
     	 	residualForChemo2(fe_values, dim+1, dt, currentTime, totalTime, R, defMap, defMapConv,  c_li_plus, c_li_plus_conv, c_li_plus_grad, T,jn, eps_l, eps_l_conv, params);
		residualForChemoThermal(fe_values, dim+2, electrolyte_fe_face_values,ULocal, cell, dt, currentTime, totalTime, R, defMap, T_conv, T, T_grad, Q, params);
     		residualForChemo4(fe_values, dim+4, electrolyte_fe_face_values, dt, currentTime, totalTime, R, defMap, c_li_plus, c_li_plus_grad, T, phi_e, phi_e_grad, jn, eps_l, params);
		double youngsModulus = 0.5e-3; // Pa/um^2
		residualForMechanics(Pyy, fe_values, 0, electrolyte_fe_face_values, ULocal, ULocalConv, R, defMap, youngsModulus, cell,c_li_conv, c_li, c_li_0, T, currentTime, totalTime, dt);
			}
			
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				std::cout<<"i="<<i<<"  "<<R[i].val()<<std::endl;
			}
			exit(-1);
	for(unsigned int q=0;q<n_q_points;q++){
	   projection_quadratureData_5[cell][q]=Pyy[q];
	  // projection_quadratureData_3[cell][q]=1;	
	}
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
     	 	for (unsigned int j=0; j<dofs_per_cell; ++j){
					// R' by AD
					local_matrix(i,j)= R[i].dx(j);
		 			 if (std::abs(local_matrix(i,j))<1.0e-16) local_matrix(i,j)=0.0;	 
      	}
      	//R
      	local_rhs(i) = -R[i].val();  
	  		if (std::abs(local_rhs(i))<1.0e-16) local_rhs(i)=0.0;
			}
    	//Global assemble with apply boundary condition and interface constraints
   	 constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
		 
		}
	 	for (unsigned int i=0; i<dof_handler.n_dofs();i++){
	 		//if(Global_index[i][0]==7)std::cout<<"system_rhs"<<system_rhs[i]<<" i="<<i<<" "<<Global_index[i][0]<<" "<<Global_index[i][1]<<std::endl;
	 	}
}


template <int dim>
void battery<dim>::solve(){
  double res=1, tol=1.0e-12, abs_tol=1.0e-13, initial_norm=0, current_norm=0;
  double machineEPS=1.0e-15;
  currentIteration=0;


	double iteration_n=params->getDouble("iteration_n");
	double iteration_ini=params->getDouble("iteration_ini");
  while (true){
    if (currentIteration>=iteration_ini) {printf ("Maximum number of iterations reached without convergence. \n"); break; exit (1);}
    if(currentTime>0 and currentIteration>=iteration_n) {printf ("Maximum number of iterations reached without convergence. \n"); break; exit (1);}
    //if (current_norm>1/std::pow(tol,2)){printf ("\nNorm is too high. \n\n"); break; exit (1);}
    assemble_system();
    current_norm=system_rhs.l2_norm();

    initial_norm=std::max(initial_norm, current_norm);
    res=current_norm/initial_norm;
    printf ("Inc:%3u (time:%10.3e, dt:%10.3e), Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); 
    if (res<tol || current_norm< abs_tol){printf ("Residual converged in %u iterations.\n\n", currentIteration); break;}

    //Direct solver	
    std::cout<<"begin solves"<<std::endl;
    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (dU, system_rhs);
									
		constraints.distribute (dU);
    U+=dU;
    ++currentIteration;
    std::cout << std::flush;
    std::cout<<"finish solves"<<std::endl;


  }
  //apply_phi_e_BC();
  Un=U;
}

template <int dim>
void battery<dim>::output_results (const unsigned int cycle) const
{
  std::cout<<"outPut"<<std::endl;
  //Write results to VTK file
  //dof_handler (triangulation);
  char filename1 [200]; sprintf (filename1, "output/output-%u.vtk", cycle); std::ofstream output1 (filename1);
  DataOut<dim,hp::DoFHandler<dim> >data_out; 
  data_out.attach_dof_handler (dof_handler);
  //Add nodal DOF data
  data_out.add_data_vector (Un, nodal_solution_names, DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
  if (cycle>=1) {
		projection_output_results_1 (cycle);
		projection_output_results_2 (cycle);
		projection_output_results_3 (cycle);
		projection_output_results_4 (cycle);
		projection_output_results_5(cycle);
	}
}

template <int dim>
void battery<dim>::run()
{
  make_grid();
  mark_boundary();
  set_active_fe_indices ();
  setup_dofs();
  setup_system();
	Global_index_map();
	//setup_system_projection_electrode();
	setup_system_projection_whole();
  apply_initial_condition();
  U=Un;U0=Un;
  //apply_phi_e_BC();
  output_results(0); //output initial state

  currentIncrement=0;
  periodTime=0;
  for (currentTime=0; currentTime<=totalTime; currentTime+=dt){
    periodTime+=dt;
    currentIncrement++;
    //if(currentTime>1) dt=1;
    if(currentIncrement>358 and period <1) {
      period++;
      //std::cout<<"period_ini="<<period<<std::endl;                                                                                                                               
      dt=0.1;
      periodTime=0;
    }
    if(periodTime>1) dt=1;
    //if(currentTime>1) dt=10;
		//if(currentTime>=2.393e+04) dt=10;
		//if(currentTime>=2.394e+04) dt=40;
		//if(currentTime>=1000)  dt=100;
		//if(currentTime>3000)  dt=1;
		clock_t t_solve;	
  	t_solve = clock();
		solve();
	  t_solve = clock() - t_solve;
		printf("************");
	  printf ("It took me %d clicks (%f seconds) for one solve.\n ",t_solve,((float)t_solve)/CLOCKS_PER_SEC);
	
	//projection_quadrature_electrode();
	   projection_quadrature_whole();
	   // if(currentIncrement%100==0){
	   // projection_quadrature_whole();
	   //output_results(currentIncrement); //output solution at current load increment
	   //}
	 	output_results(currentIncrement);	
		if(currentTime<=1) step_load();
  }
	
	//delete electrode_fe;
	//delete electrolyte_fe;
}

template <int dim>
void battery<dim>::Global_index_map()
{ 
	Global_index.reinit(dof_handler.n_dofs(),2);
	std::cout<<"!!";
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    cell->get_dof_indices (local_dof_indices);
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	int ck = fe_values.get_fe().system_to_component_index(i).first;
				Global_index[local_dof_indices[i]][0]=ck;
				//ck_file(local_dof_indices[i])=ck;
				if(cell_is_in_electrode_domain(cell))	Global_index[local_dof_indices[i]][1]=0;
				else if(cell_is_in_electrolyte_domain(cell)) Global_index[local_dof_indices[i]][1]=1;
			}
		}
}

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
 
    for (unsigned int q=0; q<n_q_points; ++q){
   
      C_surface[cell][q]=0.0;
      if(cell_is_in_electrode_domain(cell) and face_center_upper[1]<=electrode_Y1) C_surface[cell][q]=cs_100_neg*c1max_neg;
      if(cell_is_in_electrode_domain(cell) and face_center_low[1]>=electrode_Y2) C_surface[cell][q]=cs_100_pos*c1max_pos;
    }
   
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(face_center_upper[1]==bY){
	fe_face_values.reinit (cell, 3);
	const unsigned int ckf = fe_face_values.get_fe().system_to_component_index(i).first;
	if(ckf==1) Un(local_dof_indices[i])=dis_top0;

      }

			//std::cout<<"A"<<std::endl;  
			if(ck==4){	
				Un(local_dof_indices[i])=c2_ini;//C_li_plus
    	}
			if(ck==5){	
				Un(local_dof_indices[i])=T_0;//C_li_plus
    	}
 	 	 if(ck==7){
			Un(local_dof_indices[i])=-GetOpenCirculatePotential(cs_100_neg,-1).val();//Phi_e
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
					Un(local_dof_indices[i])=GetOpenCirculatePotential(cs_100_pos,1).val()-GetOpenCirculatePotential(cs_100_neg,-1).val();//Phi_s
				}
      }
    } 
  }
}

template <int dim>
void battery<dim>::step_load()
{
  double dis_top0=params->getDouble("dis_top0");
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc=dof_handler.end();
  for (;cell!=endc; ++cell){
    FEFaceValues<dim> fe_face_values (*electrode_fe, common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    cell->get_dof_indices (local_dof_indices);
    const Point<dim> face_center_upper = cell->face(3)->center();
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      //const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(face_center_upper[1]==bY){
        fe_face_values.reinit (cell, 3);
        const unsigned int ckf = fe_face_values.get_fe().system_to_component_index(i).first;
        if(ckf==1) {
	  if(currentTime>=0.1) U(local_dof_indices[i])=-0.5;
	  if(currentTime>=0.2) U(local_dof_indices[i])=-0.8;
	  if(currentTime>=0.3) U(local_dof_indices[i])=-1;
	  if(currentTime>=0.4) U(local_dof_indices[i])=-1.2;
	  if(currentTime>=0.5) U(local_dof_indices[i])=-1.4;
          if(currentTime>=0.6) U(local_dof_indices[i])=-1.6;
          if(currentTime>=0.7) U(local_dof_indices[i])=-1.8;
          if(currentTime>=0.8) U(local_dof_indices[i])=-2;
	  if(currentTime>=0.9) U(local_dof_indices[i])=-2.2;
	  // if(currentTime>=0.8) U(local_dof_indices[i])=-2.4;
	  //if(currentTime>=0.3) IpA=44;
	}
      }
    }
  }
}

template <int dim>
void battery<dim>::projection_quadrature_electrode()
{

}

template <int dim>
void battery<dim>::projection_quadrature_whole()
{
	std::cout<<"begin solving quadrature_projection whole"<<std::endl;
	projection_system_rhs_1=0;
	projection_system_rhs_2=0;
	projection_system_rhs_3=0;
	projection_system_rhs_4=0;
	projection_system_rhs_5=0;
	//projection_system_rhs_6=0;
	Vector<double>       local_rhs_1;
	Vector<double>       local_rhs_2;
	Vector<double>       local_rhs_3;
	Vector<double>       local_rhs_4;
	Vector<double>       local_rhs_5;
	//Vector<double>       local_rhs_6;
	
	  std::vector<unsigned int> local_dof_indices;
	  hp::FEValues<dim> hp_fe_values (projection_fe_collection_whole, q_collection, update_values   | update_gradients | update_quadrature_points | update_JxW_values);
	  typename hp::DoFHandler<dim>::active_cell_iterator cell = projection_dof_handler_whole.begin_active(), endc = projection_dof_handler_whole.end();
	  for (;cell!=endc; ++cell){
	    hp_fe_values.reinit (cell);
	    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();    
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			unsigned int num_quad_points= fe_values.n_quadrature_points;
			local_dof_indices.resize (cell->get_fe().dofs_per_cell);
	    cell->get_dof_indices (local_dof_indices);

			local_rhs_1.reinit (dofs_per_cell);
			local_rhs_2.reinit (dofs_per_cell);
			local_rhs_3.reinit (dofs_per_cell);
			local_rhs_4.reinit (dofs_per_cell);
			local_rhs_5.reinit (dofs_per_cell);
			//local_rhs_6.reinit (dofs_per_cell);
			
			local_dof_indices.resize (dofs_per_cell);
			local_rhs_1=0;
			local_rhs_2=0;
			local_rhs_3=0;
			local_rhs_4=0;
			//local_rhs_5=0;
			//local_rhs_6=0;
	    //Looping over quadrature points	
	    for (unsigned int q=0; q<num_quad_points; ++q){
				for (unsigned int i=0; i<dofs_per_cell; ++i){
					local_rhs_1(i) += (fe_values.shape_value(i,q)*projection_quadratureData_1[cell][q])*fe_values.JxW(q);
					local_rhs_2(i) += (fe_values.shape_value(i,q)*projection_quadratureData_2[cell][q])*fe_values.JxW(q);
					local_rhs_3(i) += (fe_values.shape_value(i,q)*projection_quadratureData_3[cell][q])*fe_values.JxW(q);
					local_rhs_4(i) += (fe_values.shape_value(i,q)*C_surface[cell][q])*fe_values.JxW(q);
					local_rhs_5(i) += (fe_values.shape_value(i,q)*projection_quadratureData_5[cell][q])*fe_values.JxW(q);
					//local_rhs_6(i) += (fe_values.shape_value(i,q)*projection_quadratureData_1[cell][q])*fe_values.JxW(q);
				}
	    }
			projection_constraints_whole.distribute_local_to_global (local_rhs_1, local_dof_indices, projection_system_rhs_1);
			projection_constraints_whole.distribute_local_to_global (local_rhs_2, local_dof_indices, projection_system_rhs_2);
			projection_constraints_whole.distribute_local_to_global (local_rhs_3, local_dof_indices, projection_system_rhs_3);
			projection_constraints_whole.distribute_local_to_global (local_rhs_4, local_dof_indices, projection_system_rhs_4);
			projection_constraints_whole.distribute_local_to_global (local_rhs_5, local_dof_indices, projection_system_rhs_5);
			//projection_constraints_whole.distribute_local_to_global (local_rhs_6, local_dof_indices, projection_system_rhs_6);
		}

    SparseDirectUMFPACK  B_direct;
    B_direct.initialize(projection_mass_matrix_whole);	
		
		B_direct.vmult (projection_nodalValues_1, projection_system_rhs_1);
	  projection_constraints_whole.distribute (projection_nodalValues_1);
	  B_direct.vmult (projection_nodalValues_2, projection_system_rhs_2);
	  projection_constraints_whole.distribute (projection_nodalValues_2);
	  B_direct.vmult (projection_nodalValues_3, projection_system_rhs_3);
	  projection_constraints_whole.distribute (projection_nodalValues_3);
	  B_direct.vmult (projection_nodalValues_4, projection_system_rhs_4);
	  projection_constraints_whole.distribute (projection_nodalValues_4);
          B_direct.vmult (projection_nodalValues_5, projection_system_rhs_5);
          projection_constraints_whole.distribute (projection_nodalValues_5);
		std::cout<<"finished solving quadrature_projection whole"<<std::endl;
}


template <int dim>
void battery<dim>::projection_output_results_1 (const unsigned int cycle) const
{
  char filename [200]; sprintf (filename, "output/eps_s_Projections-%u.vtk", cycle); std::ofstream output (filename);
  DataOut<dim,hp::DoFHandler<dim> > data_out; data_out.attach_dof_handler (projection_dof_handler_whole);
  //Add nodal DOF data
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
  data_out.add_data_vector (projection_nodalValues_1, "eps_s" , DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches(); data_out.write_vtk (output); output.close();
}
template <int dim>
void battery<dim>::projection_output_results_2 (const unsigned int cycle) const
{
  char filename [200]; sprintf (filename, "output/eps_l_Projections-%u.vtk", cycle); std::ofstream output (filename);
  DataOut<dim,hp::DoFHandler<dim> > data_out; data_out.attach_dof_handler (projection_dof_handler_whole);
  //Add nodal DOF data
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
  data_out.add_data_vector (projection_nodalValues_2, "eps_l" , DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches(); data_out.write_vtk (output); output.close();
}

template <int dim>
void battery<dim>::projection_output_results_3 (const unsigned int cycle) const
{
  char filename [200]; sprintf (filename, "output/jn_Projections-%u.vtk", cycle); std::ofstream output (filename);
  DataOut<dim,hp::DoFHandler<dim> > data_out; data_out.attach_dof_handler (projection_dof_handler_whole);
  //Add nodal DOF data
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
  data_out.add_data_vector (projection_nodalValues_3, "jn" , DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches(); data_out.write_vtk (output); output.close();
}

template <int dim>
void battery<dim>::projection_output_results_4 (const unsigned int cycle) const
{
  char filename [200]; sprintf (filename, "output/C_surface_Projections-%u.vtk", cycle); std::ofstream output(filename);
  DataOut<dim,hp::DoFHandler<dim> > data_out; data_out.attach_dof_handler (projection_dof_handler_whole);
  //Add nodal DOF data
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
  data_out.add_data_vector (projection_nodalValues_4, "C_surface" , DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches(); data_out.write_vtk (output); output.close();
}

template <int dim>
void battery<dim>::projection_output_results_5 (const unsigned int cycle) const
{
  char filename [200]; sprintf (filename, "output/Pyy_Projections-%u.vtk", cycle); std::ofstream output(filename);
  DataOut<dim,hp::DoFHandler<dim> > data_out; data_out.attach_dof_handler (projection_dof_handler_whole);
  //Add nodal DOF data                                                                                             
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
  data_out.add_data_vector (projection_nodalValues_5, "Pyy" , DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches(); data_out.write_vtk (output); output.close();
}
int main (){
  try{
	
    deallog.depth_console (0);
		//set param values;;;unit:pmol; um;
		parametersClass params;
		//constants
		params.setDouble("F", 96485.3329);
		params.setDouble("Rr", 8.3144598);
		params.setDouble("T_0", 298.0);

		//geometry
		params.setDouble("l_neg",60.0);
		params.setDouble("l_s",23.0);
		params.setDouble("l_pos",45.0);
		
		params.setDouble("w_b1",120e3);
		params.setDouble("w_b2",85e3);
		//params.setDouble("w_b1",10);
		//params.setDouble("w_b2",10);
		// particle
		params.setDouble("alpha_neg",0.5);
		params.setDouble("alpha_pos",0.5);
		//reaction rate
		params.setDouble("k_neg",8.0e-4);
		params.setDouble("k_pos",8.0e-4);
		//params.setDouble("k_pos",1);
    //R_s Radius of solid particals(m);esp_s,esp_l volume fraction active material/electrolyte
		//params.setDouble("eps_s_neg0",0.63);
		params.setDouble("eps_s_neg0",0.53);
		params.setDouble("eps_s_pos0",0.5);
		params.setDouble("eps_s_sep0",0.35);
		params.setDouble("eps_l_sep0",0.65);
		params.setDouble("eps_l_neg0",0.32);
		params.setDouble("eps_l_pos0",0.35);
		//a_s_neg=3.01;a_s_pos=0.753;
		params.setDouble("R_s_neg",8.0);//
		params.setDouble("R_s_pos",6.0);//
		//solid diffusion
		params.setDouble("D_s_neg",0.5);
		params.setDouble("D_s_pos",0.1);
		//se, silid electronic conductivity
		params.setDouble("se_neg",1.5e8);
		params.setDouble("se_pos",0.5e8);
		
		params.setDouble("kappa_l",0.42e-3);//0.4GPA
		params.setDouble("kappa_neg",4.94e-3);//3Gpa
		params.setDouble("kappa_pos",7.4e-3);//4Gpa
		params.setDouble("kappa_s",25e-3);
		
		params.setDouble("density_neg",2.5e-15);//kg/um^3
		params.setDouble("density_sep",1.1e-15);//
		params.setDouble("density_pos",2.5e-15);//
		params.setDouble("Cp",7e14);//pJ/kgK
		params.setDouble("lambda_neg",1.04e6);//kg/um^3
		params.setDouble("lambda_sep",0.33e6);//
		params.setDouble("lambda_pos",5e6);//

		params.setDouble("h",0.15);//Heat transfer coefficient
		
		
		
		
		//params.setDouble("se_tab_neg",1.5e8);
		//params.setDouble("se_tab_pos",0.5e8);
		//electrolyte conductivity and diffusion DEPEND ON C
		//params.setDouble("Ke",1.0e6);
		//params.setDouble("D_l",2.6e2);
		
		params.setDouble("t_0",0.2);
		//
		//inputs
		params.setDouble("totalTime",1e5);//70.0
		params.setDouble("dt",0.1);//25.0
		
		params.setDouble("iteration_ini",15);
		params.setDouble("iteration_n",10);
		
		//15; 45;75
		//params.setDouble("current",1e13/1.0e6);//pA
		params.setDouble("current",33.3e9);//pA C/s
	  //params.setDouble("current",0);//pA C/s
		params.setDouble("IpA",1);//pA/um^2
    params.setDouble("csmax_neg",28.7e-3);
		params.setDouble("csmax_pos",37.5e-3);
	  params.setDouble("cs_100_neg",0.915);
	  params.setDouble("cs_0_neg",0.02);
	  params.setDouble("cs_100_pos",0.022);
	  params.setDouble("cs_0_pos",0.98);
		params.setDouble("c2_ini",1.0e-3);

		params.setDouble("dis_top0",-0.24);
		
    
    battery<DIMS> problem(1,1,1, &params);
    problem.run ();
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    return 1;
  }

  return 0;
}
