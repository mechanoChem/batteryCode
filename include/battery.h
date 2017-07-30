#ifndef battery_h
#define battery_h

#include <memory>
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

#include "model.h"
#include "supplementary/dataStruct.h"
#include "supplementary/supplementaryFunctions.h"
#include "supplementary/parameters.h"
#include "supplementary/functionEvaluations.h"


using namespace dealii;

/**
*base class define generic functions of FEM for battery problem based on deal.ii
*/
template <int dim>
class battery
{
  public:
		/**
		*constructor
		*/
    battery (const unsigned int quad_electrode, const unsigned int quad_electrolyte, parametersClass& _params);
		/**
		*destructor
		*/
    ~battery();
		
		/**
		* initial function for users
		*/
    virtual void run ();
		
		/**
		*params class to store all input parameters of material and solver
		*/
		parametersClass* params;
		
		
		/**
		*equations solved for different domains for sandwitch model
		*domain id used for Fe_Nothing
		*/
    enum
    {
      electrode_domain_id,
      electrolyte_domain_id
    };
		/**
		*determine the domiain id of each cell 
		*/
    static bool
    cell_is_in_electrode_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);
    static bool
    cell_is_in_electrolyte_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell); 
    
		/**
		*set up FeSystem of deal.ii 
		*/
		virtual void setup_FeSystem();
		/**
		*set up different domain corresponding to domain id  
		*/
		void setMultDomain();
		/**
		*generate mesh
		*/
		virtual void make_grid();
		/**
		*mark_boundary:mark cell surface if they are domain boundary 
		*/
    void mark_boundary();
		
		/**
		*set_active_fe_indices:set which FeSystem the cell use (based on the what's domain the cell in) 
		*/
    void set_active_fe_indices ();
		
		/**
		*using constraints to set up dirichlet boundary condition including interface boundary condition
		*Note, for newton method, we really solve dU rather than U. So the dirichlet boundary condition is "zero"
		*the really value goes to initial condition.
		*/
    virtual void setup_constraints();
		/**
		*setup solver system: initialize sparsePattern, Jacobbian matrix and vectors
		*/
    void setup_system();
		/**
		*applying initial condition
		*Note for time independent problem, initial condition is still needed when iteration solver is used
		*/
    virtual void apply_initial_condition();
		/**
		*initialize before assemble and set multiple threads
		*/
    void assemble_system ();
    /**
		*assemble system interval
		*update all variables that needed for equations (formula class)
		*assemble residual (redisual class)
		*/
		virtual void assemble_system_interval (const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end);
    /**
		*slover
		*/
		void solve();
		/**
		*generate and output results
		*/
    virtual void output_results (const unsigned int cycle) const;


	  /**
		*essential definition of deal.ii
		*Triangulation: for mesh
		*/ 
    Triangulation<dim>    triangulation;
		
	  /**
		*FESystem of electrode domain: for evaluate basis functions
		*/ 
		std::shared_ptr<FESystem<dim>> electrode_fe;
		
	  /**
		*FESystem of separator domain: for evaluate basis functions
		*/ 
		std::shared_ptr<FESystem<dim>> electrolyte_fe;
		
	  /**
		*QGauss: rule of gauss quadrature point of electrode domain
		*/ 
    const QGauss<dim> electrode_quadrature;
		
	  /**
		*QGauss: rule of gauss quadrature point of electrode domain
		*/ 
    const QGauss<dim> electrolyte_quadrature;
	  /**
		*QGauss: rule of gauss quadrature point of surface(2D)
		*/ 
		const QGauss<dim-1> common_face_quadrature;
	  /**
		*FECollection: collection of FESystem used for FE_Nothing 
		*/ 
    hp::FECollection<dim> fe_collection;
		
	  /**
		*essential definition of deal.ii used for cell interator 
		*/ 
    hp::DoFHandler<dim>   dof_handler;
	  /**
		*QCollection: QGauss collection
		*/ 
    hp::QCollection<dim>  q_collection;
	  /**
		*ConstraintMatrix: constraints 
		*/ 
    ConstraintMatrix      constraints;
		
	  /**
		*essential definition of deal.ii used for cell interator 
		*/ 
		SparsityPattern       sparsity_pattern;
		
		/**
		*include base class for Residual (governing equations)
		*/
		model<Sacado::Fad::DFad<double>, dim>* mechanicalResidual;
    
		/**
		*system_matrix:jacobian materix
		*/
    SparseMatrix<double> system_matrix;
		/**
		*system_rhs:right hand side of linear system 
		*/
  	Vector<double> system_rhs;
		/**
		*U: current solution
		*/
		Vector<double> U;
		/**
	 	*Un: solution of previous step iteration; 
		*/
		Vector<double> Un;
		/**
		*dU:step change of during iteration used in Newton method
		*/
		Vector<double> dU;
		
		/**
		*U0:store initial U
		*/
		Vector<double> U0;
		
		/**
		*total degree of freedom
		*/
		int totalDOF;
       
    /**
		*length of X direction
		*/
    double bX;
		
    /**
		*length of Y direction
		*/
		double bY;
    /**
		*length of Z direction
		*/
		double bZ;
		
    /**
		*length of first electrode/separator interface from bottom 
		*/
    double electrode_Y1;
		
    /**
		*length of second electrode/separator interface from bottom 
		*/
		double electrode_Y2;
		
		/**
		*currentIncrement: current increment of loading
		*/
		unsigned int currentIncrement;
		/**
		*currentIteration:current time step in iteration algorithm 
		*/
		unsigned int currentIteration;
		
		/**
		*total run time
		*/
    double totalTime;
		
		/**
		*current time of simulation
		*/
		double currentTime;
		
		/**
		*time step
		*/
		double dt;
		/**
		*deal.iii threads manager
		*/
    Threads::Mutex assembler_lock;
};

#endif