#ifndef battery_h
#define battery_h




template <int dim>
class battery_particle
{
  public:
     battery_particle (const unsigned int li_degree, const unsigned int phi_degree,const unsigned int elasticity_degree,parametersClass* _params);
    ~ battery_particle();
    void run ();
		parametersClass* params;

    enum
    {
      electrode_domain_id,
      electrolyte_domain_id
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
    void assemble_system ();
    void assemble_system_interval (const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end);
		void assemble_interface_electrode_term(const typename DoFHandler<dim>::active_cell_iterator &cell, unsigned int f, dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal, Table<1, double > &ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R);
		void assemble_interface_electrolyte_term(const typename DoFHandler<dim>::active_cell_iterator &cell, unsigned int f, dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal, Table<1, double > &ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R);
		
		void solve();
    void output_results (const unsigned int cycle) const;
    void apply_Dirichlet_BC();
		void get_potential();
	
    const unsigned int li_degree;
    const unsigned int phi_degree;
    const unsigned int elasticity_degree;
      
    Triangulation<dim>    triangulation;
    FESystem<dim>         electrode_fe;
    FESystem<dim>         electrolyte_fe;
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
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    Threads::ThreadMutex assembler_lock;
};

#endif