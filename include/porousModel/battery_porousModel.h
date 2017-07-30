#ifndef battery_porousModel_h
#define battery_porousModel_h

#include "../battery.h"
#include "porousModelResidual.h"
#include "porousModelFormula.h"


template <int dim>
class battery_porousModel: public battery<dim>
{
  public:
     battery_porousModel (const unsigned int _quad_electrode, const unsigned int _quad_electrolyte, parametersClass& _params);
    ~ battery_porousModel();
    void run ();
		
		porousModelFormula<Sacado::Fad::DFad<double>,dim>* electricChemoFormula;
		porousModelResidual<Sacado::Fad::DFad<double>,dim>* porousResidual;

		void setup_FeSystem();
		void make_grid();
    void setup_constraints();
    void apply_initial_condition();
    void assemble_system_interval (const typename hp::DoFHandler<dim>::active_cell_iterator &begin, const typename hp::DoFHandler<dim>::active_cell_iterator &end);
    void output_results (const unsigned int cycle) const;
		void step_load();
		unsigned int period;
		double periodTime;
		
		int currentflag;
};

#endif