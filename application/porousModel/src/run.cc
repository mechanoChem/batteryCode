#include"porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::run()
{

	currentflag=battery<dim>::params->getInt("currentflag");
	periodTime=0;
  //running time
  battery<dim>::totalTime=battery<dim>::params->getDouble("totalTime");
	battery<dim>::totalDOF=battery<dim>::params->getInt("totalDOF");
	
  battery<dim>::currentIncrement=0;
	battery<dim>:: currentTime=0;

	double l_neg=battery<dim>::params->getDouble("l_neg");
	double l_sep=battery<dim>::params->getDouble("l_sep");
	double l_pos=battery<dim>::params->getDouble("l_pos");
	
  battery<dim>::bX=battery<dim>::params->getDouble("width"); 
	battery<dim>::bY=l_neg+l_sep+l_pos;
	battery<dim>::bZ=battery<dim>::params->getDouble("height");
  battery<dim>::electrode_Y1=l_neg;
	battery<dim>::electrode_Y2=l_neg+l_sep;
	
	porousModelFormula<Sacado::Fad::DFad<double>,dim> _electricChemoFormula;
	porousModelResidual<Sacado::Fad::DFad<double>,dim> _porousResidual;
	electricChemoFormula=&_electricChemoFormula;
	porousResidual=&_porousResidual;
	
	electricChemoFormula->reinit(*battery<dim>::params);
	porousResidual->reinit(*battery<dim>::params);
	
	setup_FeSystem();
  make_grid();
	battery<dim>::setMultDomain();
  battery<dim>::mark_boundary();
  battery<dim>::set_active_fe_indices ();
  setup_constraints();
  battery<dim>::setup_system();
  apply_initial_condition();
  battery<dim>::U=battery<dim>::Un;
	battery<dim>::U0=battery<dim>::Un;
  //apply_phi_e_BC();
  output_results(0); //output initial state
	double Fliptime=battery<dim>::params->getDouble("Fliptime");
  battery<dim>::dt=0.1;
  for (battery<dim>::currentTime=0; battery<dim>::currentTime<=battery<dim>::totalTime; battery<dim>::currentTime+=battery<dim>::dt){
    periodTime+=battery<dim>::dt;
    battery<dim>::currentIncrement++;
    if(battery<dim>::currentIncrement>Fliptime and period <1) {
      period++;                                                                                                                          
      battery<dim>::dt=0.1;
      periodTime=0;
    }
		
    if(periodTime>1) battery<dim>::dt=1;
    //if(currentTime>1) dt=10;
		//if(currentTime>=2.393e+04) dt=10;
		//if(currentTime>=2.394e+04) dt=40;
		//if(currentTime>=1000)  dt=100;
		//if(currentTime>3000)  dt=1;
		clock_t t_solve;	
  	t_solve = clock();
		battery<dim>::solve();
	  t_solve = clock() - t_solve;
		printf("************");
	  printf ("It took me %d clicks (%f seconds) for one solve.\n ",t_solve,((float)t_solve)/CLOCKS_PER_SEC);
	 	output_results(battery<dim>::currentIncrement);	
		if(battery<dim>::currentTime<=1) step_load();
	}
}


template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;