#include"../../include/battery.h"

template <int dim>
void battery<dim>::run()
{
	model<Sacado::Fad::DFad<double>,dim> _mechanicalResidual;  //running time
  totalTime=300000.0; dt=100.0;
  currentIncrement=0; currentTime=0;
	totalDOF=8;
	
	
	mechanicalResidual=&_mechanicalResidual;
	setup_FeSystem();
  make_grid();
	setMultDomain();
  mark_boundary();
  set_active_fe_indices ();
  setup_constraints();
  setup_system();
  apply_initial_condition();
  U=Un; U0=Un;
  //apply_phi_e_BC();
  output_results(0); //output initial state
  
 	
  currentIncrement=0;
  for (currentTime=0; currentTime<=totalTime; currentTime+=dt){
    currentIncrement++;
   	 solve();
		 
    //if(currentIncrement%5==0) output_results(currentIncrement); //output solution at current load increment
    output_results(currentIncrement);	
  }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;