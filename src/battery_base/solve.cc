#include"../../include/battery.h"

template <int dim>
void battery<dim>::solve(){
  double res=1, tol=1.0e-16, abs_tol=1.0e-13, initial_norm=0, current_norm=0;
  double machineEPS=1.0e-15;
  currentIteration=0;
  while (true){
    if (currentIteration>=5) {printf ("Maximum number of iterations reached without convergence. \n"); break; exit (1);}
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

template class battery<1>;
template class battery<2>;
template class battery<3>;