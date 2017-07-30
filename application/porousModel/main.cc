#include "porousModel/battery_porousModel.h"


#define DIMS 3
#define NUM_THREADS 8
using namespace dealii;


int main (){
  try{
	
    deallog.depth_console (0);
		//set param values;;;unit:pmol; um;
		parametersClass params;
		//constants
		params.setInt("totalDOF", 8);
		
		params.setDouble("F", 96485.3329);
		params.setDouble("Rr", 8.3144598);
		params.setDouble("T_0", 298.0);

		//geometry
		params.setDouble("l_neg",60.0);
		params.setDouble("l_sep",23.0);
		params.setDouble("l_pos",45.0);
		
		params.setDouble("width",120e3);
		params.setDouble("height",85e3);
		//params.setDouble("w_b1",10);
		//params.setDouble("w_b2",10);
		// particle
		params.setDouble("alpha_neg",0.5);
		params.setDouble("alpha_pos",0.5);
		//reaction rate
		params.setDouble("k_neg",8.0e-4);
		params.setDouble("k_pos",8.0e-4);
		params.setDouble("t_0",0.2);
		//params.setDouble("k_pos",1);
    //R_s Radius of solid particals(m);esp_s,esp_l volume fraction active material/electrolyte
		//params.setDouble("eps_s_neg0",0.63);
		params.setDouble("eps_s_0_neg",0.53);
		params.setDouble("eps_s_0_pos",0.5);
		params.setDouble("eps_s_0_sep",0.35);
		
		params.setDouble("eps_l_0_neg",0.32);
		params.setDouble("eps_l_0_pos",0.35);
		params.setDouble("eps_l_0_sep",0.65);
		
		params.setDouble("eps_b_0_neg",0.15);
		params.setDouble("eps_b_0_pos",0.15);
		params.setDouble("eps_b_0_sep",0);
		
		params.setDouble("R_s_0_neg",8.0);//
		params.setDouble("R_s_0_pos",6.0);//
		params.setDouble("R_s_0_sep",0.0);//
		
		params.setDouble("pb",0.0);
		params.setDouble("pl",0.0);
		
		//solid diffusion
		params.setDouble("D_s_neg",0.5);
		params.setDouble("D_s_pos",0.1);
		//se, silid electronic conductivity
		params.setDouble("se_neg",1.5e8);
		params.setDouble("se_pos",0.5e8);
		
		params.setDouble("kappa_sep",0.42e-3);//0.4GPA
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
		
	  params.setDouble("omega_neg",9.615e-6);
		params.setDouble("omega_pos",6.025e-6);
		params.setDouble("omega_sep",82.46e-6);
		params.setDouble("omega_s",6e-6);
		
		
		params.setDouble("youngModule_neg",5.93e-3);//
		params.setDouble("youngModule_sep",0.5e-3);//
		params.setDouble("youngModule_pos",8.88e-3);//
		params.setDouble("nu", 0.3);

		params.setDouble("h",0.15);//Heat transfer coefficient
		//
		//inputs
		params.setDouble("totalTime",1e5);//70.0
		params.setDouble("dt",0.1);//25.0
		
		params.setDouble("iteration_ini",15);
		params.setDouble("iteration_n",10);
		params.setDouble("Fliptime", 295);
		
		//15; 45;75
		//params.setDouble("current",1e13/1.0e6);//pA
		params.setInt("currentflag", 1);
	  //params.setDouble("current",0);//pA C/s
		params.setDouble("IpA",1);//pA/um^2
    params.setDouble("c_li_max_neg",28.7e-3);
		params.setDouble("c_li_max_pos",37.5e-3);
	  params.setDouble("c_li_100_neg",0.915);
	  params.setDouble("c_li_0_neg",0.02);
	  params.setDouble("c_li_100_pos",0.022);
	  params.setDouble("c_li_0_pos",0.98);
		params.setDouble("c_li_plus_ini",1.0e-3);

		params.setDouble("dis_top0",-0.24);
		
    MultithreadInfo::set_thread_limit	(NUM_THREADS);
    battery_porousModel<DIMS> problem(3,3, params);
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