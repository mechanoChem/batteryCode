#include"../../include/battery.h"

template <int dim>
void battery<dim>::make_grid()
{
	//geometry
	double l_neg=params->getDouble("l_neg");//120.0 negtive electrode
	double l_sep=params->getDouble("l_sep");//23.0
	double l_pos=params->getDouble("l_pos");//92.0 positive electrode
	double w_b1=params->getDouble("w_b1");//1e3.0
	double w_b2=params->getDouble("w_b2");//1e3.0
	
  //subdivided_hyper
  bool colorize = false;
  std::vector< std::vector< double > > step_sizes;
  step_sizes.resize(3);
  step_sizes[0].push_back(w_b1);
  step_sizes[2].push_back(w_b2);

  for(unsigned int cell=1;cell<=60;cell++) step_sizes[1].push_back(1);
  for(unsigned int cell=1;cell<=23;cell++) step_sizes[1].push_back(1);
  for(unsigned int cell=1;cell<=45;cell++) step_sizes[1].push_back(1);
	
	
  bX=w_b1; bY=l_neg+l_sep+l_pos; bZ=w_b2;
  electrode_Y1=l_neg;
	electrode_Y2=l_neg+l_sep;

  GridGenerator::subdivided_hyper_rectangle (triangulation, step_sizes, Point<dim>(0.0,0.0,0.0), Point<dim>(bX,bY,bZ), colorize);
}

template class battery<1>;
template class battery<2>;
template class battery<3>;