#include "porousModel/battery_porousModel.h"

template <int dim>
void battery_porousModel<dim>::make_grid()
{
	
  //subdivided_hyper
  bool colorize = false;
  std::vector< std::vector< double > > step_sizes;
  step_sizes.resize(3);
  step_sizes[0].push_back(battery<dim>::bX);
  step_sizes[2].push_back(battery<dim>::bZ);
  //for(unsigned int cell=1;cell<=10;cell++) step_sizes[0].push_back(2);
  //for(unsigned int cell=1;cell<=10;cell++) step_sizes[2].push_back(2);
	
  for(unsigned int cell=1;cell<=60;cell++) step_sizes[1].push_back(1);
  for(unsigned int cell=1;cell<=23;cell++) step_sizes[1].push_back(1);
  for(unsigned int cell=1;cell<=45;cell++) step_sizes[1].push_back(1);
	
	/*
	step_sizes[1].push_back(tab_neg);
	step_sizes[1].push_back(l_neg);
	step_sizes[1].push_back(l_s);
	step_sizes[1].push_back(l_pos);
	step_sizes[1].push_back(tab_pos);
	*/

  GridGenerator::subdivided_hyper_rectangle (battery<dim>::triangulation, step_sizes, Point<dim>(0.0,0.0,0.0), Point<dim>(battery<dim>::bX,battery<dim>::bY,battery<dim>::bZ), colorize);
}


template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;