#include"../../../include/porousModel/battery_porousModel.h"

template <int dim>
battery_porousModel<dim>::battery_porousModel(const unsigned int _quad_electrode, const unsigned int _quad_electrolyte, parametersClass& _params)
	:battery<dim>(_quad_electrode, _quad_electrolyte, _params)
{
	std::cout<<"battery_porousModel base generated"<<std::endl;
}

template <int dim>
battery_porousModel<dim>::~battery_porousModel(){}

template class battery_porousModel<1>;
template class battery_porousModel<2>;
template class battery_porousModel<3>;



