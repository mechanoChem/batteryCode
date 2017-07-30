#include"../../include/battery.h"

template <int dim>
battery<dim>::battery (const unsigned int quad_electrode, const unsigned int quad_electrolyte, parametersClass& _params)
	: params(&_params), dof_handler (triangulation),electrode_quadrature(quad_electrode), electrolyte_quadrature(quad_electrolyte),common_face_quadrature(quad_electrolyte)
{

}

template <int dim>
battery<dim>::~battery (){dof_handler.clear ();}


template class battery<1>;
template class battery<2>;
template class battery<3>;