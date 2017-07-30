#include"porousModel/porousModelResidual.h"

template <class T, int dim>
void porousModelResidual<T,dim>::IncludeStrain_chemo(dealii::Table<1, T >& _fac_cli)
{
	fac_cli=_fac_cli;
	 includeSt_chemo=true;
}


template <class T, int dim>
void porousModelResidual<T,dim>::IncludeStrain_thermal(dealii::Table<1, T >& _fac_T)
{
	fac_T=_fac_T;
	includeSt_thermal=true;
}


template <class T, int dim>
void porousModelResidual<T,dim>::evaluateStrain(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	
	TableIndices<3> indices(n_q_points, dim, dim);
	this->Fe.reinit(indices);
	this->E.reinit(indices);

	for (unsigned int q=0; q<n_q_points; ++q){
  	for (unsigned int i=0; i<dim; ++i){
    	for (unsigned int j=0; j<dim; ++j){
				this->Fe[q][i][j]=0.0;
				this->E[q][i][j]=0.0;
				
	  		this->Fe[q][i][j]=defMap.F[q][i][j];
	  		if(includeSt_thermal==true) this->Fe[q][i][j]=this->Fe[q][i][j]/fac_T[q]; //iso-
    	}
  	}  
  	if(includeSt_chemo==true) this->Fe[q][1][1]=this->Fe[q][1][1]/fac_cli[q]; //uni-axis
  	//for(unsigned int i=0;i<dim;i++) Fe[i][i]=Fe[i][i]/fac_T[q];
  
  	//E
  	for (unsigned int i=0; i<dim; ++i){
    	for (unsigned int j=0; j<dim; ++j){
				this->E[q][i][j] = -0.5*(i==j);
				for (unsigned int k=0; k<dim; ++k){
	  			this->E[q][i][j] += 0.5*this->Fe[q][k][i]*this->Fe[q][k][j];
				}
    	}
    	//E[i][i]+=-Ec;
  	}
	}
}

template class porousModelResidual<Sacado::Fad::DFad<double>, 1>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 2>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 3>;				