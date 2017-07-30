#include"porousModel/porousModelResidual.h"

template <class T, int dim>
void porousModelResidual<T,dim>::residualForChemoCli(const FEValues<dim>& fe_values,  unsigned int DOF, dealii::Table<1, T >& R, 
																										 deformationMap<T, dim>& defMap,
	 																									 deformationMap<double, dim>& defMapConv, dealii::Table<1, T >& c1,
																										 dealii::Table<1,double>& c1_conv, dealii::Table<1, T >& jn,
											 															 dealii::Table<1, T >& eps_s, dealii::Table<1,double>& eps_s_conv)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	double dt=this->dt;
	
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	 		 	R[i] += fe_values.shape_value(i, q)*eps_s[q]*(c1[q]-c1_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];
			 // R[i] += fe_values.shape_value(i, q)*c1[q]*eps_s[q]*(defMap.detF[q]-defMapConv.detF[q])/dt*fe_values.JxW(q)*defMap.detF[q];
			  R[i] += fe_values.shape_value(i, q)*c1[q]*(eps_s[q]-eps_s_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];  
      }
    }
  }
	
 	//source term
 	for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
      for (unsigned int q=0; q<n_q_points; ++q){
	  		R[i] +=  fe_values.shape_value(i, q)*jn[q]*fe_values.JxW(q)*defMap.detF[q];
			}
    }
  }
}
	
template <class T, int dim>
void porousModelResidual<T,dim>::residualForChemoCli_plus(const FEValues<dim>& fe_values,  unsigned int DOF, dealii::Table<1, T >& R,
										   																deformationMap<T, dim>& defMap, deformationMap<double, dim>& defMapConv,
											 															  dealii::Table<1, T >& c2, dealii::Table<1,double>& c2_conv, 
		        					 															  dealii::Table<2, T >& c2_j, dealii::Table<1, T >& jn,
											 															  dealii::Table<1, T >& eps_l, dealii::Table<1,double>& eps_l_conv,dealii::Table<1,T > D_l)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	double t_0=this->params->getDouble("t_0");
	double dt=this->dt;

  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	   	 	R[i] += fe_values.shape_value(i, q)*eps_l[q]*(c2[q]-c2_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];
				//R[i] += fe_values.shape_value(i, q)*c2[q]*eps_l[q]*(defMap.detF[q]-defMapConv.detF[q])/dt*fe_values.JxW(q)*defMap.detF[q];
				R[i] += fe_values.shape_value(i, q)*c2[q]*(eps_l[q]-eps_l_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];
		  	for (unsigned int j = 0; j < dim; j++){
					R[i] += fe_values.shape_grad(i, q)[j]*D_l[q]*std::pow(eps_l[q], 1.5)*c2_j[q][j]*fe_values.JxW(q)*defMap.detF[q];
		  	}
      }
    }
  }

  //source term
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck==0){
      for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=  -fe_values.shape_value(i, q)*(1.0-t_0)*jn[q]*fe_values.JxW(q)*defMap.detF[q];
      }
    }
  }																												
}


template <class T, int dim>
void porousModelResidual<T,dim>::residualForChemoPhis(const FEValues<dim>& fe_values, unsigned int DOF, 
		       																					 	dealii::Table<1,T >& R, deformationMap<T, dim>& defMap,
		       																					 	dealii::Table<1,T >& phis, dealii::Table<2,T >& phi_s_grad,
		       																					  dealii::Table<1,T >& jn, dealii::Table<1,T >& eps_s)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	
	double F=this->params->getDouble("F");
	double l_neg=this->l_neg;
	double l_sep=this->l_sep;
	double l_pos=this->l_pos;
	double se;
	
  //se pore fraction (-)
  const Point<dim> X = fe_values.quadrature_point(0);
	//-
  if (X[1]<l_neg) {
		se=this->params->getDouble("se_neg");
  }
	//+
  if (X[1]>l_neg+l_sep) {
		se=this->params->getDouble("se_pos");
  }
	
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	    	for (unsigned int j = 0; j < dim; j++){
					R[i] += -fe_values.shape_grad(i, q)[j]*se*eps_s[q]*(-phi_s_grad[q][j])*fe_values.JxW(q)*defMap.detF[q];
	    	}
      }
    }
  }
 //source term
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
      for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=  fe_values.shape_value(i, q)*F*jn[q]*fe_values.JxW(q)*defMap.detF[q]; //test -Fajn
      }
    }
  }
}
	
template <class T, int dim>
void porousModelResidual<T,dim>::residualForPhisBC(const FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values,
																									 unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap,unsigned int period, double periodTime)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	double IpA=this->params->getDouble("IpA");
	double ni;
  if(period%2==0)ni=1;
  else ni=-1;	
	if(periodTime>=0.1) IpA=5*ni;
	if(periodTime>=0.2) IpA=10*ni;
	if(periodTime>=0.3) IpA=22*ni;
	if(periodTime>=0.4) IpA=70*ni;
	if(periodTime>=0.5) IpA=100*ni;
	if(periodTime>=0.6) IpA=110*ni;
	if(periodTime>=0.7) IpA=160*ni;
	if(periodTime>=0.8) IpA=190*ni;
	if(periodTime>=0.9) IpA=220*ni;
	
	dealii::Table<1,T > detF_surface(n_q_points);
	for(unsigned int q=0;q<n_q_points;q++){
		detF_surface[q]=defMap.F[q][0][0]*defMap.F[q][2][2]-defMap.F[q][0][2]*defMap.F[q][2][0];
	}
	
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
  	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
  	if (ck==0){
    	for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
      	R[i] += fe_face_values.shape_value(i, q)*(IpA)*fe_face_values.JxW(q);              
    	}
  	}
	}
}


template <class T, int dim>
void porousModelResidual<T,dim>::residualForChemoPhie(const FEValues<dim>& fe_values, unsigned int DOF,
											 																dealii::Table<1,T >& R,  deformationMap<T, dim>& defMap, 
											 															 	dealii::Table<1,T >& c2, dealii::Table<2,T >& c2_j,
											 															  dealii::Table<1,T >& Temp, 
											 															  dealii::Table<1,T >& phie, dealii::Table<2,T >& phi_e_grad, 
											 															  dealii::Table<1,T >& jn, dealii::Table<1,T >& eps_l, dealii::Table<1,T >Ke)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	double dt=this->dt;

	double Rr=this->params->getDouble("Rr");
	double F=this->params->getDouble("F");
	double t_0=this->params->getDouble("t_0");

  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				for (unsigned int j = 0; j < dim; j++){
					R[i] += -fe_values.shape_grad(i, q)[j]*(Ke[q]*std::pow(eps_l[q], 1.5))*(-phi_e_grad[q][j]+2.0*Rr*Temp[q]/F*(1.0-t_0)/c2[q]*(c2_j[q][j]))*fe_values.JxW(q)*defMap.detF[q];
				}
      }
    }
  }

  //source term
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
      for (unsigned int q=0; q<n_q_points; ++q){
	    R[i] += -fe_values.shape_value(i, q)*F*jn[q]*fe_values.JxW(q)*defMap.detF[q]; //test Fajn
      }
    }
  }
	
}

template <class T, int dim>
void porousModelResidual<T,dim>::residualForChemoThermal(const FEValues<dim>& fe_values, unsigned int DOF,
											       														 dealii::Table<1, T >& R,  deformationMap<T, dim>& defMap, 
											 			 														 dealii::Table<1, T >& Temp, dealii::Table<1,double>& T_conv, 
											       														 dealii::Table<2,T >& T_grad,
											      														 dealii::Table<1, T >& Q)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	double dt=this->dt;
	double l_neg=this->l_neg;//120.0 negtive electrode
	double l_sep=this->l_sep;//23.0
	
	double Cp=this->params->getDouble("Cp");

	double density, lambda;
	//double length=l_neg+l_s+l_pos;
	
  const Point<dim> X = fe_values.quadrature_point(0);
	//-
  if (X[1]<l_neg) {
		density=this->params->getDouble("density_neg");
		lambda=this->params->getDouble("lambda_neg");
  }
	else if(X[1]<l_neg+l_sep){
		density=this->params->getDouble("density_sep");
		lambda=this->params->getDouble("lambda_sep");
	}
	else if (X[1]>l_neg+l_sep) {
		density=this->params->getDouble("density_pos");
		lambda=this->params->getDouble("lambda_pos");
  }
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
					R[i] +=  fe_values.shape_value(i, q)*(1/dt)*density*Cp*(Temp[q]-T_conv[q])*fe_values.JxW(q)*defMap.detF[q];	

					R[i] += fe_values.shape_grad(i, q)[0]*lambda*T_grad[q][0]*fe_values.JxW(q)*defMap.detF[q];
					R[i] += fe_values.shape_grad(i, q)[1]*lambda*T_grad[q][1]*fe_values.JxW(q)*defMap.detF[q];
					R[i] += fe_values.shape_grad(i, q)[2]*lambda*T_grad[q][2]*fe_values.JxW(q)*defMap.detF[q];
				
					R[i] -=  fe_values.shape_value(i, q)*Q[q]*fe_values.JxW(q)*defMap.detF[q];
      }
    }
  }
}


template <class T, int dim>
void porousModelResidual<T,dim>::residualForHeatConv(const FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values,
																												 unsigned int DOF, dealii::Table<1, T >& R,dealii::Table<1, T >& T_face)
{
	double h=this->params->getDouble("h");
	double T_0=this->params->getDouble("T_0");
	
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
 	for (unsigned int i=0; i<dofs_per_cell; ++i) {
   	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
   	if (ck==0){
     	for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
       	R[i] += fe_face_values.shape_value(i, q)*h*(T_face[q]-T_0)*fe_face_values.JxW(q);              
     	}
   	}
 	}																								 	
}
																													 
template class porousModelResidual<Sacado::Fad::DFad<double>, 1>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 2>;
template class porousModelResidual<Sacado::Fad::DFad<double>, 3>;																									 
																													 
																													 
																													 
																													 