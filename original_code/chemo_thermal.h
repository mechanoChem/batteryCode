#ifndef CHEMO_H_
#define CHEMO_H_
#include <deal.II/dofs/dof_handler.h>
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "parameters.h"



//1st Transport equation (Sink term) (Mechanics residual implementation)
template <int dim>
void residualForChemo1(const FEValues<dim>& fe_values, unsigned int DOF, double dt, double currentTime, double totalTime,
	 										 dealii::Table<1, Sacado::Fad::DFad<double> >& R, 
											 deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, deformationMap<double, dim>& defMapConv,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& c1, dealii::Table<1,double>& c1_conv,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& jn,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& eps_s, dealii::Table<1,double>& eps_s_conv, 
											 parametersClass* params)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	double R_s_neg=params->getDouble("R_s_neg");
	double R_s_pos=params->getDouble("R_s_pos");
	double l_s=params->getDouble("l_s");
	double l_neg=params->getDouble("l_neg");
	double eps_s_neg0=params->getDouble("eps_s_neg0");
	double eps_s_pos0=params->getDouble("eps_s_pos0");
	double eps_s_0;
	const Point<dim> X = fe_values.quadrature_point(0);
  if (X[1]<l_neg) {
		eps_s_0=params->getDouble("eps_s_neg0");
  }
	//+
  if (X[1]>l_neg+l_s) {
		eps_s_0=params->getDouble("eps_s_pos0");
  }
	
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
	//R[i] +=  fe_values.shape_value(i, q)*jn[q]/eps_s[q]*fe_values.JxW(q);
				R[i] +=  fe_values.shape_value(i, q)*jn[q]*fe_values.JxW(q)*defMap.detF[q];
			}
    }
  }
	
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			//R[i] = R[i]*1.0e1; //test -Fajn
    }
  }
	
}

//2nd Transport equation (Diffusion with reaction (source term))
template <int dim>
void residualForChemo2(const FEValues<dim>& fe_values, unsigned int DOF, double dt, double currentTime, double totalTime,			 
											 dealii::Table<1, Sacado::Fad::DFad<double> >& R,
										   deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, deformationMap<double, dim>& defMapConv,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& c2, dealii::Table<1,double>& c2_conv, 
		        					 dealii::Table<2, Sacado::Fad::DFad<double> >& c2_j, 
		       						 dealii::Table<1, Sacado::Fad::DFad<double> >& T, dealii::Table<1, Sacado::Fad::DFad<double> >& jn,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& eps_l, dealii::Table<1,double>& eps_l_conv,
											 parametersClass* params)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

	double t_0=params->getDouble("t_0");
	double l_s=params->getDouble("l_s");
	double l_neg=params->getDouble("l_neg");
	dealii::Table<1, Sacado::Fad::DFad<double> > D_l(n_q_points);
	
	const Point<dim> X = fe_values.quadrature_point(0);

	for(unsigned int q=0; q<n_q_points; ++q){
	  // D_l[q]=0.3*6.5e-10*exp(-0.7*c2[q]*1.0e3)*1.0e12;
	   D_l[q]=std::pow(10,(-4.43-54/(T[q]-229-5e3*c2[q])-2.2e2*c2[q]))*1.0e8;
	  // D_l[q]=10e3;
	}
	
	dealii::Table<2,Sacado::Fad::DFad<double> > c2_j_Spat(n_q_points, dim);
	for (unsigned int q=0; q<n_q_points; ++q){
		for(unsigned int t=0;t<dim;t++) c2_j_Spat[q][t]=0.0;
		
		for(unsigned int j=0;j<dim;j++){
			for(unsigned int k=0; k<dim;k++){
				c2_j_Spat[q][j]+=c2_j[q][k]*defMap.invF[q][k][j];
			}
		}
	}
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	   	 	R[i] += fe_values.shape_value(i, q)*eps_l[q]*(c2[q]-c2_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];
			//R[i] += fe_values.shape_value(i, q)*c2[q]*eps_l[q]*(defMap.detF[q]-defMapConv.detF[q])/dt*fe_values.JxW(q)*defMap.detF[q];
			R[i] += fe_values.shape_value(i, q)*c2[q]*(eps_l[q]-eps_l_conv[q])/dt*fe_values.JxW(q)*defMap.detF[q];
		  	for (unsigned int j = 0; j < dim; j++){
					R[i] += fe_values.shape_grad(i, q)[j]*D_l[q]*std::pow(eps_l[q], 1.5)*c2_j_Spat[q][j]*fe_values.JxW(q)*defMap.detF[q];
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
	
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			R[i] = R[i]*1.0e1; //test -Fajn
    }
  }
	
}

//3nd Transport equation (phis; phi_s)
template <int dim>
void residualForChemo3(const FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values,
		       const typename hp::DoFHandler<dim>::active_cell_iterator &cell, double dt, double currentTime, double totalTime,double periodTime, int period, 
		       dealii::Table<1, Sacado::Fad::DFad<double> >& R, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap,
		       dealii::Table<1, Sacado::Fad::DFad<double> >& phis, dealii::Table<2,Sacado::Fad::DFad<double> >& phi_s_grad,
		       dealii::Table<1, Sacado::Fad::DFad<double> >& jn, dealii::Table<1, Sacado::Fad::DFad<double> >& eps_s,
		       parametersClass* params)
{


  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	
  //evaluate gradients
	dealii::Table<2,Sacado::Fad::DFad<double> > phis_j_Spat(n_q_points, dim);
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int j=0; j<dim; j++)  phis_j_Spat[q][j]=0.0;
		for(unsigned int j=0;j<dim;j++){
			for(unsigned int k=0; k<dim;k++){
				phis_j_Spat[q][j]+=phi_s_grad[q][k]*defMap.invF[q][k][j];
			}
		}
  }
	double F=params->getDouble("F");
	double l_neg=params->getDouble("l_neg");//120.0 negtive electrode
	double l_s=params->getDouble("l_s");//23.0
	double l_pos=params->getDouble("l_pos");//92.0 positive electrode
	double w_b1=params->getDouble("w_b1");
	double w_b2=params->getDouble("w_b2");
	double se;
	
  //se pore fraction (-)
  const Point<dim> X = fe_values.quadrature_point(0);
	//-
  if (X[1]<l_neg) {
		se=params->getDouble("se_neg");
  }
	//+
  if (X[1]>l_neg+l_s) {
		se=params->getDouble("se_pos");
  }
	
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	    for (unsigned int j = 0; j < dim; j++){
	      R[i] += -fe_values.shape_grad(i, q)[j]*se*eps_s[q]*(-phis_j_Spat[q][j])*fe_values.JxW(q)*defMap.detF[q];
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
	double current=params->getDouble("current");
	
	double IpA=params->getDouble("IpA");
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
	//if(ni==-1)IpA=0;
	//if(currentTime>=1) IpA=-380;
	dealii::Table<1, Sacado::Fad::DFad<double> > detF_surface(n_q_points);
	for(unsigned int q=0;q<n_q_points;q++){
		detF_surface[q]=defMap.F[q][0][0]*defMap.F[q][2][2]-defMap.F[q][0][2]*defMap.F[q][2][0];
	}

	
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell->face(faceID)->at_boundary()){
      const Point<dim> face_center = cell->face(faceID)->center();
	/*		
      if (std::abs(face_center[1])<1.0e-3){ //flux on Y boundary  
	    	fe_face_values.reinit (cell, faceID);
	    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
        	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	      	if (ck==0){
	        	for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
							//const Point<dim> X = fe_face_values.quadrature_point(q);
							//double fac_ipa=-X[0]*(X[0]-1000)/250000;
	          	R[i] += fe_face_values.shape_value(i, q)*(-IpA)*fe_face_values.JxW(q);              
	        	}
	      	}
	    	}
      }	
	*/		
      if (std::abs(face_center[1]-(l_neg+l_s+l_pos))<1.0e-3){ //flux on Y boundary  
	    	fe_face_values.reinit (cell, faceID);
	    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
        	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	      	if (ck==0){
	        	for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
	          	R[i] += fe_face_values.shape_value(i, q)*(IpA)*fe_face_values.JxW(q);              
	        	}
	      	}
	    	}
      }
		}	
  }
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			//R[i] = R[i]*1.0e-8; //test -Fajn
    }
  }
}

//4nd Transport equation (phie; phi_e)
template <int dim>
void residualForChemo4(const FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values,
											 double dt, double currentTime, double totalTime, 
											 dealii::Table<1, Sacado::Fad::DFad<double> >& R,  deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, 
											 dealii::Table<1, Sacado::Fad::DFad<double> >& c2, dealii::Table<2, Sacado::Fad::DFad<double> >& c2_j,
											 dealii::Table<1, Sacado::Fad::DFad<double> >& T, 
											 dealii::Table<1, Sacado::Fad::DFad<double> >& phie, dealii::Table<2,Sacado::Fad::DFad<double> >& phi_e_grad, 
											 dealii::Table<1, Sacado::Fad::DFad<double> >& jn, dealii::Table<1, Sacado::Fad::DFad<double> >& eps_l,
											 parametersClass* params)
{


  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	
  //evaluate gradients
	dealii::Table<2,Sacado::Fad::DFad<double> > phie_j_Spat(n_q_points, dim);
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int j=0; j<dim; j++) phie_j_Spat[q][j]=0.0;

		for(unsigned int j=0;j<dim;j++){
			for(unsigned int k=0; k<dim;k++){
				phie_j_Spat[q][j]+=phi_e_grad[q][k]*defMap.invF[q][k][j];
			}
		}
  }

	double Rr=params->getDouble("Rr");
	double F=params->getDouble("F");
	double t_0=params->getDouble("t_0");
	double l_s=params->getDouble("l_s");
	double l_neg=params->getDouble("l_neg");

	
	dealii::Table<1, Sacado::Fad::DFad<double> > Ke(n_q_points);
	for(unsigned int q=0; q<n_q_points; ++q){
		Ke[q]=((34.5*exp(-798/T[q])*std::pow((c2[q]*1.0e3),3)-485*exp(-1080/T[q])*std::pow((c2[q]*1.0e3),2)+2440*exp(-1440/T[q])*(c2[q]*1.0e3))/10)*1.0e6;
//Ke[q]=((34.5*exp(-798/T)*std::pow((c2[q]),3)-485*exp(-1080/T)*std::pow((c2[q]),2)+2440*exp(-1440/T)*(c2[q]))/10)*1.0e6;
	}
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				for (unsigned int j = 0; j < dim; j++){
					R[i] += -fe_values.shape_grad(i, q)[j]*(Ke[q]*std::pow(eps_l[q], 1.5))*(-phie_j_Spat[q][j]+2.0*Rr*T[q]/F*(1.0-t_0)/c2[q]*(c2_j[q][j]))*fe_values.JxW(q)*defMap.detF[q];
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
	
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			//R[i] = R[i]*1.0e-5; //test -Fajn
    }
  }
}

template <int dim>
void residualForChemoThermal(const FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, Table<1, Sacado::Fad::DFad<double> >& ULocal, 
											 			 const typename hp::DoFHandler<dim>::active_cell_iterator &cell, double dt, double currentTime, double totalTime,
											       dealii::Table<1, Sacado::Fad::DFad<double> >& R,  deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, 
											 			 dealii::Table<1,double>& T_conv, dealii::Table<1, Sacado::Fad::DFad<double> >& T,
											       dealii::Table<2,Sacado::Fad::DFad<double> >& T_grad,
											       dealii::Table<1, Sacado::Fad::DFad<double> >& Q,
											       parametersClass* params)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	
	double l_neg=params->getDouble("l_neg");//120.0 negtive electrode
	double l_s=params->getDouble("l_s");//23.0
	double l_pos=params->getDouble("l_pos");//92.0 positive electrode
	double w_b1=params->getDouble("w_b1");//x
	double w_b2=params->getDouble("w_b2");//z
	
	double density_neg=params->getDouble("density_neg");
	double density_sep=params->getDouble("density_sep");
	double density_pos=params->getDouble("density_pos");
	double Cp=params->getDouble("Cp");
	double lambda_neg=params->getDouble("lambda_neg");
	double lambda_sep=params->getDouble("lambda_sep");
	double lambda_pos=params->getDouble("lambda_pos");
	double h=params->getDouble("h");
	double T_0=params->getDouble("T_0");
	double l;
	
	double density, lambda;
	//double length=l_neg+l_s+l_pos;
	
  const Point<dim> X = fe_values.quadrature_point(0);
	//-
  if (X[1]<l_neg) {
		density=params->getDouble("density_neg");
		//Cp=params->getDouble("Cp_neg");
		lambda=params->getDouble("lambda_neg");
		l=params->getDouble("l_neg");
  }
	else if(X[1]<l_neg+l_s){
		density=params->getDouble("density_sep");
		//Cp=params->getDouble("Cp_sep");
		lambda=params->getDouble("lambda_sep");
		l=params->getDouble("l_s");
	}
	else if (X[1]>l_neg+l_s) {
		density=params->getDouble("density_pos");
		//Cp=params->getDouble("Cp_pos");
		lambda=params->getDouble("lambda_pos");
		l=params->getDouble("l_pos");
  }
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
					R[i] +=  fe_values.shape_value(i, q)*(1/dt)*density*Cp*(T[q]-T_conv[q])*fe_values.JxW(q)*defMap.detF[q];	

					R[i] += fe_values.shape_grad(i, q)[0]*lambda*T_grad[q][0]*fe_values.JxW(q)*defMap.detF[q];
					R[i] += fe_values.shape_grad(i, q)[1]*lambda*T_grad[q][1]*fe_values.JxW(q)*defMap.detF[q];
					R[i] += fe_values.shape_grad(i, q)[2]*lambda*T_grad[q][2]*fe_values.JxW(q)*defMap.detF[q];
				
					R[i] -=  fe_values.shape_value(i, q)*Q[q]*fe_values.JxW(q)*defMap.detF[q];
      }
    }
  }
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell->face(faceID)->at_boundary()){ 
			double length=1;

	    fe_face_values.reinit (cell, faceID);
			const unsigned int n_face_quadrature_points = fe_face_values.n_quadrature_points;
			dealii::Table<1,Sacado::Fad::DFad<double> > T_face(n_face_quadrature_points);
			for(unsigned int q=0;q<n_face_quadrature_points;q++){
				T_face[q]=0;
			  for (unsigned int i=0; i<dofs_per_cell; ++i) {
			    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
					if(ck==0){
						T_face[q]+=fe_face_values.shape_value(i, q)*ULocal[i];
					}
				}
			}
			
	   	for (unsigned int i=0; i<dofs_per_cell; ++i) {
       	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	     	if (ck==0){
	       	for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
	         	R[i] += fe_face_values.shape_value(i, q)*h*(T_face[q]-T_0)*fe_face_values.JxW(q);              
	       	}
	     	}
     	}
		}	
  }
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			//R[i] = R[i]*1.0e-7; //test -Fajn
    }
  }
	
	
}

Sacado::Fad::DFad<double> GetOpenCirculatePotential(Sacado::Fad::DFad<double> x,int flag)
{
	Sacado::Fad::DFad<double> value=0.0;
	//negtive electrode
	if(flag==-1){
		double c1 =0.265697795660872, c2 =0.555104680448495, c3 =178.97991682123549, c4 =0.012357396765331, c5 =0.55727132360432;
		double c6 =0.028219099799268, c7 =0.011683704080029, c8 =0.239318990250894, c9 =0.048646992277392, c10 =0.012910660088849;
		double c11 =0.174909417938192, c12 =0.03483002163646, c13 =0.050098062010346, c14 =0.024505122678677, c15 =0.03529369961247;
		double c16 =0.011931381413342, c17 =0.12992241633878, c18 =0.019797869695897, c19 =0.152636640731331, c20 =0.030000057933125, c21 =0.022725508023415;
		
		value=c1+c2*exp(-c3*x)-c4*std::tanh((x-c5)/c6)-c7*std::tanh((x-c8)/c9)-c10*std::tanh((x-c11)/c12)-c13*std::tanh((x-0.99)/c14)-c15*x-c16*std::tanh((x-c17)/c18)-c19*std::tanh((x-c20)/c21);//
		//value=-0.132+1.41*std::exp(-3.52*x);
	}
	//positive electrode
	else if(flag==1){
	 double	b1= -0.0922859116552415, b2= -7.8680409125385697, b3= 50.072175799512607, b4= -122.28161948058685, b5= 82.985110649682696;
	 double b6= 140.2938943391359, b7= -374.73497214300698, b8= 403.2463575744942, b9= -221.19151490076541, b10= 49.3392659530526530;
	 double	b11= -0.0217591621507594, b12= -1.9006524442210881, b13= 11.726362513914014, b14= -28.784794013633256, b15= 27.542704665893613;
	 double b16= -8.6342730487746202;
		value = (b1+b2*x+b3*std::pow(x,2)+b4*std::pow(x,3)+b5*std::pow(x,4)+b6*std::pow(x,5)+b7*std::pow(x,6)+b8*std::pow(x,7)+b9*std::pow(x,8)+b10*std::pow(x,9))/(b11+b12*x+b13*std::pow(x,2)+b14*std::pow(x,3)+b15*std::pow(x,4)+b16*std::pow(x,5));
	}
	return value;
}

Sacado::Fad::DFad<double> GetQrev(Sacado::Fad::DFad<double> x)
{
	Sacado::Fad::DFad<double> value=0.0;
	//negtive electrode
	if(x<=0.2)  value=0.01442*x*x-0.00291*x-0.000138;
	else if(x<=0.4) value=0.00634*x*x*x-0.006625*x*x+0.002635*x-0.0004554;
	else if(x<=0.5) value=0.001059*x-0.0004793;
	else if(x<=0.7) value=0.00025*x-7.5e-5;
	else if(x<=0.8) value=-0.001*x+0.0008;
	else if(x<=0.85) value=0.0333*x*x-0.057*x+0.02427;
	else if(x<=0.95) value=0.002*x*x-0.0039*x+0.00177;
  	else if (x<=1) value=-0.0014*x+0.0012;

	return value;
}

Sacado::Fad::DFad<double> GetBeta_s(Sacado::Fad::DFad<double> UnitC,int flag)
{
	Sacado::Fad::DFad<double> value=0.0;
	value=1.496*std::pow(UnitC,3)-1.739*UnitC*UnitC+1.02*UnitC-0.03304*std::exp(2.972*UnitC)-0.04587*tanh((UnitC-0.1)/0.1)-0.003608*tanh((UnitC-0.3)/0.1)+0.0214*tanh((UnitC-0.65)/0.1);
	if(flag==1) value=value-0.0939;
	return value;
}

Sacado::Fad::DFad<double> GetBeta(Sacado::Fad::DFad<double> UnitC,int flag)
{
	Sacado::Fad::DFad<double> value=0.0;
	value=0.1315*std::pow(UnitC,4)-0.1989*std::pow(UnitC,3)+0.06481*UnitC*UnitC+0.02793*UnitC-0.000655;
	if(flag==1) value=value-0.019;
	return value;
}



#endif /* CHEMO_H_ */
