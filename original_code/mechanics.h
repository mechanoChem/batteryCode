/*
 * mechanics.h
 *  Created on: Apr 21, 2011
 */

#ifndef MECHANICS_H_
#define MECHANICS_H_
#include <deal.II/dofs/dof_handler.h>
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "quadraturePointData.h"
//#include "quadraturePointData_hp.h"

//Saint-Venant Kirchhoff constitutive model
template <int dim>
inline double SVK3D(unsigned int i, unsigned int j, unsigned int k, unsigned int l, double E){
  double nu=0.3;
  double lambda=(E*nu)/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
}

template <int dim>
inline double SVK2D(unsigned int i, unsigned int j, double E){
  double nu=0.3;
  double lambda=(E*nu)/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
  if (i==j && i<2) return lambda + 2*mu;
  else if (i==2 && j==2) return mu;
  else if ((i+j)==1) return lambda;
  else return 0.0;
}

template <int dim>
inline double SVK1D(double E){
  return E;
}

//Mechanics implementation
template <class T, int dim>
  void getDeformationMap(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap, const unsigned int iteration){
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //evaluate dx/dX
  Table<3, T> gradU(n_q_points, dim, dim);
  evaluateVectorFunctionGradient<T, dim>(fe_values, DOF, ULocal, gradU);

  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.F[q][i][j] = Fq[i][j] = (i==j) + gradU[q][i][j]; //F (as double value)
      }
    }
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.invF[q][i][j] = invFq[i][j];
      }
    }
    //detF
    /*
    if (defMap.detF[q].val()<=1.0e-15 && iteration==0){
     printf("**************Non positive jacobian detected**************. Value: %12.4e\n", defMap.detF[q].val());
    for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j) printf("%12.6e  ", defMap.F[q][i][j].val());
				printf("\n"); exit(-1);
      }
      //throw "Non positive jacobian detected";
    }
    */
  }
}

//Mechanics implementation
template <class T, int dim>
  void evaluateStress(const FEValues<dim>& fe_values,const unsigned int DOF, const Table<1, T>& ULocal, Table<3, T>& P, const deformationMap<T, dim>& defMap, double youngsModulus, typename hp::DoFHandler<dim>::active_cell_iterator& cell, dealii::Table<1,double>& c1_conv, dealii::Table<1, Sacado::Fad::DFad<double> >& c1, dealii::Table<1, Sacado::Fad::DFad<double> >& c1_0, dealii::Table<1, Sacado::Fad::DFad<double> >& Temp, double dt, double currentTime){
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  dealii::Table<1,Sacado::Fad::DFad<double> > fac(n_q_points),fac_T(n_q_points), x_fac(n_q_points), beta(n_q_points);

  double c1_max=28.7e-3;
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){

		
    //determine quadrature points mapping position
    fac[q]=1.0; // sep
    fac_T[q]=1.0;	

    //double beta = 1.e-3;

    const Point<dim> X = fe_values.quadrature_point(q);

    // -
    if (X[1]<60.) {
	Sacado::Fad::DFad<double> UnitC;
	UnitC=c1[q]/c1_max;
	fac[q]=0.1315*std::pow(UnitC,4)-0.1989*std::pow(UnitC,3)+0.06481*UnitC*UnitC+0.02793*UnitC-0.000655+1;
	fac[q]=fac[q]-0.019;
	fac[q]=1.01;
	//ac_T[q]=std::pow((Temp[q].val()-298)*9.615e-6+1,1.0/3.0);
    }
     //else if (X[1]<83) fac_T[q]=std::pow((Temp[q].val()-298)*82.46e-6+1,1.0/3.0);
    //else if (X[1]>83) fac_T[q]=std::pow((Temp[q].val()-298)*6.025e-6+1,1.0/3.0);
	if(!(fac_T[q]>0)) fac_T[q]=1;
    //Fe
    Table<2, Sacado::Fad::DFad<double> > Fe (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
		  Fe[i][j]=defMap.F[q][i][j];
		  Fe[i][j]=defMap.F[q][i][j]/fac_T[q]; //iso-
      }
    }  
    Fe[1][1]=Fe[1][1]/fac[q]; //uni-axis
    //for(unsigned int i=0;i<dim;i++) Fe[i][i]=Fe[i][i]/fac_T[q];
    
    //E
    Table<2, Sacado::Fad::DFad<double> > E (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				E[i][j] = -0.5*(i==j);
				for (unsigned int k=0; k<dim; ++k){
		  		E[i][j] += 0.5*Fe[k][i]*Fe[k][j];
				}
      }
      //E[i][i]+=-Ec;
    }
   // E[2][2]=0.5*(defMap.F[q][2][2]/fac[q]*defMap.F[q][2][2]/fac[q]);

    
    
    
    //S
    //layer 1 -
    Table<2, Sacado::Fad::DFad<double> > S (dim, dim);
    if (X[1]<60. or X[1]>83.) {
	double youngsModulus1;
      if (X[1]<60) youngsModulus1 = 5.93e-3;
      if (X[1]>83) youngsModulus1 = 8.88e-3; //E_1=1.5Gpa (NMC)
      if(dim==3){
        for (unsigned int i=0; i<dim; ++i){
		  for (unsigned int j=0; j<dim; ++j){
		    S[i][j]=0;
		    for (unsigned int k=0; k<dim; ++k){
		      for (unsigned int l=0; l<dim; ++l){
		        S[i][j] += SVK3D<dim>(i, j, k, l, youngsModulus1)*E[k][l];
		      }
		    }
		  }
		}
      }
      else if(dim==2){
        S[0][0]=SVK2D<dim>(0,0,youngsModulus1)*E[0][0]+SVK2D<dim>(0,1, youngsModulus1)*E[1][1]+SVK2D<dim>(0,2,youngsModulus1)*(E[0][1]+E[1][0]);
        S[1][1]=SVK2D<dim>(1,0,  youngsModulus1)*E[0][0]+SVK2D<dim>(1,1, youngsModulus1)*E[1][1]+SVK2D<dim>(1,2, youngsModulus1)*(E[0][1]+E[1][0]);
        S[0][1]=S[1][0]=SVK2D<dim>(2,0, youngsModulus1)*E[0][0]+SVK2D<dim>(2,1, youngsModulus1)*E[1][1]+SVK2D<dim>(2,2, youngsModulus1)*(E[0][1]+E[1][0]);
      }
      else if(dim==1){
        S[0][0]=SVK1D<dim>(youngsModulus1)*E[0][0];
      }
    }
    //layer 2 sep
    else{
      double youngsModulus2 = youngsModulus*1.0; //E_2=1.0GPa
      if(dim==3){
        for (unsigned int i=0; i<dim; ++i){
		  for (unsigned int j=0; j<dim; ++j){
		    S[i][j]=0;
		    for (unsigned int k=0; k<dim; ++k){
		      for (unsigned int l=0; l<dim; ++l){
		        S[i][j] += SVK3D<dim>(i, j, k, l, youngsModulus2)*E[k][l];
		      }
		    }
		  }
		}
      }
      else if(dim==2){
        S[0][0]=SVK2D<dim>(0,0,youngsModulus2)*E[0][0]+SVK2D<dim>(0,1, youngsModulus2)*E[1][1]+SVK2D<dim>(0,2,youngsModulus2)*(E[0][1]+E[1][0]);
        S[1][1]=SVK2D<dim>(1,0,  youngsModulus2)*E[0][0]+SVK2D<dim>(1,1, youngsModulus2)*E[1][1]+SVK2D<dim>(1,2, youngsModulus2)*(E[0][1]+E[1][0]);
        S[0][1]=S[1][0]=SVK2D<dim>(2,0, youngsModulus2)*E[0][0]+SVK2D<dim>(2,1, youngsModulus2)*E[1][1]+SVK2D<dim>(2,2, youngsModulus2)*(E[0][1]+E[1][0]);
      }
      else if(dim==1){
        S[0][0]=SVK1D<dim>(youngsModulus2)*E[0][0];
      }
    }
  
    //else throw "dim not equal to 1/2/3";
    

    //P
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				P[q][i][j]=0;
				for (unsigned int k=0; k<dim; ++k){
	 	 			P[q][i][j]+=Fe[i][k]*S[k][j];
				}
      }
    }
    Table<2, Sacado::Fad::DFad<double> > Sigma (dim, dim);
    //Sigma
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
 				Sigma[i][j]=0;
				for (unsigned int k=0; k<dim; ++k){
					Sigma[i][j]+=P[q][i][k]*defMap.F[q][j][k];
				}
				Sigma[i][j]/=defMap.detF[q];
      }
    }
    //strain[cell][q]=E[0][0].val()+E[1][1].val(); //trace(E)
    //strain[cell][q]= defMap.detF[q].val();
  }
}

//Mechanics residual implementation
template <int dim>
void residualForMechanics(Table<1, double >& Pyy, const FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, Table<1, Sacado::Fad::DFad<double> >& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double> >& R, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, double youngsModulus, typename hp::DoFHandler<dim>::active_cell_iterator& cell, dealii::Table<1,double>& c1_conv, dealii::Table<1, Sacado::Fad::DFad<double> >& c1, dealii::Table<1, Sacado::Fad::DFad<double> >& c1_0, dealii::Table<1, Sacado::Fad::DFad<double> >& Temp, double currentTime, double totalTime, double dt){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //Temporary arrays
  Table<3,Sacado::Fad::DFad<double> > P (n_q_points, dim, dim);
  //evaluate mechanics=====added 
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, ULocal, P, defMap, youngsModulus, cell, c1_conv, c1, c1_0,Temp, dt, currentTime);
  //evaluate Residual

  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<dim){
      // R = Grad(w)*P
      for (unsigned int q=0; q<n_q_points; ++q){
				for (unsigned int d = 0; d < dim; d++){
	  			R[i] +=  fe_values.shape_grad(i, q)[d]*P[q][ck][d]*fe_values.JxW(q);
				}
      }
			//R[i] = R[i]*1.0e-10;
    }
  }
	for(unsigned int q=0;q<n_q_points;q++){
		const Point<dim> X = fe_values.quadrature_point(q);
		Pyy[q]=P[q][1][1].val();
	}
//std::cout<<"Pyy[q]=="<<Pyy[0]<<std::endl;

  //input Force boundary condition
  //double t=currentTime/totalTime;
  //double force= -1.47; //10N 
  //double force= -1e-4; //unit:N/um^2, 100g 
  //double force= 0;
  //Neumann conditions
//Flux like tracking bounday coding
/*	
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell->face(faceID)->at_boundary()){
      const Point<dim> face_center = cell->face(faceID)->center();
      if (face_center[1] == 128){ //flux on Y=1 boundary
				fe_face_values.reinit (cell, faceID);
				for (unsigned int i=0; i<dofs_per_cell; ++i) {
					const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	  		 	if (ck==1){
						for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
	   					R[i] += -fe_face_values.shape_value(i, q)*force*fe_face_values.JxW(q); 
	   		 		 }
	  	 		}
		 		}
    	}
   	 	else if (face_center[1] == 200.0){ //flux on Y=0 boundary
		      fe_face_values.reinit (cell, faceID);
		      for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  			const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	  			if (ck==1){
						for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q){
	      			R[i] += fe_face_values.shape_value(i, q)*force*fe_face_values.JxW(q); 
						}
				}
		      }
		}
   	}
  }
*/	
}

#endif /* MECHANICS_H_ */
