#include"../../include/model.h"

//default evaluateStress No thermal chemical coupling

template <class T, int dim>
void model<T,dim>::evaluateStrain(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<3> indices(n_q_points, dim, dim);
	Fe.reinit(indices);
	E.reinit(indices);
	
	

	for (unsigned int q=0; q<n_q_points; ++q){   
		for (unsigned int i=0; i<dim; ++i){
    	for (unsigned int j=0; j<dim; ++j){
				Fe[q][i][j]=0.0;
				E[q][i][j]=0.0;
				
	  		Fe[q][i][j]=defMap.F[q][i][j];
	  		Fe[q][i][j]=defMap.F[q][i][j]; //iso-
    	}
  	}
		//E
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				E[q][i][j] = -0.5*(i==j);
				for (unsigned int k=0; k<dim; ++k){
					E[q][i][j] += 0.5*Fe[q][k][i]*Fe[q][k][j];
				}
			}	
		}
	}
}

template <class T, int dim>
void model<T,dim>::evaluateStress(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap){
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<3> indices(n_q_points, dim, dim);
	P.reinit(indices);

	evaluateStrain(fe_values, defMap);
  for (unsigned int q=0; q<n_q_points; ++q){
  	//S
  	Table<2, T > S (dim, dim);
		const Point<dim> X = fe_values.quadrature_point(q);
		double youngsModulus;
  	if (X[1]<l_neg) youngsModulus=youngModule_neg;
		else if (X[1]<l_neg+l_sep) youngsModulus=youngModule_sep;
		else if (X[1]<l_neg+l_sep+l_pos) youngsModulus=youngModule_pos;
  	if(dim==3){
			for (unsigned int i=0; i<dim; ++i){
		  	for (unsigned int j=0; j<dim; ++j){
		    	S[i][j]=0;
		    	for (unsigned int k=0; k<dim; ++k){
		      	for (unsigned int l=0; l<dim; ++l){
		        	S[i][j] += SVK3D(i, j, k, l, youngsModulus)*E[q][k][l];
		      	}
		    	}
		  	}
			}
		}
  	else if(dim==2){
			S[0][0]=SVK2D(0,0, youngsModulus)*E[q][0][0]+SVK2D(0,1, youngsModulus)*E[q][1][1]+SVK2D(0,2, youngsModulus)*(E[q][0][1]+E[q][1][0]);
    	S[1][1]=SVK2D(1,0, youngsModulus)*E[q][0][0]+SVK2D(1,1, youngsModulus)*E[q][1][1]+SVK2D(1,2, youngsModulus)*(E[q][0][1]+E[q][1][0]);
    	S[0][1]=S[1][0]=SVK2D(2,0, youngsModulus)*E[q][0][0]+SVK2D(2,1, youngsModulus)*E[q][1][1]+SVK2D(2,2, youngsModulus)*(E[q][0][1]+E[q][1][0]);
		}
		else if(dim==1){
   	 S[0][0]=youngsModulus*E[q][0][0];
  	}
		//P
  	for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				P[q][i][j]=0.0;
				for (unsigned int k=0; k<dim; ++k){
					P[q][i][j]+=Fe[q][i][k]*S[k][j];
				}
			}
		}
	}
}
template <class T, int dim>
void model<T,dim>::residualForMechanics(const FEValues<dim>& fe_values, unsigned int DOF, deformationMap<T, dim>& defMap, Table<1, T >& R)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //Temporary arrays
  evaluateStress(fe_values, defMap);
  //evaluate ResidualULocal
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<dim){
      for (unsigned int q=0; q<n_q_points; ++q){
				for (unsigned int d = 0; d < dim; d++){
	  			R[i] +=  fe_values.shape_grad(i, q)[d]*P[q][ck][d]*fe_values.JxW(q);
				}
      }
  	}
  }
}

template <class T, int dim>
void model<T,dim>::residualForNewmmanBC(Table<1, T >& R){
	
}

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;
