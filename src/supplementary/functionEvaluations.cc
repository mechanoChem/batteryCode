#include"../../include/supplementary/functionEvaluations.h"

template <class T, int dim>
void evaluateScalarFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;

	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		U[q]=0.0; //U
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				U[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
			}
		}
	}
}
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);


template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);


template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, bool gradientInCurrentConfiguration, deformationMap<T, dim>& defMap)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	Table<1, T> refGradU(dim);
	for (unsigned int i=0;i<dim;i++) refGradU[i]=0;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		//refGradU=0;
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
		//Transform gradient to current configuration. gradW=(F^-T)*GradW
		for (unsigned int i=0; i<dim; ++i){
			if (gradientInCurrentConfiguration==true)
				gradU[q][i]=0.0;
				for (unsigned int j=0; j<dim; ++j){
					gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
			}
		}
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU,bool gradientInCurrentConfiguration, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU,bool gradientInCurrentConfiguration, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, bool gradientInCurrentConfiguration, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU,bool gradientInCurrentConfiguration, deformationMap<double, 1>& defMap);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU,bool gradientInCurrentConfiguration, deformationMap<double, 2>& defMap);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU,bool gradientInCurrentConfiguration, deformationMap<double, 3>& defMap);



template <class T, int dim>
void evaluateVectorFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& U){
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int k=0; k<dim; ++k) U[q][k]=0; 
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q); //U
			}
		}
	}
}
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);
template void evaluateVectorFunction<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);
template void evaluateVectorFunction<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);



template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU){
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			for(unsigned int j=0;j<dim;++j){
					gradU[q][i][j]=0.0; //gradU
			}
		}
	}
	//gradU.fill(0.0);
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
	}
}
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);



template <class T, int dim>
  void getDeformationMap(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap,unsigned int iteration){
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
template void getDeformationMap<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap, unsigned int iteration);
template void getDeformationMap<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 1>& defMap, unsigned int iteration);
template void getDeformationMap<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 2>& defMap, unsigned int iteration);
template void getDeformationMap<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 3>& defMap, unsigned int iteration);

