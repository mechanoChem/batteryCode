/*
 * quadraturePointData.h
 *
 *  Created on: Jun 9, 2011
 *      Author: rudraa
 */

#ifndef QUADRATUREPOINTDATA_H_
#define QUADRATUREPOINTDATA_H_
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include "supplementaryFunctions.h"

template <int dim, int degree>
class quadPointProjection
{
	public:
		quadPointProjection (): objectInitializationCount(100000) {};
		~quadPointProjection(){};
		static void reinit(Triangulation<dim>& triangulation_temp, QGauss<dim>& quadrature_formula_temp){
			triangulation=&triangulation_temp; quadrature_formula=quadrature_formula_temp;
			dof_handler.clear(); dof_handler.initialize(*triangulation, fe);
			//generate sparsity pattern
			hanging_node_constraints.clear (); DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);  hanging_node_constraints.close ();
			CompressedSparsityPattern c_sparsity(dof_handler.n_dofs()); DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
			hanging_node_constraints.condense (c_sparsity); sparsity_pattern.copy_from(c_sparsity);
			//initialize quadratureData structures
			mass_matrix.reinit(sparsity_pattern);
			MatrixCreator::create_mass_matrix(dof_handler, quadrature_formula, mass_matrix);
			hanging_node_constraints.condense (mass_matrix);
			classInitializationCount++;
		}
		void reset(){
			if (classInitializationCount==0){printf("quadPointProjection class static members not initialized"); throw "quadPointProjection class static members not initialized";}
			quadratureData.clear(); gradientQuadratureData.clear();
			typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
			for (;cell!=endc; ++cell){
				std::vector<double> temp1(quadrature_formula.size(), 0.0); Table<2, double> temp2(quadrature_formula.size(), dim); temp2.fill(0.0);
				quadratureData.insert(std::make_pair(cell, temp1));
				gradientQuadratureData.insert(std::make_pair(cell, temp2));
			}
			system_rhs.reinit (dof_handler.n_dofs()); nodalValues.reinit (dof_handler.n_dofs()); nodalValuesGradX.reinit (dof_handler.n_dofs()); nodalValuesGradY.reinit (dof_handler.n_dofs());
			objectInitializationCount=classInitializationCount;
		}
		std::vector<double>& operator[](typename DoFHandler<dim>::active_cell_iterator& cell){
			if(objectInitializationCount!=classInitializationCount){throw("stress not reset");}
			return quadratureData[cell];
		}
		double gradient(typename DoFHandler<dim>::active_cell_iterator& cell, unsigned int q, unsigned int component){
			if(objectInitializationCount!=classInitializationCount){throw("stress not reset");}
			return gradientQuadratureData[cell][q][component];
		}
		double operator[](unsigned int index){
			return nodalValues(index);
		}
		void project(){
			if(objectInitializationCount!=classInitializationCount){reset();}
			//Temporary variables
			const unsigned int   dofs_per_cell = fe.dofs_per_cell;
			const unsigned int   num_quad_points    = quadrature_formula.size();
			Vector<double>       local_rhs (dofs_per_cell);
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
			FEValues<dim> fe_values (fe, quadrature_formula, update_values   | update_gradients | update_quadrature_points | update_JxW_values);

			//Assembly
			system_rhs=0;
			typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
			for (;cell!=endc; ++cell){
				fe_values.reinit (cell);
				local_rhs = 0;
				cell->get_dof_indices (local_dof_indices);

				//Looping over quadrature points
				for (unsigned int q=0; q<num_quad_points; ++q){
					for (unsigned int i=0; i<dofs_per_cell; ++i){
						local_rhs(i) +=  (fe_values.shape_value(i,q)*quadratureData[cell][q])*fe_values.JxW(q);
					}
				}
				hanging_node_constraints.distribute_local_to_global (local_rhs, local_dof_indices, system_rhs);
			}
			nodalValues=0; solveSystem(mass_matrix, system_rhs, nodalValues, 1.0e-15);
			hanging_node_constraints.distribute (nodalValues);

			//Compute Gradients
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			for (;cell!=endc; ++cell){
				fe_values.reinit (cell);
				cell->get_dof_indices (local_dof_indices);
				//Looping over quadrature points
				for (unsigned int q=0; q<num_quad_points; ++q){
					for (unsigned i=0; i<dim; ++i){
						gradientQuadratureData[cell][q][i]=0.0;
						for (unsigned int j=0; j<dofs_per_cell; ++j){
							gradientQuadratureData[cell][q][i]+=fe_values.shape_grad(j,q)[i]*nodalValues(local_dof_indices[j]);
						}
					}
				}
			}

			//Assembly GradX
			system_rhs=0;
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			for (;cell!=endc; ++cell){
				fe_values.reinit (cell);
				local_rhs = 0;
				cell->get_dof_indices (local_dof_indices);

				//Looping over quadrature points
				for (unsigned int q=0; q<num_quad_points; ++q){
					for (unsigned int i=0; i<dofs_per_cell; ++i){
						local_rhs(i) +=  (fe_values.shape_value(i,q)*gradientQuadratureData[cell][q][0])*fe_values.JxW(q);
					}
				}
				hanging_node_constraints.distribute_local_to_global (local_rhs, local_dof_indices, system_rhs);
			}
			nodalValuesGradX=0; solveSystem(mass_matrix, system_rhs, nodalValuesGradX, 1.0e-15);
			hanging_node_constraints.distribute (nodalValuesGradX);

			//Assembly GradY
			system_rhs=0;
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			for (;cell!=endc; ++cell){
				fe_values.reinit (cell);
				local_rhs = 0;
				cell->get_dof_indices (local_dof_indices);

				//Looping over quadrature points
				for (unsigned int q=0; q<num_quad_points; ++q){
					for (unsigned int i=0; i<dofs_per_cell; ++i){
						local_rhs(i) +=  (fe_values.shape_value(i,q)*gradientQuadratureData[cell][q][1])*fe_values.JxW(q);
					}
				}
				hanging_node_constraints.distribute_local_to_global (local_rhs, local_dof_indices, system_rhs);
			}
			nodalValuesGradY=0; solveSystem(mass_matrix, system_rhs, nodalValuesGradY, 1.0e-15);
			hanging_node_constraints.distribute (nodalValuesGradY);

		}
		void output(unsigned int cycle) const{
		  //Write results to VTK file
		  char filename [200]; sprintf (filename, "/home/zhenlin/workspace/results/batteryProject/%s_Projections-%u.vtk", outputFilePrefix.c_str(), cycle); std::ofstream output (filename);
		  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);
		  //Add nodal DOF data
		  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation(1, DataComponentInterpretation::component_is_scalar);
		  data_out.add_data_vector (nodalValues, "Jn" , DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
		  //data_out.add_data_vector (nodalValuesGradX, "gradCauchyStressX" , DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
		  //data_out.add_data_vector (nodalValuesGradY, "gradCauchyStressY" , DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
		  data_out.build_patches(); data_out.write_vtk (output); output.close();
		}
		static std::string outputFilePrefix;
	private:
		unsigned int objectInitializationCount;
		static unsigned int classInitializationCount;
		static Triangulation<dim>* triangulation;
		static FESystem<dim> fe;
		static FE_Q<dim> finiteElement;
		static DoFHandler<dim> dof_handler;
		static SparsityPattern sparsity_pattern;
		static QGauss<dim> quadrature_formula;
		static ConstraintMatrix hanging_node_constraints;
		std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector<double> > quadratureData;
		std::map<typename DoFHandler<dim>::active_cell_iterator, Table<2, double> > gradientQuadratureData;
		Vector<double> system_rhs;
		static SparseMatrix<double> mass_matrix;
		Vector<double> nodalValues, nodalValuesGradX, nodalValuesGradY;
};

template <int dim, int degree>
std::string quadPointProjection<dim, degree>::outputFilePrefix="default";

template <int dim, int degree>
unsigned int quadPointProjection<dim, degree>::classInitializationCount(0);

template <int dim, int degree>
DoFHandler<dim> quadPointProjection<dim, degree>::dof_handler;

template <int dim, int degree>
Triangulation<dim>* quadPointProjection<dim, degree>::triangulation=0;

template <int dim, int degree>
QGauss<dim> quadPointProjection<dim, degree>::quadrature_formula(2);

template <int dim, int degree>
SparsityPattern quadPointProjection<dim, degree>::sparsity_pattern;

template <int dim, int degree>
SparseMatrix<double> quadPointProjection<dim, degree>::mass_matrix;

template <int dim, int degree>
ConstraintMatrix quadPointProjection<dim, degree>::hanging_node_constraints;

template <int dim, int degree>
FE_Q<dim> quadPointProjection<dim, degree>::finiteElement(degree);

template <int dim, int degree>
FESystem<dim>  quadPointProjection<dim, degree>::fe(finiteElement, 1);

#endif /* QUADRATUREPOINTDATA_H_ */

