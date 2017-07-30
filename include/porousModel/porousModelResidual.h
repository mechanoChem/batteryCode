#ifndef porpuseModelResidual_h
#define porpuseModelResidual_h

#include "../model.h"

template <class T, int dim>
class porousModelResidual : public model<T, dim>
{
 public:
   porousModelResidual ();
   ~porousModelResidual();
	 
	 /*
	 *chemo
	 */
	 
	 /*
	 *chemical residual
	 */ 
 	 void residualForChemoCli(const FEValues<dim>& fe_values,  unsigned int DOF, dealii::Table<1, T >& R,deformationMap<T, dim>& defMap,deformationMap<double, dim>& defMapConv, dealii::Table<1, T >& c1,dealii::Table<1,double>& c1_conv, dealii::Table<1, T >& jn,dealii::Table<1, T >& eps_s, dealii::Table<1,double>& eps_s_conv);
	 void residualForChemoCli_plus(const FEValues<dim>& fe_values,  unsigned int DOF, dealii::Table<1, T >& R, deformationMap<T, dim>& defMap, deformationMap<double, dim>& defMapConv,dealii::Table<1, T >& c2, dealii::Table<1,double>& c2_conv, dealii::Table<2, T >& c2_j, dealii::Table<1, T >& jn, dealii::Table<1, T >& eps_l, dealii::Table<1,double>& eps_l_conv,dealii::Table<1,T > D_l);
	 void residualForChemoPhis(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, dealii::Table<1,T >& phis, dealii::Table<2,T >& phi_s_grad, dealii::Table<1,T >& jn, dealii::Table<1,T >& eps_s);
	 void residualForChemoPhie(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R,  deformationMap<T, dim>& defMap, dealii::Table<1,T >& c2, dealii::Table<2,T >& c2_j, dealii::Table<1,T >& Temp, dealii::Table<1,T >& phie, dealii::Table<2,T >& phi_e_grad, dealii::Table<1,T >& jn, dealii::Table<1,T >& eps_l, dealii::Table<1,T >Ke);
	 void residualForPhisBC(const FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, unsigned int period, double periodTime);
   void residualForChemoThermal(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R,  deformationMap<T, dim>& defMap, dealii::Table<1, T >& Temp, dealii::Table<1,double>& T_conv, dealii::Table<2,T >& T_grad, dealii::Table<1, T >& Q);
	 void residualForHeatConv(const FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R,dealii::Table<1, T >& T_face);
	 
	 void evaluateStrain(const FEValues<dim>& fe_values, deformationMap<T, dim>& defMap);
	 void IncludeStrain_chemo(dealii::Table<1, T >& Unit_c_li);
	 void IncludeStrain_thermal(dealii::Table<1, T >& Temp);
	 
	 Table<1, T> fac_cli, fac_T; 
	 bool includeSt_chemo=false, includeSt_thermal=false;
};

#endif
