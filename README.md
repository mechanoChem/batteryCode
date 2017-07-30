openIGA code: modelling coupled electrostatic, chemical, thermal and mechanical processes in battery materials

List of contributors:

Zhenlin Wang (Lead Developer)

Krishna Garikipati


Overview
====================================================================

batteryCode is based on deal.ii (www.dealii.org). It consists of four four modules. 1. battey_base provide generic finite element code; 2. model consists of mechanics at finite strain, reaction-diffusion equations and Possion equations; 3. ElectricChemoExpression provides electric-chemical formulations e.g. Butler-Volmer equations. 4. Application consists of specific application.

To run one application please switch to the correspding branch. Then goes to application/

Currently we have:

1.homogenzed model 
 Note: original_code is old code for paper Z. Wang, J. Siegel, K. Garikipati, “Intercalation driven porosity effects on the electro-chemo-thermo-mechanical response in continuum models for battery material electrodes”. 
2.particle model (coming soon)

License
====================================================================
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.

Acknowledgements
The work was funded in part by Advanced Research
Projects AgencyEnergy (ARPA-E), U.S. Department of Energy, under Award #DE-AR0000269. The computations have been carried out on the Flux computing cluster at University of Michigan, using hardware resources supported by the U.S. Department of Energy, Oce of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan.

Installation with cmake:
====================================================================
Install pre-required libs:
Install CMake [http://www.cmake.org/download/]

Install deal.II (version 8.4.1 and up) with Trilinos and Petsc [www.dealii.org/download.html] Deal.II OSX binaries include full packages of deal.ii with Trillions, Petsc and other useful libs.

Install Trillions [https://trilinos.org/]

Install Petsc [https://www.mcs.anl.gov/petsc/]

$ cmake CMakeLists.txt

$ make install

$ make run


Reference
====================================================================
If you write a paer using results obtained with the help of this code, please consider citing one or more of the following:

For use of homogenized model:

Z. Wang, J. Siegel, K. Garikipati, “Intercalation driven porosity effects on the electro-chemo-thermo-mechanical response in continuum models for battery material electrodes”, To appear in Journal of the Electrochemical Society, [https://arxiv.org/abs/1609.08944]


data:08/01/2017