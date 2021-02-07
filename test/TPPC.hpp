#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <cassert>
#include "PropagationPropertiesCalculator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "Node.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:
     void TestConductionBidomain3D()
	 {
		unsigned middle_index = 425;
        unsigned rhs_index = 435;
		
		Hdf5DataReader simulation_data("projects/userproject/Results/DoubleCell/NormalDoubleCell500",
                                       "NormalDoubleCell500", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
		//ppc.CalculateConductionVelocity(middle_index,rhs_index,0.1);
		std::cout << "The conduction velocity is "<<ppc.CalculateConductionVelocity(middle_index,rhs_index,0.2) << "\n";
	 }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_