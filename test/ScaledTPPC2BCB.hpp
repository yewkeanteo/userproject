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
		unsigned lhs_index = 262;
        unsigned rhs_index = 270;
		//unsigned lhs_2_index = 424;
        //unsigned rhs_2_index = 434;
		//unsigned lhs_3_index = 425;
        //unsigned rhs_3_index = 435;
		//unsigned lhs_4_index = 426;
        //unsigned rhs_4_index = 436;
		//unsigned lhs_5_index = 427;
        //unsigned rhs_5_index = 437;
		
		Hdf5DataReader simulation_data("projects/userproject/Results/Scaled2BCB/ScaledEN2BCB1000",
                                       "ScaledEN2BCB1000", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
		//ppc.CalculateConductionVelocity(middle_index,rhs_index,0.1);
		std::cout << "The normal conduction velocity between Node 262 and Node 270 is "<<ppc.CalculateConductionVelocity(lhs_index,rhs_index,0.008) << "\n";
		//std::cout << "The normal conduction velocity between Node 424 and Node 434 is "<<ppc.CalculateConductionVelocity(lhs_2_index,rhs_2_index,0.2) << "\n";
		//std::cout << "The normal conduction velocity between Node 425 and Node 435 is "<<ppc.CalculateConductionVelocity(lhs_3_index,rhs_3_index,0.2) << "\n";
		//std::cout << "The normal conduction velocity between Node 426 and Node 436 is "<<ppc.CalculateConductionVelocity(lhs_4_index,rhs_4_index,0.2) << "\n";
		//std::cout << "The normal conduction velocity between Node 427 and Node 437 is "<<ppc.CalculateConductionVelocity(lhs_5_index,rhs_5_index,0.2) << "\n";
	 }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_