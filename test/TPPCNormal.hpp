#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <cassert>
#include "PropagationPropertiesCalculator.hpp"

class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:
     void TestConductionBidomain3D()
	 {
		unsigned middle_index = 218;
        unsigned rhs_index = 222;
		
		Hdf5DataReader simulation_data("projects/userproject/Results/1000ms/NormalConductance2",
                                       "NormalConductance2", false);

        PropagationPropertiesCalculator ppc(&simulation_data);
		//ppc.CalculateConductionVelocity(middle_index,rhs_index,0.08);
		std::cout << "The normal conduction velocity is "<<ppc.CalculateConductionVelocity(middle_index,rhs_index,0.08) << "\n";
	 }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_