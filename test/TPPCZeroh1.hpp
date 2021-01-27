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
		unsigned lhs_index = 838;
        unsigned rhs_index = 842;
		unsigned lhs_2_index = 836;
        unsigned rhs_2_index = 844;
		unsigned lhs_3_index = 834;
        unsigned rhs_3_index = 846;
		unsigned lhs_4_index = 832;
        unsigned rhs_4_index = 848;
		unsigned lhs_5_index = 830;
        unsigned rhs_5_index = 850;
		Hdf5DataReader simulation_data("testoutput/ZeroConductance1100h1","ZeroConductance1100h1", false);

        PropagationPropertiesCalculator ppc(&simulation_data);
		//ppc.CalculateConductionVelocity(lhs_index,rhs_index,0.08);
		std::cout << "The normal conduction velocity between Node 219 and Node 221 is "<<ppc.CalculateConductionVelocity(lhs_index,rhs_index,0.04) << "\n";
		std::cout << "The normal conduction velocity between Node 218 and Node 222 is "<<ppc.CalculateConductionVelocity(lhs_2_index,rhs_2_index,0.08) << "\n";
		std::cout << "The normal conduction velocity between Node 217 and Node 223 is "<<ppc.CalculateConductionVelocity(lhs_3_index,rhs_3_index,0.12) << "\n";
		std::cout << "The normal conduction velocity between Node 216 and Node 224 is "<<ppc.CalculateConductionVelocity(lhs_4_index,rhs_4_index,0.16) << "\n";
		std::cout << "The normal conduction velocity between Node 215 and Node 225 is "<<ppc.CalculateConductionVelocity(lhs_5_index,rhs_5_index,0.20) << "\n";
	 }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_