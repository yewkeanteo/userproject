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
		unsigned middle_index = 218;
        unsigned rhs_index = 222;
		
		Hdf5DataReader simulation_data("projects/userproject/Results/1000ms/NormalConductance2",
                                       "NormalConductance2", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
		//ppc.CalculateConductionVelocity(middle_index,rhs_index,0.1);
		std::cout << "The conduction velocity is "<<ppc.CalculateConductionVelocity(middle_index,rhs_index,0.1) << "\n";
		ChastePoint<1> point1(0.16);
     	ChastePoint<1> point2(0.2);
		Node<1> node1(0, point1);
		TS_ASSERT_EQUALS(node1.GetIndex(), 0u);
		//n.
	 }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_