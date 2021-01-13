#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "ToRORddynClmid.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class BenchmarkCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BenchmarkCellFactory()
        : AbstractCardiacCellFactory<3>(), // <3> here as well!
          mpStimulus(new SimpleStimulus(-100000.0, 2))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        if ((x<0.1+1e-6) && (y<0.1+1e-6) && (z<0.1+1e-6))
        {
            return new CellToRORddynClmidFromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellToRORddynClmidFromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

class TestMonodomain3dExampleTutorial : public CxxTest::TestSuite
{
public:
    void TestMonodomain3d()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        double h=0.02;
        mesh.ConstructRegularSlabMesh(h, 0.8 /*length*/, 0.3 /*width*/, 0.3 /*depth*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(5); //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3dExample");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.01, 0.1);

        BenchmarkCellFactory cell_factory;

        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.SetMesh(&mesh);

        bool partial_output = true;
        if (partial_output)
        {
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(0);
            nodes_to_be_output.push_back((unsigned)round( (mesh.GetNumNodes()-1)/2 ));
            nodes_to_be_output.push_back(mesh.GetNumNodes()-1);
            monodomain_problem.SetOutputNodes(nodes_to_be_output);
        }

        monodomain_problem.SetWriteInfo();

        monodomain_problem.Initialise();                                                                                         
        monodomain_problem.Solve();

        /*ReplicatableVector voltage(monodomain_problem.GetSolution());
        TS_ASSERT_DELTA(voltage[0], 34.9032, 1e-2); */
    }
};