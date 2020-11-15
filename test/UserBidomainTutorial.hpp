#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "SimpleStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ToRORddynClmid.hpp"

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-5e5, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        if (x<0.02+1e-6 && y<0.02+1e-6) // ie if x<=0.02 and y<=0.02 (and we are assuming here x,y>=0).
        {
            return new CellToRORddynClmidFromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellToRORddynClmidFromCellML(mpSolver, mpZeroStimulus);
        }
    }

};

class TestRunningBidomainSimulationsTutorial : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation()
    {
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        //HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);

        PointStimulus2dCellFactory cell_factory;

        BidomainProblem<2> bidomain_problem( &cell_factory );

        // bidomain_problem.Initialise();
        // bidomain_problem.Solve();

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(6.2, 2.4));

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.001, 0.01);

        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for (unsigned i=0; i<res_repl.GetSize(); i++)
        {
        //    std::cout << res_repl[i] << "\n";
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
