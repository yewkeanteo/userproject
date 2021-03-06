#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
/*#include "PlaneStimulusCellFactory.hpp"*/

#include "LuoRudy1991.hpp"
#include "SimpleStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TetrahedralMesh.hpp"
#include <cmath>

class CustomCellFactory : public AbstractCardiacCellFactory<2> // <3> here
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    CustomCellFactory()
        : AbstractCardiacCellFactory<2>(), // <3> here as well!
          mpStimulus(new SimpleStimulus(-5e5, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];

        if ((x<0.1+1e-6) && (y<0.1+1e-6) )
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

class TestBidomainWithBathTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestWithBathAndElectrodes()
    {
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("CellFactoryTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CellFactoryTutorialResults");
		HeartConfig::Instance()->SetVisualizeWithVtk(true);
		
		HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.01, 0.1);

		/*Generate a Mesh Here*/
		DistributedTetrahedralMesh<2,2> mesh;
        double h=0.02;
        mesh.ConstructRegularSlabMesh(h, 0.8 /*length*/, 0.3 /*width*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        /*PlaneStimulusCellFactory<CellToRORddynClmidFromCellML,2> cell_factory(0.0);

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);*/
		
		
        
		std::set<unsigned> tissue_ids;
        static unsigned tissue_id=0;
        tissue_ids.insert(tissue_id);

        std::set<unsigned> bath_ids;
        static unsigned bath_id1=1;
        bath_ids.insert(bath_id1);
        static unsigned bath_id2=2;
        bath_ids.insert(bath_id2);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            double x = iter->CalculateCentroid()[0];
            double y = iter->CalculateCentroid()[1];
            if (sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02)
            {
                if (y<0.05)
                {
                    //Outside circle on the bottom
                    iter->SetAttribute(bath_id1);
                }
                else
                {
                    //Outside circle on the top
                    iter->SetAttribute(bath_id2);
                }
            }
            else
            {
                //IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
            }
        }

        mesh.SetMeshHasChangedSinceLoading();

        HeartConfig::Instance()->SetBathConductivity(7.0);  //bath_id1 tags will take the default value (actually 7.0 is the default)
        std::map<unsigned, double> multiple_bath_conductivities;
        multiple_bath_conductivities[bath_id2] = 6.5;  // mS/cm
		
        HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);
		
        // For default conductivities and explicit cell model -1e4 is under threshold, -1.4e4 too high - crashes the cell model
        // For heterogeneous conductivities as given, -1e4 is under threshold
        double magnitude = -14.0e3; // uA/cm^2
        double start_time = 0.0;
        double duration = 1; //ms
		
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);
		
		
		CustomCellFactory cell_factory;
		
		/*Create the problem class*/
        BidomainProblem<2> bidomain_problem( &cell_factory, true);

        bidomain_problem.SetMesh(&mesh);
		
		bool partial_output = true;
        if (partial_output)
        {
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(0);
            nodes_to_be_output.push_back((unsigned)round( (mesh.GetNumNodes()-1)/2 ));
            nodes_to_be_output.push_back(mesh.GetNumNodes()-1);
            bidomain_problem.SetOutputNodes(nodes_to_be_output);
        }

		bidomain_problem.SetWriteInfo();
		 
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        Vec solution = bidomain_problem.GetSolution(); // the Vs and phi_e's, as a PetSc vector
        ReplicatableVector solution_repl(solution);
		
		HeartEventHandler::Headings();
		HeartEventHandler::Headings();
		
        bool ap_triggered = false;
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            if (HeartRegionCode::IsRegionTissue( iter->GetRegion() ))
            {
                if (solution_repl[2*iter->GetIndex()] > 0.0) // 2*i, ie the voltage for this node (would be 2*i+1 for phi_e for this node)
                {
                    ap_triggered = true;
                }
            }
        }
        TS_ASSERT(ap_triggered);
    }
};