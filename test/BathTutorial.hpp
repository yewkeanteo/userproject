#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"

#include "ToRORddynClmid.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include <cmath>

class TestBidomainWithBathTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestWithBathAndElectrodes()
    {
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BathTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BathResults");
		HeartConfig::Instance()->SetVisualizeWithVtk(true);
		
        HeartConfig::Instance()->SetOdeTimeStep(0.0001);  //ms

        PlaneStimulusCellFactory<CellToRORddynClmidFromCellML,2> cell_factory(0.0);

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

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
        double magnitude = -15.5e3; // uA/cm^2
        double start_time = 0.0;
        double duration = 1; //ms
		
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);
		
        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        bidomain_problem.SetMesh(&mesh);

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