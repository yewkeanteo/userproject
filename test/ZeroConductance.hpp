#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"

#include "ToRORddynClmidCvode.hpp" // Imports the native CVODE version of the cell model.
#include "SimpleStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TetrahedralMesh.hpp"
#include <cmath>
#ifdef CHASTE_CVODE

class CustomCellFactory : public AbstractCardiacCellFactory<2> // <2> for 2D simulation.
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    CustomCellFactory() // Defining a cell factory as input
        : AbstractCardiacCellFactory<2>(), // <2> for 2D simulation.
          mpStimulus(new SimpleStimulus(-0, 2))
    {
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        AbstractCvodeCell* p_cell;
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];

        if ((x<0.1+1e-6) && (y<0.1+1e-6)) //Zero Stimulus to the cell as the cell is sto be stimulated from through the bath.
        {
            p_cell = new CellToRORddynClmidFromCellMLCvode(p_empty_solver, mpZeroStimulus);
        }
        else
        {
            p_cell = new CellToRORddynClmidFromCellMLCvode(p_empty_solver, mpZeroStimulus);
        }
        p_cell->SetTolerances(1e-5,1e-7);
		if (x<0.2)
		{
			//No change to conductance of cell factory (left side)
		}
		else
		{
			//Change conductance of cell factory (right side)
			p_cell->SetParameter("membrane_fast_sodium_current_conductance", 0);
		}
        return p_cell;
    }
};

#endif // CHASTE_CVODE

class TestBidomainWithBathTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestWithBathAndElectrodes()
    {
		/*Generate a Mesh Here*/
		DistributedTetrahedralMesh<2,2> mesh;
        double h=0.02;
        mesh.ConstructRegularSlabMesh(h, 0.4 /*length*/, 0.4 /*width*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
		
        HeartConfig::Instance()->SetSimulationDuration(1000.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("ZeroConductance2");
        HeartConfig::Instance()->SetOutputFilenamePrefix("ZeroConductance2");
		HeartConfig::Instance()->SetVisualizeWithVtk(true);
		
		HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.0001, 0.001); 
		// CVODE will take as many adaptive internal timesteps as it requires each time it is called, 
		// so we should just call it once per PDE timestep - i.e. set the ODE and PDE timesteps to be the same.
		
		// Original bidomain bath code
		std::set<unsigned> tissue_ids;
        static unsigned tissue_id=0;
        tissue_ids.insert(tissue_id);

        std::set<unsigned> bath_ids;
        static unsigned bath_id1=1;
        bath_ids.insert(bath_id1);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);
		// from the center, we set those more than 0.1 to be the bath elements
        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            double x = iter->CalculateCentroid()[0];
            double y = iter->CalculateCentroid()[1];
            if (sqrt((x-0.2)*(x-0.2) + (y-0.2)*(y-0.2)) > 0.1)
            {
                    //Outside circle on the bottom
                    iter->SetAttribute(bath_id1);
            }
            else
            {
                //IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
            }
        }

        mesh.SetMeshHasChangedSinceLoading(); // Inform Chaste of the fact that we have modified the mesh by setting element attributes.

        HeartConfig::Instance()->SetBathConductivity(7.0);  //bath_id1 tags will take the default value (actually 7.0 is the default)
		
        //HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);

        double magnitude = -20.0e3; // uA/cm^2
        double start_time = 1.0;
        double duration = 100; //ms
		
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);
		
		CustomCellFactory cell_factory;
		
		/*Create the problem class*/
        BidomainProblem<2> bidomain_problem( &cell_factory, true);  // Create problem class with pointer to cell factory and pass true to indicate we are solving it
        bidomain_problem.SetMesh(&mesh); // Sets mesh and electrodes
		bidomain_problem.SetWriteInfo(); 
		 
        bidomain_problem.Initialise(); // Initialise and Solve
        bidomain_problem.Solve();

        Vec solution = bidomain_problem.GetSolution(); // the Vs and phi_e's, as a PetSc vector
        ReplicatableVector solution_repl(solution);
		
		HeartEventHandler::Headings();
		HeartEventHandler::Headings();
		
        /*bool ap_triggered = false;
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
		*/
    }
};