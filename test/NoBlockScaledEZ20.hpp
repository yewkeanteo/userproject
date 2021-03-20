#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"

#include "ToRORddynClmidCvode.hpp"
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
        : AbstractCardiacCellFactory<2>(), 
          mpStimulus(new SimpleStimulus(-0e4, 2))
    {
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        AbstractCvodeCell* p_cell;
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
         if ((x<0.00801) && (y<0.00801))
			{
			if ((x>0.00799) && (y>0.00799))
				{
				p_cell = new CellToRORddynClmidFromCellMLCvode(p_empty_solver, mpZeroStimulus);
				}
			else
				{
				p_cell = new CellToRORddynClmidFromCellMLCvode(p_empty_solver, mpZeroStimulus);
				}
			}
        else
        {
            p_cell = new CellToRORddynClmidFromCellMLCvode(p_empty_solver, mpZeroStimulus);
        }
        p_cell->SetTolerances(1e-5,1e-7);
		
		if (x<0.0675)
		{
		}
		else
		{
			if (x>0.105)
			{
			}
			else
			{
			//Change conductance of cell factory (right side)
			//p_cell->SetParameter("membrane_fast_sodium_current_conductance", 0);
			}
		}
        return p_cell;
    }
};

class TestBidomainWithBathTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestWithBathAndElectrodes()
    {
		/*Generate a Mesh Here*/
		DistributedTetrahedralMesh<2,2> mesh;
        double h=0.0005;
        mesh.ConstructRegularSlabMesh(h, 0.165 /*length*/, 0.016 /*width*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(100.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("NoBlockScaledEZ20BCB100");
        HeartConfig::Instance()->SetOutputFilenamePrefix("NoBlockScaledEZ20BCB100");
		HeartConfig::Instance()->SetVisualizeWithVtk(true);
		
		HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

		// Original bidomain bath code
		std::set<unsigned> tissue_ids;
        static unsigned tissue_id=0;
        tissue_ids.insert(tissue_id);

        std::set<unsigned> bath_ids;
        static unsigned bath_id1=1;
        bath_ids.insert(bath_id1);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            double x = iter->CalculateCentroid()[0];
            double y = iter->CalculateCentroid()[1];
            if (sqrt((x-0.01125)*(x-0.01125) + (y-0.008)*(y-0.008)) < 0.004)  // cell1
            {
                if (sqrt((x-0.01125)*(x-0.01125) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.01875)*(x-0.01875) + (y-0.008)*(y-0.008)) < 0.004)  //  cell2
            {
                if (sqrt((x-0.01875)*(x-0.01875) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.02625)*(x-0.02625) + (y-0.008)*(y-0.008)) < 0.004)  // cell3
            {
                if (sqrt((x-0.02625)*(x-0.02625) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.03375)*(x-0.03375) + (y-0.008)*(y-0.008)) < 0.004)  // cell4
            {
                if (sqrt((x-0.03375)*(x-0.03375) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.04125)*(x-0.04125) + (y-0.008)*(y-0.008)) < 0.004)  // cell5
            {
                if (sqrt((x-0.04125)*(x-0.04125) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.04875)*(x-0.04875) + (y-0.008)*(y-0.008)) < 0.004)  // cell6
            {
                if (sqrt((x-0.04875)*(x-0.04875) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.05625)*(x-0.05625) + (y-0.008)*(y-0.008)) < 0.004)  // cell7
            {
                if (sqrt((x-0.05625)*(x-0.05625) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.06375)*(x-0.06375) + (y-0.008)*(y-0.008)) < 0.004)  // cell8
            {
                if (sqrt((x-0.06375)*(x-0.06375) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.07125)*(x-0.07125) + (y-0.008)*(y-0.008)) < 0.004)  // cell9
            {
                if (sqrt((x-0.07125)*(x-0.07125) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.07875)*(x-0.07875) + (y-0.008)*(y-0.008)) < 0.004)  // cell10
            {
                if (sqrt((x-0.07875)*(x-0.07875) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.08625)*(x-0.08625) + (y-0.008)*(y-0.008)) < 0.004)  // cell11
            {
                if (sqrt((x-0.08625)*(x-0.08625) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.09375)*(x-0.09375) + (y-0.008)*(y-0.008)) < 0.004)  // cell12
            {
                if (sqrt((x-0.09375)*(x-0.09375) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.10125)*(x-0.10125) + (y-0.008)*(y-0.008)) < 0.004)  // cell13
            {
                if (sqrt((x-0.10125)*(x-0.10125) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.10875)*(x-0.10875) + (y-0.008)*(y-0.008)) < 0.004)  // cell14
            {
                if (sqrt((x-0.10875)*(x-0.10875) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.11625)*(x-0.11625) + (y-0.008)*(y-0.008)) < 0.004)  // cell15
            {
                if (sqrt((x-0.11625)*(x-0.11625) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.12375)*(x-0.12375) + (y-0.008)*(y-0.008)) < 0.004)  // cell16
            {
                if (sqrt((x-0.12375)*(x-0.12375) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.13125)*(x-0.13125) + (y-0.008)*(y-0.008)) < 0.004)  // cell17
            {
                if (sqrt((x-0.13125)*(x-0.13125) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.13875)*(x-0.13875) + (y-0.008)*(y-0.008)) < 0.004)  // cell18
            {
                if (sqrt((x-0.13875)*(x-0.13875) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.14625)*(x-0.14625) + (y-0.008)*(y-0.008)) < 0.004)  // cell19
            {
                if (sqrt((x-0.14625)*(x-0.14625) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
			else if (sqrt((x-0.15375)*(x-0.15375) + (y-0.008)*(y-0.008)) < 0.004)  // cell20
            {
                if (sqrt((x-0.15375)*(x-0.15375) + (y-0.008)*(y-0.008)) < 0.0025)
				{	
				iter->SetAttribute(bath_id1);
				}
				else
				{	
				//IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
				}
			}
            else
            {
                //Outside circle on the bottom
                 iter->SetAttribute(bath_id1);
            }
        }

        mesh.SetMeshHasChangedSinceLoading();

        HeartConfig::Instance()->SetBathConductivity(7.0);  //bath_id1 tags will take the default value (actually 7.0 is the default)
        //std::map<unsigned, double> multiple_bath_conductivities;
        //multiple_bath_conductivities[bath_id2] = 6.5;  // mS/cm
		
        //HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);
		
        // For default conductivities and explicit cell model -1e4 is under threshold, -1.4e4 too high - crashes the cell model
        // For heterogeneous conductivities as given, -1e4 is under threshold
        double magnitude = -20.0e3; // uA/cm^2
        double start_time = 0.0;
        double duration = 1; //ms
		
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);
		
		
		CustomCellFactory cell_factory;
		
		/*Create the problem class*/
        BidomainProblem<2> bidomain_problem( &cell_factory, true);

        bidomain_problem.SetMesh(&mesh);
		bidomain_problem.SetWriteInfo();
		 
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        Vec solution = bidomain_problem.GetSolution(); // the Vs and phi_e's, as a PetSc vector
        ReplicatableVector solution_repl(solution);
		
		HeartEventHandler::Headings();
		HeartEventHandler::Headings();

    }
};