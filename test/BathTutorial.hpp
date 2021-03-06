#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"

#include "ToRORddynClmid.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include <cmath>

#include "UblasIncludes.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "CellProperties.hpp"
#include "Exception.hpp"
#include <sstream>
#include "HeartEventHandler.hpp"

#include "Hdf5DataReader.hpp"
#include <string>

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
		
		bool partial_output = false;
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
#ifndef _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#include "Hdf5DataReader.hpp"
#include <string>
class PropagationPropertiesCalculator
{
private:
    /** Reader to get the data from which we use to calculate properties. */
    Hdf5DataReader* mpDataReader;
    /** Name of the variable representing the membrane potential. */
    const std::string mVoltageName;
    /** Time values */
    std::vector<double> mTimes;
    /** Which node voltages have been cached for, if any */
    unsigned mCachedNodeGlobalIndex;
    /** The cached voltages vector */
    std::vector<double> mCachedVoltages;

protected:
    /**
     * @return the voltages vector for the given node and cache it, returning a reference
     * to the cached vector.  If subsequently called with the same index, will return
     * the cached vector without re-reading from file.
     *
     * Note: will only cache the last node index used.
     *
     * @param globalNodeIndex  the index of the node to cache voltages for
     */
    std::vector<double>& rGetCachedVoltages(unsigned globalNodeIndex);

public:
    /**
     * Constructor.
     *
     * @param pDataReader  Pointer to the data reader containing the simulation.
     * @param voltageName  Optionally the name of the variable representing the
     *     membrane potential.  Defaults to "V".
     */
    PropagationPropertiesCalculator(Hdf5DataReader* pDataReader,
                                    const std::string voltageName = "V");

    /** Destructor */
    //virtual ~PropagationPropertiesCalculator();

    /**
     * @return the maximum upstroke velocity at a single cell.
     * We calculate for the last upstroke found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    //double CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex);

     /**
     * @return the maximum upstroke velocity at a single cell.
     * We return all the max upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold (we use this for marking the end of an AP)
     */
    //std::vector<double> CalculateAllMaximumUpstrokeVelocities(unsigned globalNodeIndex, double threshold);

     /**
     * @return the times of upstroke at a single cell.
     * We return all the times of upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    //std::vector<double> CalculateUpstrokeTimes(unsigned globalNodeIndex, double threshold);

    /**
     * @return the conduction velocity between two cells, i.e. the time
     * taken for an AP to propagate from one to the other. It returns
     * the value of conduction velocity of the LAST action potential
     * that reached both nodes. Throws exceptions if an AP never reached
     * one of the nodes.
     *
     * @param globalNearNodeIndex  The cell to measure from.
     * @param globalFarNodeIndex  The cell to measure to.
     * @param euclideanDistance  The distance the AP travels between the cells, along the tissue.
     */
	globalNearNodeIndex=58;
	globalNearNodeIndex=62;
	euclideanDistance=0.0.4;
    double CalculateConductionVelocity(unsigned globalNearNodeIndex,
                                       unsigned globalFarNodeIndex,
                                       const double euclideanDistance);

     /**
     * @return all the conduction velocities between two cells, i.e. the time
     * taken for all APs to propagate from one to the other. It returns a vector
     * containing all the conduction velocities for each of the APs that
     * reached the two nodes (only the APs that reached both nodes).
     * Throws exceptions if an AP never reached one of the nodes.
     *
     * @param globalNearNodeIndex  The cell to measure from.
     * @param globalFarNodeIndex  The cell to measure to.
     * @param euclideanDistance  The distance the AP travels between the cells, along the tissue.
     */
     std::vector<double> CalculateAllConductionVelocities(unsigned globalNearNodeIndex,
                                                          unsigned globalFarNodeIndex,
                                                          const double euclideanDistance);
    /**
     * @return the action potential duration at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            unsigned globalNodeIndex);
    /**
     * @return the maximum transmembrane potential (maximum systolic
     * potential) at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculatePeakMembranePotential(unsigned globalNodeIndex);

     /**
     * @return all the action potentials duration at a single cell.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold  The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<double> CalculateAllActionPotentialDurations(const double percentage,
                                                             unsigned globalNodeIndex,
                                                             double threshold);

     /**
     * @return all the action potentials duration at cells [lowerNodeIndex, upperNodeIndex-1].
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param lowerNodeIndex  First cell at which to calculate.
     * @param upperNodeIndex  One past the last cell at which to calculate.
     * @param threshold  The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<std::vector<double> > CalculateAllActionPotentialDurationsForNodeRange(const double percentage,
                                                                                       unsigned lowerNodeIndex,
                                                                                       unsigned upperNodeIndex,
                                                                                       double threshold);

     /**
      * @return all the depolarisations that occur above threshold at a single cell.
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    std::vector<unsigned> CalculateAllAboveThresholdDepolarisations(unsigned globalNodeIndex,
                                                                    double threshold);

     /**
      * @return the depolarisations that occur above threshold at a single cell during the last recorded Ap
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    unsigned CalculateAboveThresholdDepolarisationsForLastAp(unsigned globalNodeIndex,
                                                             double threshold);

    /**
     * Provide a new pointer to an HDF5 data reader
     *
     * @param pDataReader  An HDF5 data reader to use (needed if the existing one is deleted and a new one opened)
     */
    void SetHdf5DataReader(Hdf5DataReader* pDataReader);
};
#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_