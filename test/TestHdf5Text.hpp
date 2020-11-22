/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTHDF5CONVERTERS_HPP_
#define TESTHDF5CONVERTERS_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "Hdf5ToTxtConverter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToXdmfConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"


#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

typedef Hdf5ToVtkConverter<3,3> VTK_3D;
typedef Hdf5ToMeshalyzerConverter<3,3> MESHA_3D;

/* HOW_TO_TAG Cardiac/Post-processing
 * Convert already generated simulation (HDF5) results to text, VTK or Meshalyzer format.
 */

class TestHdf5Converters : public CxxTest::TestSuite
{
private:
    // Copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        OutputFileHandler handler(dir);
        FileFinder file_finder(file, RelativeTo::ChasteSourceRoot);

        FileFinder copied_file = handler.CopyFileTo(file_finder);
        TS_ASSERT(copied_file.IsFile());
    }

public:

    /**
     * This tests the HDF5 to .txt converter using a 3D example
     * taken from a bidomain simulation.
     */
    void TestBidomainTxtConversion3D()
    {
        std::string working_directory = "TestHdf5ToTxtConverter_bidomain";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToTxtConverter_bidomain,
         * as that is where the reader reads from.
         */
        /*CopyToTestOutputDirectory("pde/test/data/BathResults.h5",
                                  working_directory);
								  */

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToTxtConverter<2,2> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                          "BathResults", &mesh);

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("BathResults_V_0.txt");
        files_to_compare.push_back("BathResults_V_1.txt");
        files_to_compare.push_back("BathResults_Phi_e_0.txt");
        files_to_compare.push_back("BathResults_Phi_e_1.txt");

        /*for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            std::cout << "Comparing generated and reference " << files_to_compare[i] << std::endl;
            FileFinder generated_file(working_directory +"/txt_output/" + files_to_compare[i], RelativeTo::ChasteTestOutput);
            FileFinder reference_file("pde/test/data/" + files_to_compare[i], RelativeTo::ChasteSourceRoot);
            NumericFileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }*/
    }

   
};

#endif /*TESTHDF5CONVERTERS_HPP_*/
