/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "gtest/gtest.h"
#include "Simulation.h"
#include "ComMod.h"
#include "Parameters.h"
#include "read_files.h"
#include "read_msh.h"
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace read_files_ns;
using namespace read_msh_ns;

class ReadFilesTest : public ::testing::Test {
protected:
    Simulation* simulation;
    ComMod com_mod;
    Parameters params;
    std::string test_dir;

    void SetUp() override {
        simulation = new Simulation();
        test_dir = std::filesystem::temp_directory_path().string() + "/sv_read_test";
        std::filesystem::create_directories(test_dir);
    }

    void TearDown() override {
        delete simulation;
        std::filesystem::remove_all(test_dir);
    }

    // Helper function to create a temporary XML file
    std::string createTempSolverXml(const std::string& content) {
        std::string file_path = test_dir + "/solver_test.xml";
        std::ofstream file(file_path);
        file << content;
        file.close();
        return file_path;
    }

    // Helper function to create a temporary data file
    std::string createTempDataFile(const std::string& content, const std::string& filename) {
        std::string file_path = test_dir + "/" + filename;
        std::ofstream file(file_path);
        file << content;
        file.close();
        return file_path;
    }
};

// Test face_match function
TEST_F(ReadFilesTest, TestFaceMatch) {
    // Set up test data
    faceType lFa, gFa;
    int nsd = 3;  // 3D 
    com_mod.nsd = nsd;
    
    // Initialize the face data with simple coordinates
    lFa.nNo = 4;  // Quad face
    gFa.nNo = 4;
    
    lFa.x.resize(nsd, lFa.nNo);
    gFa.x.resize(nsd, gFa.nNo);
    
    // Local face vertices
    lFa.x(0,0) = 0.0; lFa.x(1,0) = 0.0; lFa.x(2,0) = 0.0;
    lFa.x(0,1) = 1.0; lFa.x(1,1) = 0.0; lFa.x(2,1) = 0.0;
    lFa.x(0,2) = 1.0; lFa.x(1,2) = 1.0; lFa.x(2,2) = 0.0;
    lFa.x(0,3) = 0.0; lFa.x(1,3) = 1.0; lFa.x(2,3) = 0.0;
    
    // Global face vertices (same as local for direct match)
    gFa.x(0,0) = 0.0; gFa.x(1,0) = 0.0; gFa.x(2,0) = 0.0;
    gFa.x(0,1) = 1.0; gFa.x(1,1) = 0.0; gFa.x(2,1) = 0.0;
    gFa.x(0,2) = 1.0; gFa.x(1,2) = 1.0; gFa.x(2,2) = 0.0;
    gFa.x(0,3) = 0.0; gFa.x(1,3) = 1.0; gFa.x(2,3) = 0.0;
    
    Vector<int> ptr(lFa.nNo);
    
    // Call the face_match function
    face_match(com_mod, lFa, gFa, ptr);
    
    // Check that the vertices are matched correctly
    for (int i = 0; i < lFa.nNo; i++) {
        EXPECT_EQ(ptr(i), i) << "Vertex " << i << " not matched correctly";
    }
}

// Test read_bc function
TEST_F(ReadFilesTest, TestReadBC) {
    // Create a simple XML element for boundary condition
    std::string bc_xml = R"(
    <add_bcType>
        <bc name="INLET">
            <imposing_bc>STRONGLY</imposing_bc>
            <eq_type>FLUID</eq_type>
            <condition>DIRICHLET</condition>
            <value>1.0</value>
        </bc>
    </add_bcType>
    )";
    
    // Set up test data
    EquationParameters* eq_params = new EquationParameters();
    BoundaryConditionParameters* bc_params = new BoundaryConditionParameters();
    eqType lEq;
    bcType lBc;
    
    // Initialize equation type
    lEq.phys = EquationPhys::fluid;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(bc_xml.c_str());
    
    // Call read_bc with the XML element
    read_bc(simulation, eq_params, lEq, bc_params, lBc, doc.FirstChildElement("add_bcType")->FirstChildElement("bc"));
    
    // Check that the boundary condition was set correctly
    EXPECT_EQ(lBc.bType, static_cast<BCType>(BoundaryConditionType::Dirichlet));
    EXPECT_EQ(lBc.weakDir, false); // Strongly imposed
    EXPECT_DOUBLE_EQ(lBc.g, 1.0);  // Value
    
    // Clean up
    delete eq_params;
    delete bc_params;
}

// Test read_bct function
TEST_F(ReadFilesTest, TestReadBCT) {
    // Create a simple time data file for boundary condition
    std::string time_data = R"(
    # Time   Value1   Value2   Value3
    0.0      0.0      0.0      0.0
    0.1      0.1      0.2      0.3
    0.2      0.2      0.4      0.6
    0.3      0.3      0.6      0.9
    )";
    
    std::string file_path = createTempDataFile(time_data, "bct_test.dat");
    
    // Set up test data
    MBType lMB;
    faceType lFa;
    
    // Call read_bct with the test file
    read_bct(com_mod, lMB, lFa, file_path);
    
    // Check that the data was loaded correctly
    ASSERT_EQ(lMB.dof, 3);      // 3 columns of data beyond time
    ASSERT_EQ(lMB.nTP, 4);      // 4 time points
    
    // Check first and last time points
    EXPECT_DOUBLE_EQ(lMB.dt(0), 0.0);
    EXPECT_DOUBLE_EQ(lMB.dt(3), 0.3);
    
    // Check values at specific time points and DOFs
    EXPECT_DOUBLE_EQ(lMB.d(0,0), 0.0);
    EXPECT_DOUBLE_EQ(lMB.d(1,1), 0.2);
    EXPECT_DOUBLE_EQ(lMB.d(2,2), 0.6);
    EXPECT_DOUBLE_EQ(lMB.d(3,2), 0.9);
}

// Test read_bf function
TEST_F(ReadFilesTest, TestReadBF) {
    // Create a simple XML element for body force
    std::string bf_xml = R"(
    <add_bfType>
        <bf name="GRAVITY">
            <eq_type>FLUID</eq_type>
            <value>0.0 0.0 -9.81</value>
        </bf>
    </add_bfType>
    )";
    
    // Set up test data
    BodyForceParameters* bf_params = new BodyForceParameters();
    bfType lBf;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(bf_xml.c_str());
    
    // Call read_bf with the XML element
    read_bf(com_mod, bf_params, lBf, doc.FirstChildElement("add_bfType")->FirstChildElement("bf"));
    
    // Check that the body force was set correctly
    EXPECT_EQ(lBf.dof, 3);  // 3D vector
    EXPECT_DOUBLE_EQ(lBf.f(0), 0.0);
    EXPECT_DOUBLE_EQ(lBf.f(1), 0.0);
    EXPECT_DOUBLE_EQ(lBf.f(2), -9.81);
    
    // Clean up
    delete bf_params;
}

// Test read_domain function
TEST_F(ReadFilesTest, TestReadDomain) {
    // Create a simple XML element for domain
    std::string domain_xml = R"(
    <add_dmnType>
        <domain name="FLUID_DOMAIN">
            <eq_type>FLUID</eq_type>
            <material>BLOOD</material>
        </domain>
    </add_dmnType>
    )";
    
    // Set up test data
    EquationParameters* eq_params = new EquationParameters();
    EquationProps propList;
    eqType lEq;
    
    // Initialize equation type
    lEq.phys = EquationPhys::fluid;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(domain_xml.c_str());
    
    // Call read_domain with the XML element
    read_domain(simulation, eq_params, lEq, propList, lEq.phys, 
                doc.FirstChildElement("add_dmnType")->FirstChildElement("domain"));
    
    // Check that the domain was set correctly
    EXPECT_EQ(propList.dmnId, 1);  // First domain should have ID 1
    EXPECT_EQ(propList.dmnName, "FLUID_DOMAIN");
    
    // Clean up
    delete eq_params;
}

// Test read_eq function
TEST_F(ReadFilesTest, TestReadEq) {
    // Create a simple XML element for equation
    std::string eq_xml = R"(
    <add_eqType>
        <equation name="FLUID">
            <type>NAVIER-STOKES</type>
            <coupled>true</coupled>
            <min_iter>3</min_iter>
            <max_iter>50</max_iter>
            <tol_res>1.0e-5</tol_res>
        </equation>
    </add_eqType>
    )";
    
    // Set up test data
    EquationParameters* eq_params = new EquationParameters();
    eqType lEq;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(eq_xml.c_str());
    
    // Call read_eq with the XML element
    read_eq(simulation, eq_params, lEq, doc.FirstChildElement("add_eqType")->FirstChildElement("equation"));
    
    // Check that the equation parameters were set correctly
    EXPECT_EQ(lEq.phys, EquationPhys::fluid);
    EXPECT_EQ(lEq.coupled, true);
    EXPECT_EQ(lEq.minItr, 3);
    EXPECT_EQ(lEq.maxItr, 50);
    EXPECT_DOUBLE_EQ(lEq.tol, 1.0e-5);
    
    // Clean up
    delete eq_params;
}

// Test read_files function
TEST_F(ReadFilesTest, TestReadFiles) {
    // Create a minimal solver XML file
    std::string solver_xml = R"(
    <svFSIFile>
        <Simulation>
            <add_eqType>
                <equation name="FLUID">
                    <type>NAVIER-STOKES</type>
                    <coupled>true</coupled>
                    <min_iter>3</min_iter>
                    <max_iter>50</max_iter>
                    <tol_res>1.0e-5</tol_res>
                </equation>
            </add_eqType>
            <output>
                <format>VTK</format>
                <freq>10</freq>
            </output>
        </Simulation>
    </svFSIFile>
    )";
    
    std::string file_path = createTempSolverXml(solver_xml);
    
    // Call read_files with the test file
    read_files(simulation, file_path);
    
    // Check that the simulation parameters were set correctly
    EXPECT_EQ(simulation->com_mod.nEq, 1);  // One equation defined
    
    // Get the first equation
    auto& lEq = simulation->com_mod.eq[0];
    
    // Check equation parameters
    EXPECT_EQ(lEq.phys, EquationPhys::fluid);
    EXPECT_EQ(lEq.coupled, true);
    EXPECT_EQ(lEq.minItr, 3);
    EXPECT_EQ(lEq.maxItr, 50);
    EXPECT_DOUBLE_EQ(lEq.tol, 1.0e-5);
    
    // Check output parameters
    EXPECT_EQ(simulation->com_mod.o.fmt, OutputFormat::vtk);
    EXPECT_EQ(simulation->com_mod.o.freq, 10);
}

// Test read_mat_model function
TEST_F(ReadFilesTest, TestReadMatModel) {
    // Create a simple XML element for material model
    std::string mat_xml = R"(
    <add_mater>
        <material name="BLOOD">
            <type>NEWTONIAN</type>
            <density>1.06</density>
            <viscosity>0.04</viscosity>
        </material>
    </add_mater>
    )";
    
    // Set up test data
    EquationParameters* eq_params = new EquationParameters();
    DomainParameters* domain_params = new DomainParameters();
    dmnType lDmn;
    
    // Initialize domain type
    lDmn.phys = EquationPhys::fluid;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(mat_xml.c_str());
    
    // Call read_mat_model with the XML element
    read_mat_model(simulation, eq_params, domain_params, lDmn, 
                   doc.FirstChildElement("add_mater")->FirstChildElement("material"));
    
    // Check that the material model was set correctly
    EXPECT_EQ(lDmn.mID, static_cast<MaterialModelType>(Material_Fluid::newtonian));
    EXPECT_DOUBLE_EQ(lDmn.prop.at(PhysicalProperyType::FluidDensity), 1.06);
    EXPECT_DOUBLE_EQ(lDmn.prop.at(PhysicalProperyType::FluidViscosity), 0.04);
    
    // Clean up
    delete eq_params;
    delete domain_params;
}

// Test read_outputs function
TEST_F(ReadFilesTest, TestReadOutputs) {
    // Create a simple XML element for outputs
    std::string output_xml = R"(
    <output>
        <format>VTK</format>
        <freq>10</freq>
        <restart>100</restart>
        <folder>results</folder>
    </output>
    )";
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(output_xml.c_str());
    
    // Call read_outputs with the XML element
    read_outputs(simulation, doc.FirstChildElement("output"));
    
    // Check that the output parameters were set correctly
    EXPECT_EQ(simulation->com_mod.o.fmt, OutputFormat::vtk);
    EXPECT_EQ(simulation->com_mod.o.freq, 10);
    EXPECT_EQ(simulation->com_mod.o.rfrq, 100);
    EXPECT_EQ(simulation->com_mod.o.dir, "results");
}

// Test read_rmsh function
TEST_F(ReadFilesTest, TestReadRmsh) {
    // Create a simple XML element for remeshing
    std::string rmsh_xml = R"(
    <remeshing>
        <frequency>10</frequency>
        <max_edge_size>0.1</max_edge_size>
        <min_edge_size>0.01</min_edge_size>
    </remeshing>
    )";
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(rmsh_xml.c_str());
    
    // Call read_rmsh with the XML element
    read_rmsh(simulation, doc.FirstChildElement("remeshing"));
    
    // Check that the remeshing parameters were set correctly
    EXPECT_EQ(simulation->com_mod.rmsh.isRmsh, true);
    EXPECT_EQ(simulation->com_mod.rmsh.freq, 10);
    EXPECT_DOUBLE_EQ(simulation->com_mod.rmsh.maxEdgeSize, 0.1);
    EXPECT_DOUBLE_EQ(simulation->com_mod.rmsh.minEdgeSize, 0.01);
}

// Test read_temporal_values function
TEST_F(ReadFilesTest, TestReadTemporalValues) {
    // Create a simple time data file
    std::string time_data = R"(
    # Time   Value
    0.0      0.0
    0.1      0.1
    0.2      0.2
    0.3      0.3
    )";
    
    std::string file_path = createTempDataFile(time_data, "time_values.dat");
    
    // Set up test data
    Vector<double> t, v;
    
    // Call read_temporal_values with the test file
    read_temporal_values(file_path, t, v);
    
    // Check that the data was loaded correctly
    ASSERT_EQ(t.size(), 4);  // 4 time points
    ASSERT_EQ(v.size(), 4);  // 4 values
    
    // Check first and last time points and values
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(v(0), 0.0);
    EXPECT_DOUBLE_EQ(t(3), 0.3);
    EXPECT_DOUBLE_EQ(v(3), 0.3);
}

// Test read_spatial_values function
TEST_F(ReadFilesTest, TestReadSpatialValues) {
    // Create a simple spatial data file
    std::string spatial_data = R"(
    # X   Y   Z   Value
    0.0   0.0   0.0   0.0
    1.0   0.0   0.0   0.1
    1.0   1.0   0.0   0.2
    0.0   1.0   0.0   0.3
    )";
    
    std::string file_path = createTempDataFile(spatial_data, "spatial_values.dat");
    
    // Set up test data
    Array<double> x;
    Vector<double> v;
    int nsd = 3;  // 3D
    
    // Call read_spatial_values with the test file
    read_spatial_values(file_path, x, v, nsd);
    
    // Check that the data was loaded correctly
    ASSERT_EQ(x.nrows(), nsd);  // 3D coordinates
    ASSERT_EQ(x.ncols(), 4);    // 4 points
    ASSERT_EQ(v.size(), 4);     // 4 values
    
    // Check coordinates and values for specific points
    EXPECT_DOUBLE_EQ(x(0,0), 0.0);
    EXPECT_DOUBLE_EQ(x(1,0), 0.0);
    EXPECT_DOUBLE_EQ(x(2,0), 0.0);
    EXPECT_DOUBLE_EQ(v(0), 0.0);
    
    EXPECT_DOUBLE_EQ(x(0,2), 1.0);
    EXPECT_DOUBLE_EQ(x(1,2), 1.0);
    EXPECT_DOUBLE_EQ(x(2,2), 0.0);
    EXPECT_DOUBLE_EQ(v(2), 0.2);
}

// Test read_temp_spat_values function
TEST_F(ReadFilesTest, TestReadTempSpatValues) {
    // Create a simple temporal-spatial data file
    std::string temp_spat_data = R"(
    # Time=0.0
    # X   Y   Z   Value
    0.0   0.0   0.0   0.00
    1.0   0.0   0.0   0.10
    # Time=0.1
    # X   Y   Z   Value
    0.0   0.0   0.0   0.01
    1.0   0.0   0.0   0.11
    )";
    
    std::string file_path = createTempDataFile(temp_spat_data, "temp_spat_values.dat");
    
    // Set up test data
    Vector<double> t;
    Array3<double> x;
    Array<double> v;
    int nsd = 3;  // 3D
    
    // Call read_temp_spat_values with the test file
    read_temp_spat_values(file_path, t, x, v, nsd);
    
    // Check that the data was loaded correctly
    ASSERT_EQ(t.size(), 2);     // 2 time points
    ASSERT_EQ(x.dim1(), 2);     // 2 time points
    ASSERT_EQ(x.dim2(), nsd);   // 3D coordinates
    ASSERT_EQ(x.dim3(), 2);     // 2 spatial points per time
    ASSERT_EQ(v.nrows(), 2);    // 2 time points
    ASSERT_EQ(v.ncols(), 2);    // 2 values per time
    
    // Check time points
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(1), 0.1);
    
    // Check coordinates and values for specific points and times
    EXPECT_DOUBLE_EQ(x(0,0,0), 0.0);  // Time 0, Point 0, X
    EXPECT_DOUBLE_EQ(x(0,1,0), 0.0);  // Time 0, Point 0, Y
    EXPECT_DOUBLE_EQ(x(0,2,0), 0.0);  // Time 0, Point 0, Z
    EXPECT_DOUBLE_EQ(v(0,0), 0.00);   // Value at Time 0, Point 0
    
    EXPECT_DOUBLE_EQ(x(1,0,1), 1.0);  // Time 1, Point 1, X
    EXPECT_DOUBLE_EQ(x(1,1,1), 0.0);  // Time 1, Point 1, Y
    EXPECT_DOUBLE_EQ(x(1,2,1), 0.0);  // Time 1, Point 1, Z
    EXPECT_DOUBLE_EQ(v(1,1), 0.11);   // Value at Time 1, Point 1
}

// Test read_fluid_visc_model function
TEST_F(ReadFilesTest, TestReadFluidViscModel) {
    // Create a simple XML element for fluid viscosity model
    std::string visc_xml = R"(
    <viscosity_model>
        <type>CARREAU-YASUDA</type>
        <viscosity_inf>0.0035</viscosity_inf>
        <viscosity_0>0.16</viscosity_0>
        <lambda>8.2</lambda>
        <a>0.64</a>
        <n>0.2128</n>
    </viscosity_model>
    )";
    
    // Set up test data
    dmnType lDmn;
    
    // Initialize domain type
    lDmn.phys = EquationPhys::fluid;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(visc_xml.c_str());
    
    // Call read_fluid_visc_model with the XML element
    read_fluid_visc_model(simulation, lDmn, doc.FirstChildElement("viscosity_model"));
    
    // Check that the viscosity model was set correctly
    EXPECT_EQ(lDmn.mID, static_cast<MaterialModelType>(Material_Fluid::carreau_yasuda));
    EXPECT_DOUBLE_EQ(lDmn.viscosity.visc_inf, 0.0035);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.visc_0, 0.16);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.lambda, 8.2);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.a, 0.64);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.n, 0.2128);
}

// Test read_solid_visc_model function
TEST_F(ReadFilesTest, TestReadSolidViscModel) {
    // Create a simple XML element for solid viscosity model
    std::string visc_xml = R"(
    <viscosity_model>
        <type>NEWTONIAN</type>
        <viscoelastic_modulus>1000.0</viscoelastic_modulus>
        <viscosity>10.0</viscosity>
    </viscosity_model>
    )";
    
    // Set up test data
    dmnType lDmn;
    
    // Initialize domain type
    lDmn.phys = EquationPhys::struct_; // Solid/structural domain
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(visc_xml.c_str());
    
    // Call read_solid_visc_model with the XML element
    read_solid_visc_model(simulation, lDmn, doc.FirstChildElement("viscosity_model"));
    
    // Check that the viscosity model was set correctly
    EXPECT_EQ(lDmn.viscType, ConstitutiveViscoType::newtonian);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.beta, 10.0);
    EXPECT_DOUBLE_EQ(lDmn.viscosity.visco_mod, 1000.0);
}

// Test read_wall_props_ff function
TEST_F(ReadFilesTest, TestReadWallPropsFF) {
    // Create a simple wall properties data file
    std::string wall_data = R"(
    # X   Y   Z   Thickness   Young's Modulus   Poisson's Ratio
    0.0   0.0   0.0   0.1   1000000.0   0.3
    1.0   0.0   0.0   0.1   1000000.0   0.3
    1.0   1.0   0.0   0.1   1000000.0   0.3
    0.0   1.0   0.0   0.1   1000000.0   0.3
    )";
    
    std::string file_path = createTempDataFile(wall_data, "wall_props.dat");
    
    // Set up test data
    faceType lFa;
    int nsd = 3;  // 3D
    com_mod.nsd = nsd;
    
    // Initialize face data
    lFa.nNo = 4;
    lFa.x.resize(nsd, lFa.nNo);
    
    // Set face coordinates
    lFa.x(0,0) = 0.0; lFa.x(1,0) = 0.0; lFa.x(2,0) = 0.0;
    lFa.x(0,1) = 1.0; lFa.x(1,1) = 0.0; lFa.x(2,1) = 0.0;
    lFa.x(0,2) = 1.0; lFa.x(1,2) = 1.0; lFa.x(2,2) = 0.0;
    lFa.x(0,3) = 0.0; lFa.x(1,3) = 1.0; lFa.x(2,3) = 0.0;
    
    // Call read_wall_props_ff with the test file
    read_wall_props_ff(com_mod, lFa, file_path);
    
    // Check that the wall properties were loaded correctly
    ASSERT_TRUE(lFa.varWallProps);
    ASSERT_EQ(lFa.thickness.size(), 4);
    ASSERT_EQ(lFa.E.size(), 4);
    ASSERT_EQ(lFa.nu.size(), 4);
    
    // Check values for specific points
    EXPECT_DOUBLE_EQ(lFa.thickness(0), 0.1);
    EXPECT_DOUBLE_EQ(lFa.E(0), 1000000.0);
    EXPECT_DOUBLE_EQ(lFa.nu(0), 0.3);
    
    EXPECT_DOUBLE_EQ(lFa.thickness(3), 0.1);
    EXPECT_DOUBLE_EQ(lFa.E(3), 1000000.0);
    EXPECT_DOUBLE_EQ(lFa.nu(3), 0.3);
}

// Test set_cmm_bdry function
TEST_F(ReadFilesTest, TestSetCmmBdry) {
    // Create a simple XML element for CMM boundary
    std::string cmm_xml = R"(
    <cmm_bc>
        <face>INLET</face>
        <value>1.0</value>
    </cmm_bc>
    )";
    
    // Set up test data
    eqType lEq;
    ChnlMod& chnl = com_mod.chnl;
    mshType& lM = com_mod.msh[0];
    
    // Initialize equation type
    lEq.phys = EquationPhys::CMM;
    
    // Create a face in the mesh
    com_mod.msh.resize(1);
    lM.fa.resize(1);
    lM.fa[0].name = "INLET";
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(cmm_xml.c_str());
    
    // Call set_cmm_bdry with the XML element
    set_cmm_bdry(com_mod, lEq, doc.FirstChildElement("cmm_bc"));
    
    // Check that the CMM boundary condition was set correctly
    EXPECT_EQ(chnl.bcdof, 1);         // 1 DOF for the boundary condition
    EXPECT_DOUBLE_EQ(chnl.bcinp, 1.0); // Value of 1.0
}

// Test set_equation_properties function
TEST_F(ReadFilesTest, TestSetEquationProperties) {
    // Create minimal equation properties XML
    std::string eq_xml = R"(
    <Equations>
        <fluid>
            <viscosity>1.0</viscosity>
            <density>1.06</density>
        </fluid>
    </Equations>
    )";
    
    // Set up test data
    eqType lEq;
    
    // Initialize equation type
    lEq.phys = EquationPhys::fluid;
    
    // Parse XML to a structured element
    tinyxml2::XMLDocument doc;
    doc.Parse(eq_xml.c_str());
    
    // Call set_equation_properties with the XML element
    set_equation_properties(com_mod, lEq, doc.FirstChildElement("Equations"));
    
    // Check that the equation properties were set correctly
    EXPECT_DOUBLE_EQ(lEq.visc, 1.0);
    EXPECT_DOUBLE_EQ(lEq.dmn[0].prop.at(PhysicalProperyType::FluidDensity), 1.06);
}

// Test read_cplbc_initialization_file function
TEST_F(ReadFilesTest, TestReadCplbcInitializationFile) {
    // Create a simple coupling initialization data file
    std::string cpl_data = R"(
    # X   Y   Z   Value
    0.0   0.0   0.0   1.0
    1.0   0.0   0.0   1.1
    1.0   1.0   0.0   1.2
    0.0   1.0   0.0   1.3
    )";
    
    std::string file_path = createTempDataFile(cpl_data, "cpl_init.dat");
    
    // Set up test data
    cplFaceType lCpl;
    int nsd = 3;  // 3D
    com_mod.nsd = nsd;
    
    // Initialize coupling face data
    lCpl.fa.nNo = 4;
    lCpl.fa.x.resize(nsd, lCpl.fa.nNo);
    
    // Set face coordinates
    lCpl.fa.x(0,0) = 0.0; lCpl.fa.x(1,0) = 0.0; lCpl.fa.x(2,0) = 0.0;
    lCpl.fa.x(0,1) = 1.0; lCpl.fa.x(1,1) = 0.0; lCpl.fa.x(2,1) = 0.0;
    lCpl.fa.x(0,2) = 1.0; lCpl.fa.x(1,2) = 1.0; lCpl.fa.x(2,2) = 0.0;
    lCpl.fa.x(0,3) = 0.0; lCpl.fa.x(1,3) = 1.0; lCpl.fa.x(2,3) = 0.0;
    
    // Call read_cplbc_initialization_file with the test file
    read_cplbc_initialization_file(com_mod, lCpl, file_path);
    
    // Check that the coupling data was loaded correctly
    ASSERT_EQ(lCpl.y.size(), 4);
    
    // Check values for specific points
    EXPECT_DOUBLE_EQ(lCpl.y(0), 1.0);
    EXPECT_DOUBLE_EQ(lCpl.y(1), 1.1);
    EXPECT_DOUBLE_EQ(lCpl.y(2), 1.2);
    EXPECT_DOUBLE_EQ(lCpl.y(3), 1.3);
}

// Test read_fourier_coeff_values_file function
TEST_F(ReadFilesTest, TestReadFourierCoeffValuesFile) {
    // Create a simple Fourier coefficients data file
    std::string fourier_data = R"(
    # n   a_n   b_n
    0   1.0   0.0
    1   0.5   0.3
    2   0.2   0.1
    )";
    
    std::string file_path = createTempDataFile(fourier_data, "fourier_coeffs.dat");
    
    // Set up test data
    bcType lBc;
    faceType lFa;
    
    // Call read_fourier_coeff_values_file with the test file
    read_fourier_coeff_values_file(com_mod, lBc, lFa, file_path);
    
    // Check that the Fourier coefficients were loaded correctly
    ASSERT_EQ(lBc.fourier_a.size(), 3);
    ASSERT_EQ(lBc.fourier_b.size(), 3);
    
    // Check coefficients
    EXPECT_DOUBLE_EQ(lBc.fourier_a(0), 1.0);
    EXPECT_DOUBLE_EQ(lBc.fourier_b(0), 0.0);
    
    EXPECT_DOUBLE_EQ(lBc.fourier_a(1), 0.5);
    EXPECT_DOUBLE_EQ(lBc.fourier_b(1), 0.3);
    
    EXPECT_DOUBLE_EQ(lBc.fourier_a(2), 0.2);
    EXPECT_DOUBLE_EQ(lBc.fourier_b(2), 0.1);
}

// Test read_trac_bcff function
TEST_F(ReadFilesTest, TestReadTracBCFF) {
    // Create a simple traction boundary condition data file
    std::string trac_data = R"(
    # X   Y   Z   Value_X   Value_Y   Value_Z
    0.0   0.0   0.0   1.0   0.0   0.0
    1.0   0.0   0.0   1.0   0.0   0.0
    1.0   1.0   0.0   1.0   0.0   0.0
    0.0   1.0   0.0   1.0   0.0   0.0
    )";
    
    std::string file_path = createTempDataFile(trac_data, "trac_bc.dat");
    
    // Set up test data
    faceType lFa;
    int nsd = 3;  // 3D
    com_mod.nsd = nsd;
    
    // Initialize face data
    lFa.nNo = 4;
    lFa.x.resize(nsd, lFa.nNo);
    
    // Set face coordinates
    lFa.x(0,0) = 0.0; lFa.x(1,0) = 0.0; lFa.x(2,0) = 0.0;
    lFa.x(0,1) = 1.0; lFa.x(1,1) = 0.0; lFa.x(2,1) = 0.0;
    lFa.x(0,2) = 1.0; lFa.x(1,2) = 1.0; lFa.x(2,2) = 0.0;
    lFa.x(0,3) = 0.0; lFa.x(1,3) = 1.0; lFa.x(2,3) = 0.0;
    
    // Call read_trac_bcff with the test file
    MBType lMB;
    read_trac_bcff(com_mod, lMB, lFa, file_path);
    
    // Check that the traction data was loaded correctly
    ASSERT_EQ(lMB.nTP, 1);  // 1 time point (assumed steady state since no time column)
    ASSERT_EQ(lMB.dof, 3);  // 3 components of traction vector
    
    // Check traction vector components for all nodes
    for (int i = 0; i < lFa.nNo; i++) {
        EXPECT_DOUBLE_EQ(lMB.d(0,i*nsd+0), 1.0);  // X component
        EXPECT_DOUBLE_EQ(lMB.d(0,i*nsd+1), 0.0);  // Y component
        EXPECT_DOUBLE_EQ(lMB.d(0,i*nsd+2), 0.0);  // Z component
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}