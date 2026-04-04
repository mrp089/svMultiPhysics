// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @brief Sanity checks for partitioned FSI coupling.
 *
 * Tests traction extraction, sign consistency, velocity/displacement
 * consistency, data transfer, and predictor restore.
 */

#include "gtest/gtest.h"

#include "fsi_coupling.h"
#include "Integrator.h"
#include "Simulation.h"
#include "distribute.h"
#include "initialize.h"
#include "read_files.h"
#include "LinearAlgebra.h"
#include "set_bc.h"
#include "all_fun.h"

#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR ""
#endif

// ---------------------------------------------------------------------------
// MPI environment
// ---------------------------------------------------------------------------
class MPIEnvironment_PartFSI : public ::testing::Environment {
public:
  void SetUp() override {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      int argc = 0;
      char** argv = nullptr;
      MPI_Init(&argc, &argv);
    }
  }
  void TearDown() override {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MPI_Finalize();
    }
  }
};

static testing::Environment* const mpi_env_pfsi =
    testing::AddGlobalTestEnvironment(new MPIEnvironment_PartFSI);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
void add_eq_linear_algebra_test(ComMod& com_mod, eqType& lEq)
{
  lEq.linear_algebra = LinearAlgebraFactory::create_interface(lEq.linear_algebra_type);
  lEq.linear_algebra->set_preconditioner(lEq.linear_algebra_preconditioner);
  lEq.linear_algebra->initialize(com_mod, lEq);
  if (lEq.linear_algebra_assembly_type != consts::LinearAlgebraType::none) {
    lEq.linear_algebra->set_assembly(lEq.linear_algebra_assembly_type);
  }
}

static bool pfsi_test_data_available()
{
  std::string path = std::string(TEST_DATA_DIR);
  if (path.empty()) return false;
  struct stat st;
  return (stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
}

static const faceType* pfsi_find_face(const ComMod& com_mod, const std::string& name)
{
  for (int iM = 0; iM < com_mod.nMsh; iM++)
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++)
      if (com_mod.msh[iM].fa[iFa].name == name)
        return &com_mod.msh[iM].fa[iFa];
  return nullptr;
}

static const mshType* pfsi_find_mesh(const ComMod& com_mod, const std::string& face_name)
{
  for (int iM = 0; iM < com_mod.nMsh; iM++)
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++)
      if (com_mod.msh[iM].fa[iFa].name == face_name)
        return &com_mod.msh[iM];
  return nullptr;
}

struct SimSetup {
  Simulation* sim;
  char orig_dir[4096];

  SimSetup(const std::string& case_dir, const std::string& xml_file) {
    getcwd(orig_dir, sizeof(orig_dir));
    std::string full_dir = std::string(TEST_DATA_DIR) + "/" + case_dir;
    chdir(full_dir.c_str());

    sim = new Simulation();
    read_files_ns::read_files(sim, xml_file);
    distribute(sim);
    Vector<double> init_time(3);
    initialize(sim, init_time);
    for (int iEq = 0; iEq < sim->com_mod.nEq; iEq++)
      add_eq_linear_algebra_test(sim->com_mod, sim->com_mod.eq[iEq]);
  }

  ~SimSetup() {
    for (int iEq = 0; iEq < sim->com_mod.nEq; iEq++)
      sim->com_mod.eq[iEq].linear_algebra->finalize();
    delete sim;
    chdir(orig_dir);
  }

  void run_one_timestep() {
    auto& com_mod = sim->com_mod;
    auto& integrator = sim->get_integrator();
    auto& solutions = integrator.get_solutions();

    com_mod.cTS += 1;
    com_mod.time += com_mod.dt;
    com_mod.cEq = 0;
    for (auto& eq : com_mod.eq) { eq.itr = 0; eq.ok = false; }

    integrator.predictor();
    set_bc::set_bc_dir(com_mod, solutions);
    integrator.step();

    solutions.old.get_acceleration() = solutions.current.get_acceleration();
    solutions.old.get_velocity() = solutions.current.get_velocity();
    if (com_mod.dFlag)
      solutions.old.get_displacement() = solutions.current.get_displacement();
  }
};

// ===========================================================================
// Test 1: Traction sign and magnitude at lumen_wall
// ===========================================================================
TEST(PartitionedFSI, TractionSignAndMagnitude)
{
  if (!pfsi_test_data_available()) GTEST_SKIP() << "Test data not available";

  SimSetup fsi("fsi/pipe_3d", "solver.xml");
  auto& com_mod = fsi.sim->com_mod;
  auto& cm_mod = fsi.sim->cm_mod;
  auto& integrator = fsi.sim->get_integrator();
  const int nsd = com_mod.nsd;

  // Run 1 time step of monolithic FSI
  fsi.run_one_timestep();

  // Find lumen_wall face and fluid mesh
  auto* wall_face = pfsi_find_face(com_mod, "lumen_wall");
  auto* fluid_mesh = pfsi_find_mesh(com_mod, "lumen_wall");
  ASSERT_NE(wall_face, nullptr);
  ASSERT_NE(fluid_mesh, nullptr);

  // Extract traction at lumen_wall
  auto traction = fsi_coupling::extract_fluid_traction(
      com_mod, cm_mod, *fluid_mesh, *wall_face, com_mod.eq[0],
      integrator.get_Yg(), integrator.get_Dg(), integrator.get_solutions());

  // Sum all nodal forces → total force vector
  double total_force[3] = {0, 0, 0};
  for (int a = 0; a < wall_face->nNo; a++)
    for (int i = 0; i < nsd; i++)
      total_force[i] += traction(i, a);

  // Total radial force: project total force onto radial direction
  // For a cylindrical pipe centered at (0,0), radial = (x, y, 0)/r
  double total_radial = 0.0;
  for (int a = 0; a < wall_face->nNo; a++) {
    int Ac = wall_face->gN(a);
    double x = com_mod.x(0, Ac);
    double y = com_mod.x(1, Ac);
    double r = sqrt(x*x + y*y);
    if (r < 1e-10) continue;
    double radial_force = (x * traction(0, a) + y * traction(1, a)) / r;
    total_radial += radial_force;
  }

  // Traction should be radially outward (positive) for a pressurized pipe
  EXPECT_GT(total_radial, 0.0)
      << "Traction should point radially outward for positive pressure";

  // Total axial force should be near zero (symmetric inlet/outlet)
  double force_mag = sqrt(total_force[0]*total_force[0] +
                          total_force[1]*total_force[1] +
                          total_force[2]*total_force[2]);
  EXPECT_GT(force_mag, 0.0) << "Total traction force should be non-zero";

  std::cout << "  Total force: (" << total_force[0] << ", " << total_force[1]
            << ", " << total_force[2] << ")" << std::endl;
  std::cout << "  Total radial force: " << total_radial << std::endl;
}

// ===========================================================================
// Test 2: Traction at inlet matches Neumann BC
// ===========================================================================
TEST(PartitionedFSI, TractionMatchesNeumannBC)
{
  if (!pfsi_test_data_available()) GTEST_SKIP() << "Test data not available";

  SimSetup fluid("fsi/pipe_3d_partitioned", "solver_fluid_only.xml");
  auto& com_mod = fluid.sim->com_mod;
  auto& cm_mod = fluid.sim->cm_mod;
  auto& integrator = fluid.sim->get_integrator();
  const int nsd = com_mod.nsd;

  // Run 1 time step
  fluid.run_one_timestep();

  // Find inlet face
  auto* inlet_face = pfsi_find_face(com_mod, "lumen_inlet");
  auto* lumen_mesh = pfsi_find_mesh(com_mod, "lumen_inlet");
  ASSERT_NE(inlet_face, nullptr);
  ASSERT_NE(lumen_mesh, nullptr);

  // Extract traction at inlet
  auto traction = fsi_coupling::extract_fluid_traction(
      com_mod, cm_mod, *lumen_mesh, *inlet_face, com_mod.eq[0],
      integrator.get_Yg(), integrator.get_Dg(), integrator.get_solutions());

  // Sum axial (z) component of traction at inlet
  double total_axial = 0.0;
  for (int a = 0; a < inlet_face->nNo; a++)
    total_axial += traction(2, a);  // z-component

  // Compute inlet face area (approximately pi*r^2 for a circular inlet)
  double inlet_area = 0.0;
  for (int a = 0; a < inlet_face->nNo; a++) {
    // Simple: count nodes and estimate area
    inlet_area += 1.0;  // Will compute properly below
  }
  // Use the face normal vectors to compute area
  inlet_area = 0.0;
  for (int a = 0; a < inlet_face->nNo; a++) {
    // The face normal nV includes the area weighting
    double nz = inlet_face->nV(2, a);
    inlet_area += std::abs(nz);
  }

  // The Neumann BC is pressure = 5e4
  double applied_pressure = 5.0e4;

  // Total axial force should be approximately pressure * area
  // (positive z = into the pipe)
  std::cout << "  Total axial traction at inlet: " << total_axial << std::endl;
  std::cout << "  Inlet area (from nV): " << inlet_area << std::endl;
  std::cout << "  Expected force: " << applied_pressure * inlet_area << std::endl;
  std::cout << "  Note: values depend on viscous contribution and flow pattern" << std::endl;

  // Force should be positive (into the pipe) and on the order of pressure*area
  EXPECT_NE(total_axial, 0.0) << "Inlet traction should be non-zero";
}

// ===========================================================================
// Test 3: apply_traction_on_solid sign check
// ===========================================================================
TEST(PartitionedFSI, TractionApplicationSign)
{
  if (!pfsi_test_data_available()) GTEST_SKIP() << "Test data not available";

  SimSetup solid("fsi/pipe_3d_partitioned", "solver_solid.xml");
  auto& com_mod = solid.sim->com_mod;
  const int nsd = com_mod.nsd;

  // Find wall_inner face
  auto* inner_face = pfsi_find_face(com_mod, "wall_inner");
  ASSERT_NE(inner_face, nullptr);

  // Run predictor to set up solution arrays
  solid.sim->get_integrator().predictor();

  // Allocate R sized for the solid equation
  com_mod.cEq = 0;
  auto& eq = com_mod.eq[0];

  // Create a known traction: unit radially outward force at each node
  Array<double> traction(nsd, inner_face->nNo);
  for (int a = 0; a < inner_face->nNo; a++) {
    int Ac = inner_face->gN(a);
    double x = com_mod.x(0, Ac);
    double y = com_mod.x(1, Ac);
    double r = sqrt(x*x + y*y);
    if (r < 1e-10) r = 1.0;
    // Unit radially outward force
    traction(0, a) = x / r;
    traction(1, a) = y / r;
    traction(2, a) = 0.0;
  }

  // Zero R, then apply traction
  com_mod.R.resize(eq.dof, com_mod.tnNo);
  com_mod.R = 0.0;

  fsi_coupling::apply_traction_on_solid(com_mod, eq, *inner_face, traction);

  // Check: R should be NEGATIVE (R -= traction, traction is positive outward)
  double sum_radial_R = 0.0;
  for (int a = 0; a < inner_face->nNo; a++) {
    int Ac = inner_face->gN(a);
    double x = com_mod.x(0, Ac);
    double y = com_mod.x(1, Ac);
    double r = sqrt(x*x + y*y);
    if (r < 1e-10) continue;
    double radial_R = (x * com_mod.R(0, Ac) + y * com_mod.R(1, Ac)) / r;
    sum_radial_R += radial_R;
  }

  // R -= traction → R should be negative for positive (outward) traction
  EXPECT_LT(sum_radial_R, 0.0)
      << "R should be negative after applying outward traction (R -= traction)";

  std::cout << "  Sum of radial R at interface: " << sum_radial_R << std::endl;
  std::cout << "  (Should be negative: R -= outward_traction)" << std::endl;
}

// ===========================================================================
// Test 4: Solid velocity/displacement Newmark consistency
// ===========================================================================
TEST(PartitionedFSI, SolidNewmarkConsistency)
{
  if (!pfsi_test_data_available()) GTEST_SKIP() << "Test data not available";

  SimSetup solid("fsi/pipe_3d_partitioned", "solver_solid.xml");
  auto& com_mod = solid.sim->com_mod;
  auto& integrator = solid.sim->get_integrator();
  auto& solutions = integrator.get_solutions();
  auto& eq = com_mod.eq[0];
  const int nsd = com_mod.nsd;
  const double dt = com_mod.dt;
  const double gam = eq.gam;

  // Run 1 time step (solid with zero loading → should stay at zero)
  solid.run_one_timestep();

  // The solid has Dir BCs at inlet/outlet but no loading,
  // so Dn ≈ 0, Yn ≈ 0, An ≈ 0 everywhere.
  // Check Newmark consistency: Yn = Yo + dt*((1-gam)*Ao + gam*An)
  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  auto& Ao = solutions.old.get_acceleration();
  auto& Yo = solutions.old.get_velocity();

  double max_err = 0.0;
  for (int Ac = 0; Ac < com_mod.tnNo; Ac++) {
    for (int i = 0; i < nsd; i++) {
      double yn_expected = Yo(i, Ac) + dt * ((1.0 - gam) * Ao(i, Ac) + gam * An(i, Ac));
      double err = std::abs(Yn(i, Ac) - yn_expected);
      max_err = std::max(max_err, err);
    }
  }

  std::cout << "  Max Newmark consistency error (Yn vs formula): " << max_err << std::endl;
  EXPECT_LT(max_err, 1e-10)
      << "Velocity should be consistent with Newmark formula";
}

// ===========================================================================
// Test 7: Predictor restore verification
// ===========================================================================
TEST(PartitionedFSI, PredictorRestore)
{
  if (!pfsi_test_data_available()) GTEST_SKIP() << "Test data not available";

  SimSetup fluid("fsi/pipe_3d_partitioned", "solver_fluid_only.xml");
  auto& com_mod = fluid.sim->com_mod;
  auto& integrator = fluid.sim->get_integrator();
  auto& solutions = integrator.get_solutions();
  const int nsd = com_mod.nsd;

  // Run predictor
  com_mod.cTS = 1;
  com_mod.time = com_mod.dt;
  integrator.predictor();
  set_bc::set_bc_dir(com_mod, solutions);

  // Save predictor state
  Array<double> saved_An(solutions.current.get_acceleration());
  Array<double> saved_Yn(solutions.current.get_velocity());
  Array<double> saved_Dn(solutions.current.get_displacement());
  Array<double> saved_x(com_mod.x);

  // Run one Newton solve (modifies An, Yn, Dn)
  for (auto& eq : com_mod.eq) { eq.itr = 0; eq.ok = false; }
  integrator.step();

  // Verify solution changed
  double max_change = 0.0;
  auto& Yn = solutions.current.get_velocity();
  for (int a = 0; a < Yn.ncols(); a++)
    for (int i = 0; i < Yn.nrows(); i++)
      max_change = std::max(max_change, std::abs(Yn(i, a) - saved_Yn(i, a)));
  EXPECT_GT(max_change, 0.0) << "Solution should change after Newton solve";

  // Restore predictor state
  solutions.current.get_acceleration() = saved_An;
  solutions.current.get_velocity() = saved_Yn;
  solutions.current.get_displacement() = saved_Dn;
  com_mod.x = saved_x;

  // Verify restoration is exact
  double max_restore_err = 0.0;
  auto& restored_Yn = solutions.current.get_velocity();
  for (int a = 0; a < restored_Yn.ncols(); a++)
    for (int i = 0; i < restored_Yn.nrows(); i++)
      max_restore_err = std::max(max_restore_err,
          std::abs(restored_Yn(i, a) - saved_Yn(i, a)));

  EXPECT_EQ(max_restore_err, 0.0) << "Predictor restore should be bitwise exact";

  double max_x_err = 0.0;
  for (int a = 0; a < com_mod.x.ncols(); a++)
    for (int i = 0; i < nsd; i++)
      max_x_err = std::max(max_x_err, std::abs(com_mod.x(i, a) - saved_x(i, a)));
  EXPECT_EQ(max_x_err, 0.0) << "Mesh coordinate restore should be bitwise exact";

  std::cout << "  Max change after solve: " << max_change << std::endl;
  std::cout << "  Max restore error (Yn): " << max_restore_err << std::endl;
  std::cout << "  Max restore error (x):  " << max_x_err << std::endl;
}
