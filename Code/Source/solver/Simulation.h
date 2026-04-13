// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIMULATION_H 
#define SIMULATION_H 

#include "ComMod.h"
#include "SolutionStates.h"
#include "Parameters.h"
#include "SimulationLogger.h"
#include "LinearAlgebra.h"

#include <string>
#include <memory>

// Forward declarations
class Integrator;
class PartitionedFSI;

class Simulation {

  public:
    Simulation();
    ~Simulation();

    const mshType& get_msh(const std::string& name);

    CepMod& get_cep_mod() { return cep_mod; };
    ChnlMod& get_chnl_mod() { return chnl_mod; };
    ComMod& get_com_mod() { return com_mod; };
    Integrator& get_integrator();
    PartitionedFSI* get_partitioned_fsi();
    void initialize_partitioned_fsi(const std::string& xml_file_path);
};

#endif

