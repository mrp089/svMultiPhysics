// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions defined here are used to run a simulation from the command line.
//
// Usage:
//
//   svMultiPhysics XML_FILE_NAME
//
#include "Simulation.h"
#include "Integrator.h"
#include "PartitionedFSI.h"
    integrator.step();

    #ifdef debug_iterate_solution
    dmsg << ">>> End of Newton iteration" << std::endl;
    #endif

    // IB treatment: interpolate flow data on IB mesh from background
    // fluid mesh for explicit coupling, update old solution for implicit
    // coupling
    //
    /* [NOTE] Not implemented.
    if (ibFlag) {
      CALL IB_INTERPYU(Yn, Dn)
      if (ib.cpld == ibCpld_I) {
        ib.Auo = ib.Aun
        ib.Ubo = ib.Ubn
      }
    }
    */

    if (com_mod.risFlag) {
      ris::ris_meanq(com_mod, cm_mod, solutions);
      ris::ris_status(com_mod, cm_mod);
      if (cm.mas(cm_mod)) {
        std::cout << "Iteration: " << com_mod.cTS << std::endl;
        for (int iProj = 0; iProj < com_mod.ris.nbrRIS; iProj++) {
          std::cout << "Status for RIS projection: " << iProj << std::endl;
          std::cout << "            RIS iteration: " << com_mod.ris.nbrIter(iProj) << std::endl;
          std::cout << "       Is the valve close? " << com_mod.ris.clsFlg[iProj] << std::endl;
          std::cout << "            The status is: " << com_mod.ris.status[iProj] << std::endl;
        }
      }

      if (!std::all_of(com_mod.ris.status.begin(), com_mod.ris.status.end(), [](bool s) { return s; })) {
        if (std::any_of(com_mod.ris.nbrIter.begin(), com_mod.ris.nbrIter.end(), [](int iter) { return iter <= 1; })) {
          if (cm.mas(cm_mod)) {
            std::cout << "Valve status just changed. Do not update" << std::endl;
          }
        } else {
            ris::ris_updater(com_mod, cm_mod, solutions);
        }
        // goto label_11;
      }
    }

    // Saving the TXT files containing average and fluxes (or ECGs)
    #ifdef debug_iterate_solution
    dmsg << "Saving the TXT files containing average and fluxes ..." << std::endl;
    dmsg << "Saving the TXT files containing ECGs ..." << std::endl;
    #endif

    txt_ns::txt(simulation, false, solutions);

    // If remeshing is required then save current solution.
    //
    if (com_mod.rmsh.isReqd) {
      l1 = ((cTS % com_mod.rmsh.cpVar) == 0);
      if (l1) {
        #ifdef debug_iterate_solution
        dmsg << "Saving last solution for remeshing." << std::endl; 
        #endif
        com_mod.rmsh.rTS = cTS - 1;
        com_mod.rmsh.time = time - dt;
        for (int i = 0; i < com_mod.rmsh.iNorm.size(); i++) {
          com_mod.rmsh.iNorm(i) = com_mod.eq[i].iNorm;
        }

        com_mod.rmsh.A0 = Ao;
        com_mod.rmsh.Y0 = Yo;
        com_mod.rmsh.D0 = Do;
      }
    }

    // Look for a file containg a time step to stop the simulation.
    //
    // stopTrigName = "STOP_SIM"
    //
    auto& stopTrigName = com_mod.stopTrigName;
    bool l1 = false;
    int stopTS = 0;
    int count = -1;

    if (cm.mas(cm_mod)) {
      if (FILE *fp = fopen(stopTrigName.c_str(), "r")) {
        l1 = true;
        count = fscanf(fp, "%d", &stopTS);

        if (count == 0) {
          stopTS = cTS;
        }
        fclose(fp);

      } else {
        stopTS = nTS;
      }
    }

    #ifdef debug_iterate_solution
    dmsg << "cm.bcast(cm_mod, &stopTS)  ..." << std::endl; 
    #endif

    cm.bcast(cm_mod, &stopTS);

    l1 = (cTS >= stopTS);
    l2 = ((cTS % com_mod.stFileIncr) == 0);

    #ifdef debug_iterate_solution
    dmsg; 
    dmsg << "stFileIncr: " << com_mod.stFileIncr; 
    dmsg << "l1: " << l1; 
    dmsg << "l2: " << l2; 
    #endif

    // Saving the result to restart bin file
    if (l1 || l2) {
       output::write_restart(simulation, com_mod.timeP, solutions);
    }

    // Writing results into the disk with VTU format
    //
    #ifdef debug_iterate_solution
    dmsg; 
    dmsg << "saveVTK: " << com_mod.saveVTK; 
    #endif

    if (com_mod.saveVTK) {
      l2 = ((cTS % com_mod.saveIncr) == 0);
      l3 = (cTS >= com_mod.saveATS);
      #ifdef debug_iterate_solution
      dmsg << "l2: " << l2; 
      dmsg << "l3: " << l3; 
      #endif

      if (l2 && l3) {
        output::output_result(simulation, com_mod.timeP, 3, iEqOld);
        bool lAvg = false;
        vtk_xml::write_vtus(simulation, solutions, lAvg);
      } else {
        output::output_result(simulation, com_mod.timeP, 2, iEqOld);
      }

    } else {
      output::output_result(simulation, com_mod.timeP, 2, iEqOld);
    }

    // [NOTE] Not implemented.
    //
    if (com_mod.pstEq) {
      //CALL OUTDNORM()
    }

    if (com_mod.ibFlag) {
      //CALL IB_OUTCPUT()
    }

    // [HZ] Part related to RIS0D
    if (cEq == 0 && com_mod.ris0DFlag) {
      ris::ris0d_status(com_mod, cm_mod, solutions);
    }

    // [HZ] Part related to unfitted RIS
    // If the valve is active, look at the pressure difference 
    if (com_mod.urisFlag) {
      for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
        com_mod.uris[iUris].cnt++;
        if (com_mod.uris[iUris].clsFlg) {
          uris::uris_meanp(com_mod, cm_mod, iUris, solutions);
          // if (com_mod.uris[iUris].cnt == 1) {
          //   // GOTO 11 // The GOTO Statement in the Fortran code
          // }
        } else {
          uris::uris_meanv(com_mod, cm_mod, iUris, solutions);
        }
        if (cm.mas(cm_mod)) {
          std::cout << " URIS surface: " << com_mod.uris[iUris].name << ", count: " << com_mod.uris[iUris].cnt << std::endl;
        }
      }

      if (com_mod.mvMsh) {
        uris::uris_update_disp(com_mod, cm_mod, solutions);
      }

      if (cm.mas(cm_mod)) {
        if (l2 && l3) {
          uris::uris_write_vtus(com_mod);
        }
      }
    }
    // end RIS/URIS stuff 

    // Exiting outer loop if l1
    if (l1) {
      break;
    }

    // Solution is stored here before replacing it at next time step
    //
    Ao = An;
    Yo = Yn;

    if (com_mod.dFlag) {
      Do = Dn;
    }
    com_mod.cplBC.xo = com_mod.cplBC.xn;

  } // End of outer loop

  #ifdef debug_iterate_solution
  dmsg << "End of outer loop" << std::endl;
  #endif

  //#ifdef debug_iterate_solution
  //dmsg << "=======  Simulation Finished   ========== " << std::endl;
  //#endif
}


void run_simulation(Simulation* simulation)
{
  auto* partitioned_fsi = simulation->get_partitioned_fsi();
  if (partitioned_fsi) {
    partitioned_fsi->run();
  } else {
    iterate_solution(simulation);
  }
}


/// @brief Run a simulation from the command line using the name of a solver input 
/// XML file as an argument.
//
int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "[svMultiPhysics] ERROR: The svMultiPhysics program requires the solver input XML file name as an argument." << std::endl;
    exit(1);
  }

  // Process extra arguments for XML parameter substitution.
  for (int i = 2; i < argc; i++) {
    std::string str(argv[i]);
    int pos = str.find("=");
    auto name = str.substr(0,pos);
    auto value = str.substr(pos+1,str.size());
  }

  std::cout << std::scientific << std::setprecision(16);

  // Initialize MPI.
  //
  int mpi_rank, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

#ifdef ENABLE_ARRAY_INDEX_CHECKING
  if (mpi_rank == 0) {
    std::cout << "WARNING: Index checking is enabled" << std::endl;
  }
#endif

  // Create a Simulation object that stores all data structures for a simulation.
  //
  // The MPI prociess rank is set in the cmType::new_cm() method called
  // from the Simulation constructor. 
  //
  auto simulation = new Simulation();
  auto& cm = simulation->com_mod.cm;
  std::string file_name(argv[1]);

  #define n_debug_main
  #ifdef debug_main
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  // Iterate for restarting a simulation after remeshing. 
  //
  while (true) {

    // Read in the solver commands .xml file.
    //
    #ifdef debug_main
    dmsg << "Read files " << " ... ";
    #endif
    read_files(simulation, file_name);
    
    // Distribute data to processors.
    #ifdef debug_main
    dmsg << "Distribute data to processors " << " ... ";
    #endif
    distribute(simulation);

    // Initialize simulation data.
    //
    Vector<double> init_time(3);

    #ifdef debug_main
    dmsg << "Initialize " << " ... ";
    #endif
    initialize(simulation, init_time);

    // Create LinearAlgebra objects for each equation.
    //
    for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
      auto& eq = simulation->com_mod.eq[iEq];
      add_eq_linear_algebra(simulation->com_mod, eq);
    }

    // Initialize partitioned FSI coupling if configured
    simulation->initialize_partitioned_fsi(file_name);

    #ifdef debug_main
    for (int iM = 0; iM < simulation->com_mod.nMsh; iM++) {
      dmsg << "---------- iM " << iM;
      dmsg << "msh[iM].nNo: " << simulation->com_mod.msh[iM].nNo;
      dmsg << "msh[iM].gnNo: " << simulation->com_mod.msh[iM].gnNo;
      dmsg << "msh[iM].nEl: " << simulation->com_mod.msh[iM].nEl;
      dmsg << "msh[iM].gnEl: " << simulation->com_mod.msh[iM].gnEl;
    }
    #endif

    // Run the simulation.
    run_simulation(simulation);

    #ifdef debug_main
    dmsg << "resetSim: " << simulation->com_mod.resetSim;
    #endif

    // Remesh and continue the simulation.
    //
    if (simulation->com_mod.resetSim) {
      #ifdef debug_main
      dmsg << "Calling remesh_restart" << " ..."; 
      #endif
      remesh::remesh_restart(simulation);
      #ifdef debug_main
      dmsg << "Continue the simulation " << " ";
      #endif
    } else {
      break;
    }
  }

   for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
      auto& eq = simulation->com_mod.eq[iEq];
      finalize_linear_algebra(eq);
    }

  MPI_Finalize();
}
