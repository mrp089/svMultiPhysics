// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOLUTION_STATES_H
#define SOLUTION_STATES_H

#include "Array.h"

/**
 * @brief Represents solution variables at a single time level
 *
 * Contains the three primary solution arrays used in time integration:
 * A (time derivative), D (integrated variable), and Y (variable)
 */
struct Solution {
  Array<double> A;  ///< Time derivative (acceleration in structural mechanics)
  Array<double> D;  ///< Integrated variable (displacement in structural mechanics)
  Array<double> Y;  ///< Variable (velocity in structural mechanics)

  // Semantic getters for improved readability
  Array<double>& get_acceleration() { return A; }
  const Array<double>& get_acceleration() const { return A; }

  Array<double>& get_velocity() { return Y; }
  const Array<double>& get_velocity() const { return Y; }

  Array<double>& get_displacement() { return D; }
  const Array<double>& get_displacement() const { return D; }
};

/**
 * @brief Holds solution state at old and current time levels
 *
 * Contains solution arrays at two time levels for time integration:
 * - old: Previous converged solution at time n
 * - current: Current solution being computed at time n+1
 */
struct SolutionStates {
  Solution old;      ///< Previous converged solution at time n (Ao, Do, Yo)
  Solution current;  ///< Current solution being computed at time n+1 (An, Dn, Yn)
};

#endif
