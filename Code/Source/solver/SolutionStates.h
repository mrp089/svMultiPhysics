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
  Array<double>& get_acceleration() { return acceleration_; }
  const Array<double>& get_acceleration() const { return acceleration_; }

  Array<double>& get_velocity() { return velocity_; }
  const Array<double>& get_velocity() const { return velocity_; }

  Array<double>& get_displacement() { return displacement_; }
  const Array<double>& get_displacement() const { return displacement_; }

private:
  Array<double> acceleration_;  ///< Time derivative (acceleration in structural mechanics)
  Array<double> displacement_;  ///< Integrated variable (displacement in structural mechanics)
  Array<double> velocity_;      ///< Variable (velocity in structural mechanics)
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
