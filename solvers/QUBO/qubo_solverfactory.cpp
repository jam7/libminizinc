/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Kazushi (Jam) Marukawa <jam@pobox.com>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/solvers/QUBO/qubo_solverfactory.hh>
#include <minizinc/solvers/QUBO/qubo_solverinstance.hh>

namespace MiniZinc {
namespace {
void get_wrapper() { static QuboSolverFactory _qubo_solverfactory; }
}  // namespace
QuboSolverFactoryInitialiser::QuboSolverFactoryInitialiser() { get_wrapper(); }
}  // namespace MiniZinc
