/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Kazushi (Jam) Marukawa <jam@pobox.com>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <minizinc/flattener.hh>
#include <minizinc/solver.hh>

// #include <geas/solver/solver.h>

// #define QUBO_USE_MZN_DIRECTLY

namespace MiniZinc {

class QuboOptions : public SolverInstanceBase::Options {
public:
  bool allSolutions = false;
  int conflicts = 0;
  bool freeSearch = false;
  int nrSolutions = 1;
  int objProbeLimit = 0;
  bool statistics = false;
  bool verbose = false;
  std::chrono::milliseconds time = std::chrono::milliseconds(0);
};

class QuboVariable {
public:
  enum Type { BOOL_TYPE, FLOAT_TYPE, INT_TYPE };

protected:
  /// Type of the variable
  Type _t;
  /// the index in the iv/bv/fv array in the space, depending the type _t
  unsigned int _index;

public:
  explicit QuboVariable(unsigned int i) : _t(BOOL_TYPE), _index(i) {};

  QuboVariable(const QuboVariable& gv) : _t(gv._t) {
    _index = gv._index;
  }

  bool isBool() const { return _t == BOOL_TYPE; }
  bool isFloat() const { return _t == FLOAT_TYPE; }
  bool isInt() const { return _t == INT_TYPE; }

  unsigned int index() const { return _index; }
};

class QuboTypes {
public:
  typedef QuboVariable Variable;
  typedef MiniZinc::Statistics Statistics;
};

class QuboSolverInstance : public SolverInstanceImpl<QuboTypes> {
public:
  QuboSolverInstance(Env& env, std::ostream& log, SolverInstanceBase::Options* opt);
  ~QuboSolverInstance() override = default;
  void processFlatZinc() override;
#if 0
  geas::solver_data* solverData() const { return _solver.data; }
  geas::solver& solver() { return _solver; }
#endif

  Status solve() override;
  Status next() override { return SolverInstance::ERROR; }  // TODO: Implement
  void resetSolver() override;

  Expression* getSolutionValue(Id* id) override;
  void printStatistics() override;

  // MiniZinc to VA conversions
#if 0
  bool asBool(Expression* e) { return eval_bool(env().envi(), e); }
  vec<bool> asBool(ArrayLit* al);
  geas::patom_t asBoolVar(Expression* e);
  vec<geas::patom_t> asBoolVar(ArrayLit* al);
  vec<int> asInt(ArrayLit* al);
  int asInt(Expression* e) { return static_cast<int>(eval_int(env().envi(), e).toInt()); }
  geas::intvar asIntVar(Expression* e);
  vec<geas::intvar> asIntVar(ArrayLit* al);

  // TODO: create only when necessary or use VA internal
  geas::intvar zero;
#endif

protected:
#if 0
  geas::solver _solver;
#endif

  SolveI::SolveType _objType = SolveI::ST_MIN;
  std::unique_ptr<QuboTypes::Variable> _objVar;

  QuboTypes::Variable& resolveVar(Expression* e);
};

class QuboSolverFactory : public SolverFactory {
public:
  QuboSolverFactory();
  SolverInstanceBase::Options* createOptions() override;
  SolverInstanceBase* doCreateSI(Env& env, std::ostream& log,
                                 SolverInstanceBase::Options* opt) override;

  std::string getDescription(SolverInstanceBase::Options* opt) override {
    return "VA - Another Lazy Clause Generation Solver";
  };
  std::string getVersion(SolverInstanceBase::Options* opt) override { return "0.0.1"; }
#ifdef QUBO_USE_MZN_DIRECTLY
  std::string getId() override { return "org.minizinc.mzn-mzn"; }
#else
  std::string getId() override { return "org.minizinc.qubo"; }
#endif

  bool processOption(SolverInstanceBase::Options* opt, int& i, std::vector<std::string>& argv,
                     const std::string& workingDir = std::string()) override;
  void printHelp(std::ostream& os) override;
};

}  // namespace MiniZinc
