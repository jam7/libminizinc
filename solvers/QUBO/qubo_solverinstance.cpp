/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Kazushi (Jam) Marukawa <jam@pobox.com>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/solvers/QUBO/qubo_solverinstance.hh>

namespace MiniZinc {
QuboSolverInstance::QuboSolverInstance(Env& env, std::ostream& log,
                                       SolverInstanceBase::Options* opt)
    : SolverInstanceImpl<QuboTypes>(env, log, opt), _model(env.model()) {
}

void QuboSolverInstance::processFlatZinc() {
  std::cerr << "hello QuboSolverInstance::processFlatZinc()\n";
}

SolverInstanceBase::Status MiniZinc::QuboSolverInstance::solve() {
  std::cerr << "hello QuboSolverInstance::solve()\n";
  SolverInstanceBase::Status status = SolverInstance::ERROR;
  auto _opt = static_cast<QuboOptions&>(*_options);
#if 0
  auto remaining_time = [_opt] {
    if (_opt.time == std::chrono::milliseconds(0)) {
      return 0.0;
    }
    using geas_time = std::chrono::duration<double>;
    static auto timeout = std::chrono::high_resolution_clock::now() + _opt.time;
    return geas_time(timeout - std::chrono::high_resolution_clock::now()).count();
  };
#endif
  // Set objective
  SolveI* si = _model->solveItem();
  if (si->e() != nullptr) {
    _objType = si->st();
    _objVar = std::unique_ptr<QuboTypes::Variable>(new QuboTypes::Variable(resolveVar(si->e())));
    // Support only MIN and MAX.
    if (_objType != SolveI::ST_MIN || _objType != SolveI::ST_MAX) {
      return status;
    }
  }
#if 0
  if (_objType == SolveI::ST_MIN) {
    // TODO: Add float objectives
    assert(_objVar->isInt());
    geas::intvar obj = _objVar->intVar();
    geas::solver::result res;
    while (true) {
      res = _solver.solve({remaining_time(), _opt.conflicts - _solver.data->stats.conflicts});
      geas::intvar::val_t obj_val;
      if (res != geas::solver::SAT) {
        break;
      }
      status = SolverInstance::SAT;
      if (_opt.allSolutions) {
        printSolution();
      }
      obj_val = _solver.get_model()[obj];

      int step = 1;
      while (_opt.objProbeLimit > 0) {
        geas::intvar::val_t assumed_obj;
        assumed_obj = obj_val - step;
        assumed_obj = obj.lb(_solver.data) > assumed_obj ? obj.lb(_solver.data) : assumed_obj;
        if (!_solver.assume(obj == assumed_obj)) {
          _solver.retract();
          break;
        }
        res = _solver.solve({remaining_time(), _opt.objProbeLimit});
        _solver.retract();
        if (res != geas::solver::SAT) {
          break;
        }
        step *= 2;
        if (_opt.allSolutions) {
          printSolution();
        }
        obj_val = _solver.get_model()[obj];
      }
      _solver.post(obj < obj_val);
    }
    if (status == SolverInstance::ERROR) {
      switch (res) {
        case geas::solver::UNSAT:
          status = SolverInstance::UNSAT;
          break;
        case geas::solver::UNKNOWN:
          status = SolverInstance::UNKNOWN;
          break;
        default:
          assert(false);
          status = SolverInstance::ERROR;
          break;
      }
    } else {
      if (res == geas::solver::UNSAT) {
        status = SolverInstance::OPT;
      }
      if (!_opt.allSolutions) {
        printSolution();
      }
    }
  }
  if (_opt.statistics) {
    printStatistics();
  }
#endif
  return status;
}

Expression* QuboSolverInstance::getSolutionValue(Id* id) {
#if 0
  id = id->decl()->id();
  if (id->type().isvar()) {
    QuboVariable& var = resolveVar(id->decl()->id());
    geas::model solution = _solver.get_model();
    switch (id->type().bt()) {
      case Type::BT_BOOL:
        assert(var.isBool());
        return constants().boollit(solution.value(var.boolVar()));
      case Type::BT_FLOAT:
        assert(var.isFloat());
        return FloatLit::a(solution[var.floatVar()]);
      case Type::BT_INT:
        assert(var.isInt());
        return IntLit::a(solution[var.intVar()]);
      default:
        return nullptr;
    }
  } else {
    return id->decl()->e();
  }
#endif
  return nullptr;
}

void QuboSolverInstance::resetSolver() { assert(false); }

QuboTypes::Variable& QuboSolverInstance::resolveVar(Expression* e) {
  if (auto* id = e->dynamicCast<Id>()) {
    return _variableMap.get(id->decl()->id());
  }
  if (auto* vd = e->dynamicCast<VarDecl>()) {
    return _variableMap.get(vd->id()->decl()->id());
  }
  if (auto* aa = e->dynamicCast<ArrayAccess>()) {
    auto* ad = aa->v()->cast<Id>()->decl();
    auto idx = aa->idx()[0]->cast<IntLit>()->v().toInt();
    auto* al = eval_array_lit(_env.envi(), ad->e());
    return _variableMap.get((*al)[idx]->cast<Id>());
  }
  std::stringstream ssm;
  ssm << "Expected Id, VarDecl or ArrayAccess instead of \"" << *e << "\"";
  throw InternalError(ssm.str());
}

#if 0
vec<bool> QuboSolverInstance::asBool(ArrayLit* al) {
  vec<bool> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = asBool((*al)[i]);
  }
  return vec;
}
#endif

#if 0
geas::patom_t QuboSolverInstance::asBoolVar(Expression* e) {
  if (e->type().isvar()) {
    QuboVariable& var = resolveVar(follow_id_to_decl(e));
    assert(var.isBool());
    return var.boolVar();
  }
  if (auto* bl = e->dynamicCast<BoolLit>()) {
    return bl->v() ? geas::at_True : geas::at_False;
  }
  std::stringstream ssm;
  ssm << "Expected bool or int literal instead of: " << *e;
  throw InternalError(ssm.str());
}

vec<geas::patom_t> QuboSolverInstance::asBoolVar(ArrayLit* al) {
  vec<geas::patom_t> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = this->asBoolVar((*al)[i]);
  }
  return vec;
}

vec<int> QuboSolverInstance::asInt(ArrayLit* al) {
  vec<int> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = this->asInt((*al)[i]);
  }
  return vec;
}

geas::intvar QuboSolverInstance::asIntVar(Expression* e) {
  if (e->type().isvar()) {
    QuboVariable& var = resolveVar(follow_id_to_decl(e));
    assert(var.isInt());
    return var.intVar();
  }
  IntVal i;
  if (auto* il = e->dynamicCast<IntLit>()) {
    i = il->v().toInt();
  } else if (auto* bl = e->dynamicCast<BoolLit>()) {
    i = static_cast<long long>(bl->v());
  } else {
    std::stringstream ssm;
    ssm << "Expected bool or int literal instead of: " << *e;
    throw InternalError(ssm.str());
  }
  if (i == 0) {
    return zero;
  }
  return _solver.new_intvar(static_cast<geas::intvar::val_t>(i.toInt()),
                            static_cast<geas::intvar::val_t>(i.toInt()));
}

vec<geas::intvar> QuboSolverInstance::asIntVar(ArrayLit* al) {
  vec<geas::intvar> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = this->asIntVar((*al)[i]);
  }
  return vec;
}
#endif

void QuboSolverInstance::printStatistics() {
#if 0
  auto& st = _solver.data->stats;
  auto& out = getSolns2Out()->getOutput();

  out << "%%%mzn-stat: failures=" << st.conflicts << std::endl;  // TODO: Statistic name
  out << "%%%mzn-stat: solveTime=" << st.time << std::endl;
  out << "%%%mzn-stat: solutions=" << st.solutions << std::endl;
  out << "%%%mzn-stat: restarts=" << st.restarts << std::endl;
  out << "%%%mzn-stat: nogoods=" << st.num_learnts << std::endl;             // TODO: Statistic name
  out << "%%%mzn-stat: learntLiterals=" << st.num_learnt_lits << std::endl;  // TODO: Statistic name
#endif
}

QuboSolverFactory::QuboSolverFactory() {
  std::cerr << "QuboSovlerFactory\n";
  SolverConfig sc("org.minizinc.qubo", getVersion(nullptr));
  sc.name("Qubo");
  sc.mznlibVersion(1);
  // Doesn't support fzn.
  sc.supportsFzn(false);
  // Support mzn to qubo conversion.
  sc.supportsMzn(true);
  sc.description(getDescription(nullptr));
  // FIXME: Update tags
  sc.tags({
      "api",
      "cp",
      "float",
      "int",
      "lcg",
  });
  sc.stdFlags({"-a", "-f", "-n", "-s", "-t"});
  sc.extraFlags({
      SolverConfig::ExtraFlag("--conflicts",
                              "Limit the maximum number of conflicts to be used during solving.",
                              SolverConfig::ExtraFlag::FlagType::T_INT, {}, "0"),
      SolverConfig::ExtraFlag(
          "--obj-probe",
          "Number of conflicts to use to probe for better solutions after a new solution is found.",
          SolverConfig::ExtraFlag::FlagType::T_INT, {}, "0"),
  });
  SolverConfigs::registerBuiltinSolver(sc);
};

SolverInstanceBase::Options* QuboSolverFactory::createOptions() { return new QuboOptions; }

SolverInstanceBase* QuboSolverFactory::doCreateSI(Env& env, std::ostream& log,
                                                  SolverInstanceBase::Options* opt) {
  return new QuboSolverInstance(env, log, opt);
}

bool QuboSolverFactory::processOption(SolverInstanceBase::Options* opt, int& i,
                                      std::vector<std::string>& argv,
                                      const std::string& workingDir) {
  std::cerr << "hello QuboSolverFactory::processOption()\n";

  auto& _opt = static_cast<QuboOptions&>(*opt);
  CLOParser cop(i, argv);
  int num = -1;

#if 0
  if (cop.getOption("-a --all-solutions")) {
    _opt.allSolutions = true;
  } else if (cop.getOption("--conflicts", &num)) {
    if (num >= 0) {
      _opt.conflicts = num;
    }
  } else if (cop.getOption("-f")) {
    _opt.freeSearch = true;
  } else if (cop.getOption("-n", &num)) {
    if (num >= 0) {
      _opt.nrSolutions = num;
    }
  } else if (cop.getOption("--obj-probe", &num)) {
    if (num >= 0) {
      _opt.objProbeLimit = num;
    }
  } else if (cop.getOption("--solver-statistics -s")) {
    _opt.statistics = true;
  } else if (cop.getOption("--solver-time-limit -t", &num)) {
    if (num >= 0) {
      _opt.time = std::chrono::milliseconds(num);
    }
  } else
#endif
  if (cop.getOption("--verbose-solving")) {
      _opt.verbose = true;
  } else {
    return false;
  }
  return true;
}

void QuboSolverFactory::printHelp(std::ostream& os) {
  os << "VA solver plugin options:" << std::endl
     << "  --conflicts <int>" << std::endl
     << "    Limit the maximum number of conflicts to be used during solving." << std::endl
     << "  --obj-probe <int>" << std::endl
     << "    Number of conflicts to use to probe for better solutions after a new solution is "
        "found."
     << std::endl
     << std::endl;
}
}  // namespace MiniZinc
