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
    : SolverInstanceImpl<QuboTypes>(env, log, opt) {
}

// Normalize a given model
void QuboSolverInstance::normalizeModel() {
}

// Flatten a given model
void QuboSolverInstance::flattenModel() {
  auto _opt = static_cast<QuboOptions&>(*_options);
  FlatteningOptions fopts = FlatteningOptions();
  fopts.verbose = _opt.verbose;
  // TODO: get other options from MznSolver.
  //   _flt.setFlagStatistics(flagCompilerStatistics);
  //   _flt.setFlagTimelimit(flagOverallTimeLimit);
  try {
    flatten(env(), fopts);
  } catch (LocationException& e) {
    if (_opt.verbose) {
      _log << std::endl;
    }
    std::ostringstream errstream;
    errstream << e.what() << ": " << std::endl;
    env().dumpErrorStack(errstream);
    errstream << "  " << e.msg() << std::endl;
    throw Error(errstream.str());
  }

#if 0
  {
    Printer p(std::cerr);
    std::cerr << "Flat--------\n";
    p.print(env().flat());
    std::cerr << "SolveItem--------\n";
    p.print(env().flat()->solveItem());
  }
#endif

  // Doesn't perform MIP domain.

  // Optimize constraint chains.
  optimize(env(), true);

  // Warning checks.
  for (const auto& i : env().warnings()) {
    _log << "\n  WARNING: " << i;
  }
#if 0
  if (_compflags.werror && !env()->warnings().empty()) {
    throw Error("errors encountered");
  }
#endif
  env().clearWarnings();

  // Convert to old fzn.
  oldflatzinc(env());
}

// Process variables.
//
// Parse vardecls and register them through insertVar.
void QuboSolverInstance::processVariables() {
  auto _opt = static_cast<QuboOptions&>(*_options);

  SolveI* solveItem = getEnv()->flat()->solveItem();
  VarDecl* objVd = nullptr;

  if (solveItem->st() != SolveI::SolveType::ST_SAT) {
    if (Id* id = solveItem->e()->dynamicCast<Id>()) {
      objVd = id->decl();
    } else {
      std::cerr << "Objective must be Id: " << solveItem->e() << std::endl;
      throw InternalError("Objective must be Id");
    }
  }

  for (VarDeclIterator it = getEnv()->flat()->vardecls().begin();
       it != getEnv()->flat()->vardecls().end(); ++it) {
    if (it->removed()) {
      continue;
    }
    VarDecl* vd = it->e();
    if (!vd->ann().isEmpty()) {
      if (vd->ann().containsCall(constants().ann.output_array) ||
          vd->ann().contains(constants().ann.output_var)) {
        _varsWithOutput.push_back(vd);
        //         std::cerr << (*vd);
        //         if ( vd->e() )
        //           cerr << " = " << (*vd->e());
        //         cerr << endl;
      }
    }
    if (vd->type().dim() == 0 && it->e()->type().isvar() && !it->removed()) {
      MiniZinc::TypeInst* ti = it->e()->ti();
      QuboTypes::Variable::Type vType = QuboTypes::Variable::Type::FLOAT_TYPE;
      if (ti->type().isvarint() || ti->type().isint()) {
        vType = QuboTypes::Variable::Type::INT_TYPE;
      } else if (ti->type().isvarbool() || ti->type().isbool()) {
        vType = QuboTypes::Variable::Type::BOOL_TYPE;
      } else if (!(ti->type().isvarfloat() || ti->type().isfloat())) {
        std::stringstream ssm;
        ssm << "This type of var is not handled by MIP: " << *it << std::endl;
        ssm << "  VarDecl flags (ti, bt, st, ot): " << ti->type().ti() << ti->type().bt()
            << ti->type().st() << ti->type().ot() << ", dim == " << ti->type().dim()
            << "\nRemove the variable or add a constraint so it is redefined." << std::endl;
        throw InternalError(ssm.str());
      }
      double lb = 0.0;
      double ub = 1.0;  // for bool
      if (ti->domain() != nullptr) {
        if (QuboTypes::Variable::Type::FLOAT_TYPE == vType) {
          FloatBounds fb = compute_float_bounds(getEnv()->envi(), it->e()->id());
          if (fb.valid) {
            lb = fb.l.toDouble();
            ub = fb.u.toDouble();
          } else {
            lb = 1.0;
            ub = 0.0;
          }
        } else if (QuboTypes::Variable::Type::INT_TYPE == vType) {
          IntBounds ib = compute_int_bounds(getEnv()->envi(), it->e()->id());
          if (ib.valid) {  // Normally should be
            lb = static_cast<double>(ib.l.toInt());
            ub = static_cast<double>(ib.u.toInt());
          } else {
            lb = 1;
            ub = 0;
          }
        }
      } else if (QuboTypes::Variable::Type::BOOL_TYPE != vType) {
        lb = -INFINITY;  // if just 1 bound inf, using MZN's default?  TODO
        ub = -lb;
      }

      //       IntSetVal* dom = eval_intset(env,vdi->e()->ti()->domain());
      //       if (dom->size() > 1)
      //         throw runtime_error("MIPSolverinstance: domains with holes ! supported, use
      //         --MIPdomains");

      Id* id = it->e()->id();
      MZN_ASSERT_HARD(id == id->decl()->id());  // Assume all unified
      MZN_ASSERT_HARD(it->e() == id->decl());   // Assume all unified
      double obj = vd == objVd ? 1.0 : 0.0;
      auto* decl00 = follow_id_to_decl(it->e());
      MZN_ASSERT_HARD(decl00->isa<VarDecl>());
      {
        auto* vd00 = decl00->dynamicCast<VarDecl>();
        if (nullptr != vd00->e()) {
          // Should be a const
          auto dRHS = exprToConst(vd00->e());
          lb = std::max(lb, dRHS);
          ub = std::min(ub, dRHS);
        }
#if 0
        if (it->e() != vd00) {                             // A different vardecl
          res = exprToVar(vd00->id());                     // Assume FZN is sorted.
          MZN_ASSERT_HARD(!getMIPWrapper()->fPhase1Over);  // Still can change colUB, colObj
          /// Tighten the ini-expr's bounds
          lb = getMIPWrapper()->colLB.at(res) = std::max(getMIPWrapper()->colLB.at(res), lb);
          ub = getMIPWrapper()->colUB.at(res) = std::min(getMIPWrapper()->colUB.at(res), ub);
          if (0.0 != obj) {
            getMIPWrapper()->colObj.at(res) = obj;
          }
        } else {
          res = getMIPWrapper()->addVar(obj, lb, ub, vType, id->str().c_str());
        }
#endif
      }
      /// Test infeasibility
      if (lb > ub) {
        _status = SolverInstance::UNSAT;
        if (_opt.verbose) {
          std::cerr << "  VarDecl '" << *(it->e()) << "' seems infeasible: computed bounds [" << lb
                    << ", " << ub << ']' << std::endl;
        }
      }
      if (0.0 != obj) {
#if 0
        dObjVarLB = lb;
        dObjVarUB = ub;
        getMIPWrapper()->output.nObjVarIndex = res;
        if (_opt.verbose) {
          std::cerr << "  MIP: objective variable index (0-based): " << res << std::endl;
        }
#endif
      }
      insertVar(id, SolverVariable(id));
    }
  }
#if 0
  if (_opt.verbose && (!_mipWrapper->sLitValues.empty())) {
    std::cerr << "  MIPSolverinstance: during Phase 1,  " << _mipWrapper->nLitVars
              << " literals with " << _mipWrapper->sLitValues.size() << " values used."
              << std::endl;
  }
#endif
#if 0
  if (!getMIPWrapper()->fPhase1Over) {
    getMIPWrapper()->addPhase1Vars();
  }
#endif
}

// Process contaraints.
//
// Parse constraints...
void QuboSolverInstance::processConstraints() {
  auto _opt = static_cast<QuboOptions&>(*_options);

  // Constraints
#if 0
  for (ConstraintIterator it = env().flat()->constraints().begin(); it != env().flat()->constraints().end();
       ++it) {
    if (!it->removed()) {
      if (auto* c = it->e()->dynamicCast<Call>()) {
        auto* eVD = get_annotation(c->ann(), constants().ann.defines_var);
        _log << "constraint: " << *c << ", " << *eVD << "\n";
        // _constraintRegistry.post(c);
        Call* pC = eVD->dynamicCast<Call>();
        _log << "call: " << *pC << "\n";
        _log << "call arg0: " << *pC->arg(0) << "\n";
        Id* id = pC->arg(0)->dynamicCast<Id>();
        VarDecl* vd = id->decl();
        _log << "call arg0 vd: " << *vd << "\n";
      }
    }
  }
#endif
}

void QuboSolverInstance::processFlatZinc() {
  auto _opt = static_cast<QuboOptions&>(*_options);

#if 0
  {
    Printer p(std::cerr);
    std::cerr << "Model--------\n";
    p.print(env().model());
    std::cerr << "SolveItem--------\n";
    p.print(env().model()->solveItem());
  }
#endif

  // Process normalization here.
  normalizeModel();
  // Process flatten here.
  flattenModel();
  // Process variables.
  processVariables();
  // Process constraints.
  processConstraints();

#if 0
  {
    Printer p(std::cerr);
    std::cerr << "Flat--------\n";
    p.print(env().flat());
    std::cerr << "SolveItem--------\n";
    p.print(env().flat()->solveItem());
  }
#endif
}

// Calculate QUBO matrix.
bool MiniZinc::QuboSolverInstance::calcQubo(Expression* e) {

#if 1
  {
    Printer p(std::cerr);
    std::cerr << "solveItem\n";
    p.print(e);
  }
#endif
  VarDecl* objVd;
  if (Id* id = e->dynamicCast<Id>()) {
    objVd = id->decl();
  } else {
    std::cerr << "Objective must be Id: " << e << std::endl;
    throw InternalError("Objective must be Id");
  }
  _log << *objVd << "\n";
#if 1
  dumpVar();
#endif
  if (_objType == SolveI::ST_MIN) {
  } else {
    assert(_objType == SolveI::ST_MAX);
  }
#if 1
  {
    Printer p(std::cerr);
    std::cerr << "Flat--------\n";
    p.print(env().flat());
    std::cerr << "SolveItem--------\n";
    p.print(env().flat()->solveItem());
  }
#endif
  return true;
}

SolverInstanceBase::Status MiniZinc::QuboSolverInstance::solve() {
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
  SolveI* si = env().flat()->solveItem();
  if (si->e() != nullptr) {
    _objType = si->st();
    _objVar = std::unique_ptr<QuboTypes::Variable>(new QuboTypes::Variable(resolveVar(si->e())));
    // Support only MIN and MAX.
    if (_objType != SolveI::ST_MIN && _objType != SolveI::ST_MAX) {
      return status;
    }
  }
  if (!calcQubo(si->e())) {
    return status;
  }
  status = SolverInstance::SAT;

#if 0
  {
    Printer p(std::cerr);
    std::cerr << "Flat--------\n";
    p.print(env().flat());
    std::cerr << "SolveItem--------\n";
    p.print(env().flat()->solveItem());
    std::cerr << "solveItem\n";
    p.print(si->e());
  }
#endif
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

unsigned long QuboVariable::global_index = 1;

SolverVariable::SolverVariable(Id* id) : _qv(nullptr) {
  _t = QuboTypes::Variable::Type::BOOL_TYPE;
  // std::cerr << "SolverVariable id: " << *id << "\n";
  if (id->decl()->ann().contains(constants().ann.is_defined_var)) {
    // std::cerr << "  is_defined_var\n";
  } else {
    _qv = new QuboVariable();
    // std::cerr << "  is QuboVariable\n";
  }
}

inline void QuboSolverInstance::insertVar(Id* id, QuboTypes::Variable var) {
  // std::cerr << *id << ": " << id->decl() << std::endl;
  _variableMap.insert(id->decl()->id(), var);
}

inline void QuboSolverInstance::dumpVar(void) {
  _log << "dumpVar----------------\n";
  for (auto const& pair : _variableMap) {
    _log << *pair.first << ": " << pair.second.index() << "\n";
  }
}

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
    auto* al = eval_array_lit(env().envi(), ad->e());
    return _variableMap.get((*al)[idx]->cast<Id>());
  }
  std::stringstream ssm;
  ssm << "Expected Id, VarDecl or ArrayAccess instead of \"" << *e << "\"";
  throw InternalError(ssm.str());
}

std::pair<double, bool> QuboSolverInstance::exprToConstEasy(Expression* e) {
  std::cerr << "----exprToConst(" << *e << ")\n";
  std::pair<double, bool> res{0.0, true};
  if (auto* il = e->dynamicCast<IntLit>()) {
    res.first = (static_cast<double>(il->v().toInt()));
  } else if (auto* fl = e->dynamicCast<FloatLit>()) {
    res.first = (fl->v().toDouble());
  } else if (auto* bl = e->dynamicCast<BoolLit>()) {
    res.first = static_cast<double>(bl->v());
  } else {
    res.second = false;
  }
  return res;
}

double QuboSolverInstance::exprToConst(Expression* e) {
  const auto e2ce = exprToConstEasy(e);
  if (!e2ce.second) {
    std::ostringstream oss;
    oss << "ExprToConst: expected a numeric/bool literal, getting " << *e;
    throw InternalError(oss.str());
  }
  return e2ce.first;
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
#ifdef QUBO_USE_MZN_DIRECTLY
  SolverConfig sc("org.minizinc.mzn-mzn", getVersion(nullptr));
#else
  SolverConfig sc("org.minizinc.qubo", getVersion(nullptr));
#endif
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
