/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/flat_exp.hh>

namespace MiniZinc {

EE flatten_fieldaccess(EnvI& env, const Ctx& ctx, Expression* e, VarDecl* r, VarDecl* b) {
  auto* fa = e->cast<FieldAccess>();
  assert(fa->v()->type().bt() == Type::BT_TUPLE);  // TODO: Support for Records

  // TODO: Does this have the correct behaviour?
  EE ret = flat_exp(env, ctx, fa->v(), r, b);
  auto* al = eval_array_lit(env, ret.r())->cast<ArrayLit>();

  assert(fa->field()->isa<IntLit>());
  IntVal i = fa->field()->cast<IntLit>()->v();
  ret.r = (*al)[i.toInt() - 1];

  return ret;
}

}  // namespace MiniZinc
