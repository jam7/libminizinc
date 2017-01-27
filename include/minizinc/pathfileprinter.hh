/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Guido Tack <guido.tack@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __MINIZINC_HTMLPRINTER_HH
#define __MINIZINC_HTMLPRINTER_HH

#include <iostream>
#include <minizinc/ast.hh>
#include <minizinc/stl_map_set.hh>
#include <utility>

namespace MiniZinc {

  class Model;
  
  class PathFilePrinter {

  typedef std::pair<std::string, std::string> NamePair;
  typedef UNORDERED_NAMESPACE::unordered_map<Id*, NamePair > NameMap;

  NameMap betternames;
  private:
    std::ostream& os;
    EnvI& ei;
    bool remove_paths;

    void addBetterName(Id* id, std::string name, std::string path, bool overwrite);

  public:
    PathFilePrinter(std::ostream& o, EnvI& envi, bool rem = true);
    void print(Model* m);
    void print(Item* i);
  };
  
}

#endif