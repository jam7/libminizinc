### MiniZinc QUBO Solver Target

### Compile target for the QUBO interface
add_library(minizinc_qubo OBJECT
  solvers/QUBO/qubo_solverfactory.cpp
  solvers/QUBO/qubo_solverinstance.cpp

  include/minizinc/solvers/QUBO/qubo_solverfactory.hh
  include/minizinc/solvers/QUBO/qubo_solverinstance.hh
)
target_include_directories(minizinc_qubo PRIVATE "${QUBO_INCLUDE_DIRS}")
add_dependencies(minizinc_qubo minizinc_parser)

### Setup correct compilation into the MiniZinc library
target_compile_definitions(mzn PRIVATE HAS_QUBO)
target_sources(mzn PRIVATE $<TARGET_OBJECTS:minizinc_qubo>)
# No new objects to link
# target_link_libraries(mzn Qubo)
