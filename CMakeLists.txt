cmake_minimum_required(VERSION 3.25)

set(NAME cfd-unstructured-parallel)
project(${NAME} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)

include(FetchContent)

find_package(MPI REQUIRED)
find_package(OpenMP)


include(cmake/Dependencies.cmake)
#include(cmake/Documentation.cmake)

# Explicitly list your sources
add_executable(${NAME}
        src/main.cpp

        src/RunManager.cpp
        src/Simulation.cpp

        src/solver/SolverFactory.cpp
        src/solver/FiniteVolumeSolver.cpp
        include/solver/Solver.hpp
        src/solver/TimeStepCalculator.cpp
        src/solver/PositivityLimiter.cpp

        src/time/ForwardEulerTimeIntegrator.cpp
        src/time/SSPRK3TimeIntegrator.cpp
        src/time/SSPRK2TimeIntegrator.cpp
        include/time/TimeIntegrator.hpp

        src/spatial/GodunovSpatialOperator.cpp
        src/spatial/GodunovKolganRodionovSpatialOperator.cpp
        include/spatial/SpatialOperator.hpp

        src/reconstruction/P0Reconstruction.cpp
        src/reconstruction/P1Reconstruction.cpp
        include/reconstruction/Reconstruction.hpp

        src/riemann/RusanovRiemannSolver.cpp
        src/riemann/HLLRiemannSolver.cpp
        src/riemann/HLLCRiemannSolver.cpp
        #    src/riemann/RoeRiemannSolver.cpp
        #    src/riemann/ExactIdealGasRiemannSolver.cpp

        src/output/WriterFactory.cpp
        src/output/VTKWriter.cpp
        src/output/VTKRecomposer.cpp
        src/output/RunStatistics.cpp

        src/data/DataLayer.cpp
        src/data/Variables.cpp
        src/data/Workspace.cpp

        src/geometry/GmshMeshBuilder.cpp
        src/geometry/MortonOrder.cpp
        #    src/geometry/DelaunayMeshBuilder.cpp
        #    src/geometry/DelaunayTriangulator.cpp
        #    src/geometry/PointCloudBuilder.cpp
        #    src/geometry/StructuredMeshBuilder.cpp
        src/geometry/Mesh.cpp
        include/geometry/Node.hpp
        include/geometry/Cell.hpp
        src/geometry/Face.cpp

        src/parallel/MPIContext.cpp
        src/parallel/HaloExchange.cpp
        src/parallel/DomainDecomposition.cpp
        src/parallel/METISDecomposition.cpp
        src/parallel/RCBDecomposition.cpp
        src/parallel/MeshDistribution.cpp
        src/parallel/MeshSerialization.cpp
        src/parallel/HaloSerialization.cpp
        include/parallel/StateSynchronizer.hpp

        src/config/ConfigParser.cpp
        src/config/YamlConfigParser.cpp
        src/config/CliParser.cpp
        include/config/Settings.hpp
        include/config/InitialConditions.hpp
        src/config/InitialConditionInitializer.cpp

        src/bc/BoundaryFactory.cpp
        src/bc/BoundaryManager.cpp
        src/bc/FreeStreamBoundary.cpp
        src/bc/InletBoundary.cpp
        src/bc/OutletBoundary.cpp
        src/bc/PeriodicBoundary.cpp
        src/bc/ReflectiveBoundary.cpp
        src/bc/NonReflectiveBoundary.cpp
        src/bc/SymmetryBoundary.cpp
        src/bc/WallBoundary.cpp
        include/bc/BoundaryCondition.hpp

        src/utils/StringUtils.cpp
)

target_include_directories(${NAME}
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/external
)

target_link_libraries(${NAME}
        PRIVATE
        xtensor
        yaml-cpp::yaml-cpp
        cxxopts::cxxopts

        VTK::CommonCore
        VTK::CommonDataModel
        VTK::FiltersCore
        VTK::FiltersGeneral
        VTK::IOCore
        VTK::IOLegacy
        VTK::IOXML

        MPI::MPI_CXX

        OpenMP::OpenMP_CXX

        gmsh::gmsh
	
	GKlib
	
	metis::metis
)

target_compile_definitions(${NAME} PRIVATE XTENSOR_USE_OPENMP)
