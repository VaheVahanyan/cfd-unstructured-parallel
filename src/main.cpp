#include <exception>
#include <iostream>

#include "RunManager.hpp"
#include "parallel/MPIContext.hpp"

namespace
{
    class MPIInitializerGuard {
    public:
        MPIInitializerGuard(int& argc, char**& argv) {
            MPIContext::Initialize(argc, argv);
        }

        ~MPIInitializerGuard() {
            MPIContext::Finalize();
        }

        MPIInitializerGuard(const MPIInitializerGuard&) = delete;
        MPIInitializerGuard& operator=(const MPIInitializerGuard&) = delete;
    };

    void AbortMPIIfAlive() {
        if (MPIContext::IsInitialized() && !MPIContext::IsFinalized()) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}  // namespace

auto main(int argc, char* argv[]) -> int {
    MPIInitializerGuard mpi_guard(argc, argv);

    try {
        RunManager run_manager;
        return run_manager.Run(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        AbortMPIIfAlive();
        return 1;
    }
    catch (...) {
        std::cerr << "Fatal error: unknown exception\n";
        AbortMPIIfAlive();
        return 1;
    }
}
