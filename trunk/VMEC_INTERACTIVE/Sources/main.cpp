#include <iostream>

#include <unistd.h> // for sleep()

#include <mpi.h>


#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <string>

#include <future>
#include <thread>
#include <csignal>

#include <grpc/grpc.h>
#include <grpcpp/server.h>
#include <grpcpp/server_builder.h>
#include <grpcpp/server_context.h>
#include <grpcpp/security/server_credentials.h>

#include "vmec_interactive.grpc.pb.h"

using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::ServerReader;
using grpc::ServerReaderWriter;
using grpc::ServerWriter;
using grpc::Status;

using vmec_interactive::InitialInputFile;
using vmec_interactive::StatusMsg;
using vmec_interactive::SetParameterMsg;
using vmec_interactive::GetParameterMsg;
using vmec_interactive::ParameterMsg;
using vmec_interactive::ShutdownMsg;
using vmec_interactive::VmecInteractive;

/** setup the interactive VMEC */
extern "C" void vi_setup(int* mpi_comm);

/** set an input parameter for the interactive VMEC */
extern "C" void vi_set(double* centralPressure);

/** reconverge and get an input or output parameter from the interactive VMEC */
extern "C" void vi_get(int* error);

/** exit the interactive VMEC */
extern "C" void vi_exit(void);

extern "C" void get_world_comm_from_fortran(int* world_comm);

int cancel = 0;
MPI_Comm mpi_comm;
int mpi_comm_f;

// Somewhere in the global scope :/
std::promise<void> exit_requested;

class VmecInteractiveImpl final : public VmecInteractive::Service {
  
public:
  
  explicit VmecInteractiveImpl(void) {

  }

  Status Setup(ServerContext* context, const InitialInputFile* request, StatusMsg* reply) override {

    cancel = 1;
    MPI_Bcast(&cancel, 1, MPI_INTEGER, 0, mpi_comm);

    return Status::OK;
  }


  Status SetParameter(ServerContext* context, const SetParameterMsg* request, StatusMsg* reply) override {

    return Status::OK;
  }
  
  Status GetParameter(ServerContext* context, const GetParameterMsg* request, ParameterMsg* reply) override {
      
    return Status::OK;
  }
  
  Status Shutdown(ServerContext* context, const ShutdownMsg* request, StatusMsg* reply) override {
      
      cancel = 1;
      MPI_Bcast(&cancel, 1, MPI_INTEGER, 0, mpi_comm);
      

      printf("got shutdown request\n");
      exit_requested.set_value();
      
      return Status::OK;
  }
};





void RunServer() {
  std::string server_address("0.0.0.0:50051");
  VmecInteractiveImpl* service = new VmecInteractiveImpl();
  
  ServerBuilder builder;
  builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
  builder.RegisterService(service);
  std::unique_ptr<Server> server(builder.BuildAndStart());
  std::cout << "Interactive VMEC listening on " << server_address << std::endl;
  //server->Wait();
  
  auto handler = [](int s) {
    exit_requested.set_value();
    };
    std::signal(SIGINT, handler);
    std::signal(SIGTERM, handler);
    std::signal(SIGQUIT, handler);
  
  auto serveFn = [&]() {
    server->Wait();
  };

  std::thread serving_thread(serveFn);

  auto f = exit_requested.get_future();
  f.wait();
  server->Shutdown();
  serving_thread.join();
  
  delete service;
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
   
    mpi_comm = MPI_COMM_WORLD;
    
    // translate communicator to Fortran
    mpi_comm_f = MPI_Comm_c2f(mpi_comm);

    int mpi_size;
    int mpi_rank;
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    // printf("MPI size=%d rank=%d\n", mpi_size, mpi_rank);
    
    if (mpi_rank == 0) {
        // run webservice in main rank
        RunServer();
        
    } else {
        // all other ranks wait for commands and forward them to other VMEC tasks
        int saidReady = 0;
        while(!cancel) {
            if (!saidReady) {
                printf("rank %d ready and awaiting commands\n", mpi_rank);
                saidReady = 1;
            }
            
            
            
            
            
            
            MPI_Bcast(&cancel, 1, MPI_INTEGER, 0, mpi_comm);
            if (cancel) {
                printf("%d got cancel!\n", mpi_rank);
            }
        }
    }
    
    
    
    
    
    
    
    
//    
//     // setup the interactive VMEC service and converge initially
//     vi_setup(&mpi_comm_f);
//     
//     int error = 0;
//     double centralPressure = 0.0;
//     for (int i=0; i<3; ++i) {
//     
//         centralPressure = (i+1.0)*1000.0;
//         vi_set(&centralPressure);
//         
//         vi_get(&error);
//     }
// 
//   
//     vi_exit();
    
    MPI_Finalize();

    return 0;
}
