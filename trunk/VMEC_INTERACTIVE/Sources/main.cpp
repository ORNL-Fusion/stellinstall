#include <iostream>

#include <unistd.h> // for sleep()

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <string>

#include <mpi.h>

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
using vmec_interactive::VmecInteractive;

/** setup the interactive VMEC */
extern "C" void vi_setup(void);

/** set an input parameter for the interactive VMEC */
extern "C" void vi_set(double* centralPressure);

/** reconverge and get an input or output parameter from the interactive VMEC */
extern "C" void vi_get(int* error);

/** exit the interactive VMEC */
extern "C" void vi_exit(void);

extern "C" void get_world_comm_from_fortran(int* world_comm);


class VmecInteractiveImpl final : public VmecInteractive::Service {
  
public:
  
  explicit VmecInteractiveImpl(void) {

  }

  Status Setup(ServerContext* context, const InitialInputFile* request, StatusMsg* reply) override {



    return Status::OK;
  }


  Status SetParameter(ServerContext* context, const SetParameterMsg* request, StatusMsg* reply) override {

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
  server->Wait();
  
  delete service;
}

int main(int argc, char** argv) {
  
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  
  MPI_Comm comm_world = MPI_COMM_WORLD;
  
  // Get the number of processes
  int world_size;
  MPI_Comm_size(comm_world, &world_size);
  
  // Get the rank of the process
  int world_rank;
    MPI_Comm_rank(comm_world, &world_rank);
  
    std::cout << "world size: " << world_size << " rank: " << world_rank << std::endl;
    
    if (world_rank == 0) {
        RunServer();
    }
    
    MPI_Barrier(comm_world);
    
    
    MPI_Finalize();
    
//     // setup the interactive VMEC service and converge initially
//     vi_setup();
//     
//     int world_comm_from_fortran = -1;
//     get_world_comm_from_fortran(&world_comm_from_fortran);
//   
//     MPI_Comm comm_world = MPI_Comm_f2c(world_comm_from_fortran);
//   
//     
//     
//     
//     
//         
//     int error = 0;
//     double centralPressure = 0.0;
//     for (int i=0; i<10; ++i) {
//     
//         centralPressure = 1000.0 + i;
//         vi_set(&centralPressure);
//         
//         vi_get(&error);
//         
//         MPI_Barrier(comm_world);
//     }
// 
//   
//     vi_exit();
  
    return 0;
}
