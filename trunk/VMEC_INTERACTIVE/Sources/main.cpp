#include <iostream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <string>
#include <future>
#include <thread>
#include <csignal>

#include <unistd.h> // for sleep()
#include <mpi.h>

// setup the interactive VMEC
extern "C" void vi_setup(int* mpi_comm);

// set an input parameter for the interactive VMEC
extern "C" void vi_set(double* centralPressure);

// reconverge and get an input or output parameter from the interactive VMEC
extern "C" void vi_get(int* error);

// exit the interactive VMEC
extern "C" void vi_exit(void);


int command = 0;
const int CMD_SETUP    = 1;
const int CMD_SET      = 2;
const int CMD_GET      = 3;
const int CMD_SHUTDOWN = 4;

double centralPressure = 0.0;

int mpi_comm_f;

// Somewhere in the global scope :/
std::promise<void> exit_requested;

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

class VmecInteractiveImpl final : public VmecInteractive::Service {
public:
    explicit VmecInteractiveImpl(MPI_Comm new_mpi_comm) {
        mpi_comm = new_mpi_comm;
        MPI_Comm_size(mpi_comm, &mpi_size);
        MPI_Comm_rank(mpi_comm, &mpi_rank);
        
        mpi_comm_f = MPI_Comm_c2f(mpi_comm);
    }
    
    Status Setup(ServerContext* context, const InitialInputFile* request, StatusMsg* reply) override {
        command = CMD_SETUP;
        MPI_Bcast(&command, 1, MPI_INTEGER, 0, mpi_comm);
        
        printf("setup %d\n", mpi_rank);
        
        // setup the interactive VMEC service and converge initially
        vi_setup(&mpi_comm_f);
        
        return Status::OK;
    }
    
    Status SetParameter(ServerContext* context, const SetParameterMsg* request, StatusMsg* reply) override {
        command = CMD_SET;
        MPI_Bcast(&command, 1, MPI_INTEGER, 0, mpi_comm);
        
        printf("set %d\n", mpi_rank);
        
        centralPressure = request->parametervalue();
        MPI_Bcast(&centralPressure, 1, MPI_DOUBLE, 0, mpi_comm);
        
        vi_set(&centralPressure);
        
        return Status::OK;
    }
    
    Status GetParameter(ServerContext* context, const GetParameterMsg* request, ParameterMsg* reply) override {
        command = CMD_GET;
        MPI_Bcast(&command, 1, MPI_INTEGER, 0, mpi_comm);
        
        printf("get %d\n", mpi_rank);
        
        int error = 0;
        vi_get(&error);
        
        return Status::OK;
    }
    
    Status Shutdown(ServerContext* context, const ShutdownMsg* request, StatusMsg* reply) override {
        command = CMD_SHUTDOWN;
        MPI_Bcast(&command, 1, MPI_INTEGER, 0, mpi_comm);
        
        printf("shutdown %d\n", mpi_rank);
        
        vi_exit();
        
        exit_requested.set_value();
        return Status::OK;
    }
    
private:
    MPI_Comm mpi_comm;
    int mpi_rank;
    int mpi_size;
    int mpi_comm_f;
};


void RunServer(MPI_Comm mpi_comm) {
    
    std::string server_address("0.0.0.0:50051");
    
    VmecInteractiveImpl* service = new VmecInteractiveImpl(mpi_comm);
    
    ServerBuilder builder;
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
    builder.RegisterService(service);
    std::unique_ptr<Server> server(builder.BuildAndStart());
    
    std::cout << "Interactive VMEC listening on " << server_address << std::endl;
    //server->Wait();
    
    // shutdown gRPC service from client side; see https://stackoverflow.com/a/36182429
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
   
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    
    // translate communicator to Fortran
    mpi_comm_f = MPI_Comm_c2f(mpi_comm);

    int mpi_size;
    int mpi_rank;
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    // printf("MPI size=%d rank=%d\n", mpi_size, mpi_rank);
    
    if (mpi_rank == 0) {
        // run webservice in main rank
        RunServer(mpi_comm);
        
    } else {
        // all other ranks wait for commands and forward them to other VMEC tasks
        int cancel = 0;
        while(!cancel) {
            
            MPI_Bcast(&command, 1, MPI_INTEGER, 0, mpi_comm);
            
            int error = 0;
            
            switch(command) {
                case CMD_SETUP:
                    printf("setup %d\n", mpi_rank);
                    
                    // setup the interactive VMEC service and converge initially
                    vi_setup(&mpi_comm_f);
                    
                    break;
                case CMD_SET:
                    printf("set %d\n", mpi_rank);
                    
                    MPI_Bcast(&centralPressure, 1, MPI_DOUBLE, 0, mpi_comm);
                    
                    vi_set(&centralPressure);
                    
                    break;
                case CMD_GET:
                    printf("get %d\n", mpi_rank);
                    
                    vi_get(&error);
                    
                    break;
                case CMD_SHUTDOWN:
                    printf("shutdown %d\n", mpi_rank);
                    
                    vi_exit();
                    
                    cancel = 1;
                    break;
                default:
                    break;
            }
        }
    }
    
    MPI_Finalize();

    return 0;
}
