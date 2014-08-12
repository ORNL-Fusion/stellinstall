cmake_minimum_required (VERSION 2.8)

# Add subdirectories for all the sources.
add_subdirectory (General)
add_subdirectory (Hessian)
add_subdirectory (Initialization_Cleanup)
add_subdirectory (Input_Output)
add_subdirectory (NESTOR_vacuum)
add_subdirectory (Reconstruction)
add_subdirectory (Splines)
add_subdirectory (TimeStep)

set (vmec_sources "${vmec_sources}" PARENT_SCOPE)