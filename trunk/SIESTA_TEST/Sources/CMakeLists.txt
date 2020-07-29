# Append all source files to variable siesta_test_sources. As new files are
# added this must be updated.
list (APPEND siesta_test_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/siesta_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/grid_quantity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/grid_quantity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/siesta_quantities.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vmec_quantities.hpp
)
set (siesta_test_sources "${siesta_test_sources}" PARENT_SCOPE)
