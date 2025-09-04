# Enable ctest framework and add tests
enable_testing()

# Some exemples requires NetCDF
#find_package(NetCDF REQUIRED Fortran)

# Find NetCDF using nc-config
execute_process(
    COMMAND nc-config --prefix
    OUTPUT_VARIABLE NetCDF_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND nc-config --includedir
    OUTPUT_VARIABLE NetCDF_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND nc-config --libs
    OUTPUT_VARIABLE NetCDF_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(NetCDF_LIBRARIES "${NetCDF_LIBRARIES} -lnetcdff")


list ( APPEND example_srcs
      low_pass_filter_on_the_sphere.f90
      #mcmc_ensemble_update.F90
      #mcmc_ensemble_update_obs100err20loc06.h90
      #mcmc_ensemble_update_obs100err20loc06x2.h90
      #mcmc_ensemble_update_obs100err20loc06x4.h90
      #mcmc_ensemble_update_parameters.h90
      #observation_regression_on_the_sphere.f90
      random_field_on_the_sphere.f90
      scores_idealized_example.F90
)

set ( LIB_TARGET ensdam_static )

foreach( src ${example_srcs} )

    if ( WIN32 )
       string( REPLACE ".f90" "" exec_name ${src} )
    else()
       string( REPLACE ".F90" ".x" exec_name ${src} )
       string( REPLACE ".f90" ".x" exec_name ${src} )
    endif()
    add_executable( ${exec_name} EXCLUDE_FROM_ALL ${src} )
    add_dependencies( ${exec_name} ${LIB_TARGET} )
    # target_link_libraries( ${exec_name} PUBLIC ${LIB_TARGET} )
    target_link_libraries(${exec_name} PUBLIC ${LIB_TARGET} ${NetCDF_LIBRARIES})
    # Include NetCDF headers
    target_include_directories(${exec_name} PRIVATE ${NetCDF_INCLUDE_DIRS})

    add_test( NAME ${exec_name} COMMAND ${exec_name} )
    list(APPEND examples ${exec_name})

    #set_target_properties( ${exec_name} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE )
    #set_target_properties( ${exec_name} PROPERTIES SKIP_BUILD_RPATH         FALSE )

    string( REPLACE ".f90" ""   test_name ${src} )

endforeach()
add_custom_target(examples)
set_property( TARGET examples PROPERTY FOLDER examples )
add_dependencies(examples ${examples})
