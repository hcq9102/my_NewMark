#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "HPX::component_storage_component" for configuration "Release"
set_property(TARGET HPX::component_storage_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::component_storage_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libhpx_component_storage.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_component_storage.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::component_storage_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::component_storage_component "${_IMPORT_PREFIX}/lib64/libhpx_component_storage.so.1.9.0" )

# Import target "HPX::unordered_component" for configuration "Release"
set_property(TARGET HPX::unordered_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::unordered_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libhpx_unordered.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_unordered.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::unordered_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::unordered_component "${_IMPORT_PREFIX}/lib64/libhpx_unordered.so.1.9.0" )

# Import target "HPX::partitioned_vector_component" for configuration "Release"
set_property(TARGET HPX::partitioned_vector_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::partitioned_vector_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libhpx_partitioned_vector.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_partitioned_vector.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::partitioned_vector_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::partitioned_vector_component "${_IMPORT_PREFIX}/lib64/libhpx_partitioned_vector.so.1.9.0" )

# Import target "HPX::iostreams_component" for configuration "Release"
set_property(TARGET HPX::iostreams_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::iostreams_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libhpx_iostreams.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_iostreams.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::iostreams_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::iostreams_component "${_IMPORT_PREFIX}/lib64/libhpx_iostreams.so.1.9.0" )

# Import target "HPX::parcel_coalescing" for configuration "Release"
set_property(TARGET HPX::parcel_coalescing APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::parcel_coalescing PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/hpx/libhpx_parcel_coalescing.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_parcel_coalescing.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::parcel_coalescing )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::parcel_coalescing "${_IMPORT_PREFIX}/lib64/hpx/libhpx_parcel_coalescing.so.1.9.0" )

# Import target "HPX::io_counters_component" for configuration "Release"
set_property(TARGET HPX::io_counters_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::io_counters_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/hpx/libhpx_io_counters.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_io_counters.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::io_counters_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::io_counters_component "${_IMPORT_PREFIX}/lib64/hpx/libhpx_io_counters.so.1.9.0" )

# Import target "HPX::memory_counters_component" for configuration "Release"
set_property(TARGET HPX::memory_counters_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::memory_counters_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/hpx/libhpx_memory_counters.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_memory_counters.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::memory_counters_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::memory_counters_component "${_IMPORT_PREFIX}/lib64/hpx/libhpx_memory_counters.so.1.9.0" )

# Import target "HPX::process_component" for configuration "Release"
set_property(TARGET HPX::process_component APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HPX::process_component PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libhpx_process.so.1.9.0"
  IMPORTED_SONAME_RELEASE "libhpx_process.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS HPX::process_component )
list(APPEND _IMPORT_CHECK_FILES_FOR_HPX::process_component "${_IMPORT_PREFIX}/lib64/libhpx_process.so.1.9.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
