# Standalone RunStorage read/write-mixing reproducer (debug tool)
add_executable(rs_rwtest rs_rwtest.cpp)
target_link_libraries(rs_rwtest
  rm_serial
  rm_yamr
  rm_external
)
