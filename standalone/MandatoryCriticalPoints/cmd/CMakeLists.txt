cmake_minimum_required(VERSION 3.2)

project(ttkMandatoryCriticalPointsCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(ttkMandatoryCriticalPointsCmd main.cpp)
target_link_libraries(ttkMandatoryCriticalPointsCmd ttkMandatoryCriticalPoints
	ttkProgramBase)

install(TARGETS ttkMandatoryCriticalPointsCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

