cmake_minimum_required(VERSION 3.10)
project(func2svg)

# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找OpenCV库
find_package(OpenCV REQUIRED)

# 查找Eigen库
find_package(Eigen3 REQUIRED)

find_package(PkgConfig REQUIRED)

# Use pkg-config to find Graphviz (libgvc)
pkg_check_modules(GRAPHVIZ REQUIRED libgvc)

# Include Graphviz header directories
include_directories(${GRAPHVIZ_INCLUDE_DIRS})
# dot
#find_program(DOT_EXEC dot)
#if(NOT DOT_EXEC)
#    message(FATAL_ERROR "Graphviz not found. Please install Graphviz.")
#endif()
#
#set(GRAPHICS_DIR ${CMAKE_BINARY_DIR}/graphics)
#add_custom_command(
#        OUTPUT ${GRAPHICS_DIR}/diagram.svg
#        COMMAND ${DOT_EXEC} -Tsvg ${CMAKE_CURRENT_SOURCE_DIR}/diagram.dot -o ${GRAPHICS_DIR}/diagram.svg
#        DEPENDS diagram.dot
#)
#add_custom_target(generate_graphics ALL DEPENDS ${GRAPHICS_DIR}/diagram.svg)

# 假设ALGLIB是本地库，手动指定路径
set(ALGLIB_DIR ${CMAKE_SOURCE_DIR}/src)
include_directories(${ALGLIB_DIR}/alglib)
file(GLOB ALGLIB_SOURCES ${ALGLIB_DIR}/alglib/src/*.cpp)
add_library(alglib STATIC ${ALGLIB_SOURCES})


# 包含本地头文件
include_directories(${PROJECT_SOURCE_DIR}/include)

# 添加可执行文件
add_executable(func2svg
    src/QuasiRegularModel.cpp
    src/ImageProcess.cpp
    src/topology.cpp
    src/ImageRestore.cpp
    main.cpp
)

# 链接OpenCV库
target_link_libraries(func2svg ${OpenCV_LIBS})

# 链接Eigen库
target_link_libraries(func2svg Eigen3::Eigen)

# 链接ALGLIB库
target_link_libraries(func2svg alglib)

target_link_libraries(func2svg ${GRAPHVIZ_LIBRARIES})