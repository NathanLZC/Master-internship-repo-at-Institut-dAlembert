cmake_minimum_required(VERSION 3.10)

# 设置项目名和版本
project(Master-internship-repo-at-Institut-dAlembert VERSION 1.0)

# 指定C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 导入FFTW库
set(FFTW_DIR "/gagarine/temporaires/zli/fftw")
set(OpenBLAS_DIR "/gagarine/temporaires/zli/OpenBLAS")
set(Eigen_DIR "/gagarine/temporaires/zli/eigen-3.4.0")
#find_library(FFTW_LIBRARY NAMES fftw3 PATHS "${FFTW_DIR}/lib" NO_DEFAULT_PATH)

# 添加包含目录#eigen is a header only library, so no need to link it
include_directories("${FFTW_DIR}/include" "${OpenBLAS_DIR}/include" "${Eigen_DIR}")
link_directories("${FFTW_DIR}/lib" "${OpenBLAS_DIR}/lib")
link_libraries(fftw3 openblas)

# 定义目标可执行文件和相关源文件
#add_executable(our_test hello_world_for_cpp_configure/test_for_fftw.cc)
add_executable(test_eigen Week5/test_eigen.cc)

# 链接FFTW库到目标可执行文件
#target_link_libraries(our_test fftw3)

