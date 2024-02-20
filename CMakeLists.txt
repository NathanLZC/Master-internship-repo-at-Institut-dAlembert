cmake_minimum_required(VERSION 3.10)

# 设置项目名和版本
project(Master-internship-repo-at-Institut-dAlembert VERSION 1.0)

# 指定C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 导入FFTW库
set(FFTW_DIR "/gagarine/temporaires/zli/fftw")
find_library(FFTW_LIBRARY NAMES fftw3 PATHS "${FFTW_DIR}/lib" NO_DEFAULT_PATH)

# 添加包含目录
include_directories("${FFTW_DIR}/include")

# 定义目标可执行文件和相关源文件
add_executable(fftw test_for_fftw.cc)

# 链接FFTW库到目标可执行文件
target_link_libraries(fftw ${FFTW_LIBRARY})

