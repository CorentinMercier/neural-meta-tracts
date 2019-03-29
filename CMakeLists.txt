cmake_minimum_required(VERSION 3.9)
project(fiberreductor)

#Module dependencies
find_package(Qt4 REQUIRED QtGui QtOpenGL)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenMP)

find_package(VTK)
if(VTK_FOUND)
    message("VTK Found")
    include(${VTK_USE_FILE})
    set(USE_VTK 1)
else()
    message("VTK not found : you can use tck files for fibers")
    message("The cortical surface is not usable without vtk")
endif()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/sources/libraries.h.in ${CMAKE_CURRENT_BINARY_DIR}/libraries.h)

include_directories(lib/eigen)

#Source files
file(GLOB_RECURSE source_files sources/*.cpp sources/*.h sources/*.cxx)
QT4_WRAP_CPP(MOC_GENERATED sources/interface/mainwindow.h sources/interface/widgetopengl.h)
QT4_WRAP_UI(UI_GENERATED sources/interface/mainwindow.ui)

#Target to build
add_executable(fiberreductor ${source_files} ${UI_GENERATED} ${MOC_GENERATED})
set_property(TARGET fiberreductor PROPERTY CXX_STANDARD 11)
target_include_directories(fiberreductor
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR} # ui & moc files will be generated here
)

target_link_libraries(fiberreductor
        PUBLIC OpenMP::OpenMP_CXX
        PUBLIC OpenGL::GL
        PUBLIC OpenGL::GLU
        PUBLIC Qt4::QtGui
        PUBLIC Qt4::QtOpenGL
        PUBLIC GLEW::GLEW
)
if (VTK_FOUND)
    target_link_libraries(fiberreductor PUBLIC ${VTK_LIBRARIES})
endif()

# Additional flags
if(NOT MSVC)
        target_compile_options(fiberreductor PUBLIC -Wno-unused-result -O3)
endif()


