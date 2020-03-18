cmake_minimum_required(VERSION 3.12)
project(fiberreductor LANGUAGES CXX)

# Dependencies
# -------------------------------------------------------------------------------------------------
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenMP)
add_subdirectory(lib/glad)
#Module dependencies
find_package(VTK)
if(VTK_FOUND)
	message("VTK Found")
	include(${VTK_USE_FILE})
	set(USE_VTK 1)
else()
	message("VTK not found : you can use tck files for fibers")
	message("The cortical surface is not usable without vtk")
endif()

# Source files
# -------------------------------------------------------------------------------------------------
file(GLOB_RECURSE source_files sources/*.cpp sources/*.h sources/*.cxx)
QT5_WRAP_CPP(MOC_GENERATED sources/interface/mainwindow.h sources/interface/widgetopengl.h)
QT5_WRAP_UI(UI_GENERATED sources/interface/mainwindow.ui)
source_group(TREE ${CMAKE_CURRENT_LIST_DIR} FILES ${source_files} ${UI_GENERATED} ${MOC_GENERATED})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/sources/libraries.h.in ${CMAKE_CURRENT_BINARY_DIR}/libraries.h)

# fiberreductor main target
# -------------------------------------------------------------------------------------------------
add_executable(fiberreductor ${source_files} ${UI_GENERATED} ${MOC_GENERATED})
target_compile_features(fiberreductor PUBLIC cxx_std_11)
target_include_directories(fiberreductor
	PUBLIC lib/glm
	PUBLIC lib/eigen
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR} # ui & moc files will be generated here
)

# convert the fpic option in Qt5::Core over to INTERFACE_POSITION_INDEPENDENT_CODE
# -------------------------------------------------------------------------------------------------
if ( TARGET Qt5::Core )
	get_property( core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS )
	string( REPLACE "-fPIC" "" new_core_options "${core_options}" )
	set_property( TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options} )
	set_property( TARGET Qt5::Core PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE "ON" )
	set( CMAKE_CXX_COMPILE_OPTIONS_PIE "-fPIC" )
endif()

#Target to build
target_link_libraries(fiberreductor
		PUBLIC OpenMP::OpenMP_CXX
		PUBLIC OpenGL::GL
		PUBLIC OpenGL::GLU
		PUBLIC Qt5::Widgets
		PUBLIC Qt5::OpenGL
		PUBLIC glad
)
if (VTK_FOUND)
	target_link_libraries(fiberreductor PUBLIC ${VTK_LIBRARIES})
endif()

# Additional flags
if(NOT MSVC)
		target_compile_options(fiberreductor PUBLIC -Wno-unused-result -O3)
endif()
