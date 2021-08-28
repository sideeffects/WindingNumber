cmake_minimum_required(VERSION 3.16)

project(WindingNumber)

list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/recipes)

# WindingNumber library
file(GLOB SRC_FILES *.cpp *.h )
add_library(WindingNumber ${SRC_FILES})
add_library(WindingNumber::WindingNumber ALIAS WindingNumber)

target_include_directories(WindingNumber PUBLIC .)

include(onetbb)
include(simde)
target_link_libraries(WindingNumber PUBLIC TBB::tbb simde::simde)

set_target_properties(WindingNumber PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_features(WindingNumber PUBLIC cxx_std_11)
