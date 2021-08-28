if(TARGET simde::simde)
    return()
endif()

message(STATUS "Third-party: creating target 'simde::simde'")

include(FetchContent)
FetchContent_Declare(
    simde
    GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
    GIT_TAG 48edfa906d835525e2061fbf6062b7c326d66840
)
FetchContent_MakeAvailable(simde)

add_library(simde::simde INTERFACE IMPORTED GLOBAL)
target_include_directories(simde::simde INTERFACE "${simde_SOURCE_DIR}")

# Enables native aliases. Not ideal but makes it easier to convert old code.
target_compile_definitions(simde::simde INTERFACE SIMDE_ENABLE_NATIVE_ALIASES)

# Uncomment this line to ensure code can be compiled without native SIMD (i.e. emulates everything)
# target_compile_definitions(simde::simde INTERFACE SIMDE_NO_NATIVE)
