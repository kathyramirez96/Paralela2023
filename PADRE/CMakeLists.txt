# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define sources
file(GLOB_RECURSE SOURCES_FILES
    "classification/*.cpp"
    "classification/*.cxx"
)

file(GLOB_RECURSE METRICS_SOURCES_FILES
    "metrics/*.cpp"
)

# Define headers
file(GLOB_RECURSE HEADERS_FILES
    "*.h"
    "*.hpp"
)

# Set sources
add_library(classification STATIC "${SOURCES_FILES}")
add_library(classificationMetrics STATIC "${METRICS_SOURCES_FILES}")


# Include dirs
include_directories(classification)
include_directories(metrics)



# C++ flags
set(USE_PTHREAD ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
