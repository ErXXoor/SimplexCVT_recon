include(FetchContent)
FetchContent_Declare(
        geogram
        GIT_REPOSITORY https://github.com/BrunoLevy/geogram.git
        GIT_TAG v1.9.4
)
FetchContent_MakeAvailable(geogram)
target_compile_features(geogram PUBLIC cxx_std_11)
target_compile_options(geogram PUBLIC -Wimplicit-int -Wimplicit-function-declaration)