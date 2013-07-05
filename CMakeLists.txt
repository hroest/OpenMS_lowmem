### example CMakeLists.txt to develop programs using OpenMS
project("Main")
cmake_minimum_required(VERSION 2.6)

## list all your executables here (a corresponding .C file should exist, e.g. Main.C)
set(my_executables
  LowMemGaussFilter
  LowMemPeakPickerHiRes
  LowMemFileCacher
)

## find OpenMS configuration and register target "OpenMS" (our library)
find_package(OpenMS)
## if the above fails you can try:
#find_package(OpenMS PATHS "</c/path/to/OpenMS>/cmake")

# check whether the OpenMS package was found
if (OpenMS_FOUND)

  ## include directories for OpenMS headers (and contrib)
  include_directories(${OPENMS_INCLUDE_DIRS})

  
  ## append precompiler macros and compiler flags specific to OpenMS
  ## Warning: this could be harmful to your project. Check this if problems occur.
  ## Also, use this to add your own compiler flags, e.g. for OpenMP support.
  ## e.g. for Visual Studio use /openmp
  ## set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS} /openmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS}")
  add_definitions(${OPENMS_DEFINITIONS})
  
  ## add targets for the executables
  foreach(i ${my_executables})
    add_executable(${i} ${i}.C)
    ## link executables against OpenMS
  target_link_libraries(${i} OpenMS)
  endforeach(i)

  
else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)

set(TMPDIR /tmp/)
set(DIFFCMD FuzzyDiff)

ADD_TEST("LowMemGaussFilter_test_1" LowMemGaussFilter -in testdata/NoiseFilterGaussian_1_input.mzML -ini testdata/NoiseFilterGaussian_1_parameters.ini -out ${TMPDIR}/NoiseFilterGaussian_1_output.mzML.tmp -test)
ADD_TEST("LowMemGaussFilter_test_1_out" ${DIFFCMD} -in1 ${TMPDIR}/NoiseFilterGaussian_1_output.mzML.tmp -in2 testdata/NoiseFilterGaussian_1_output.mzML -ini testdata/FuzzyDiff.ini)
ADD_TEST("LowMemPeakPickerHiRes_test_1" LowMemPeakPickerHiRes -in testdata/PeakPickerHiRes_input.mzML -ini testdata/PeakPickerHiRes_parameters.ini -out ${TMPDIR}/PeakPickerHiRes_output.mzML.tmp -test)
ADD_TEST("LowMemPeakPickerHiRes_test_1_out" ${DIFFCMD} -in1 ${TMPDIR}/PeakPickerHiRes_output.mzML.tmp -in2 testdata/PeakPickerHiRes_output.mzML -ini testdata/FuzzyDiff.ini)
ADD_TEST("LowMemFileCacher_test_1" LowMemFileCacher -in testdata/NoiseFilterGaussian_1_input.mzML -out ${TMPDIR}/FileCacher_1_output.mzML.tmp -test)

ENABLE_TESTING()

