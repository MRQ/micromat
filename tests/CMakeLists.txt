find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(micromat-tests
	TestEigen.cpp
	TestMicromat.cpp
	TestRunner.cpp
)

target_link_libraries( micromat-tests PUBLIC
	micromat
	Eigen3::Eigen
	Boost::unit_test_framework
)
