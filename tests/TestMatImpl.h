#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <boost/format.hpp>

inline void ThrowingAssert(const bool test, const int line,
		const char* file, const char* expr)
{
	if(test)
		return;

	throw std::runtime_error((boost::format("Assert \"%s\" failed in %s:%d.")
			% expr % file % line).str());
};

#if defined TEST_MAT_EIGEN
#define eigen_assert(expr) ThrowingAssert((expr), __LINE__, __FILE__, #expr)
#include <Eigen/Core>
#define TEST_MAT_NAMESPACE Eigen
#define TEST_MAT_CLASS TestEigen
#define TEST_MAT_PREFIX(CASE) Eigen ## CASE

#elif defined TEST_MAT_MICROMAT
#define micromat_assert(expr) ThrowingAssert((expr), __LINE__, __FILE__, #expr)
#include <micromat/micromat.h>
#define TEST_MAT_NAMESPACE umat
#define TEST_MAT_CLASS TestMicromat
#define TEST_MAT_PREFIX(CASE) Micromat ## CASE

#else
#error This file may only be used by means of TestEigen or TestMicromat
#endif

BOOST_AUTO_TEST_SUITE(TEST_MAT_CLASS)


BOOST_AUTO_TEST_CASE(TEST_MAT_PREFIX(Blocks))
{
	namespace mat = TEST_MAT_NAMESPACE;
	mat::Array<float, 3, 3> m33 = mat::Array<float, 3, 3>::Constant(1.5);
	BOOST_CHECK_EQUAL(m33.sum(), 1.5 * 9.);
	//m33.block<3,4>(0, 0); does not compile
	//m33.block<4,3>(0, 0); does not compile
	m33.block<3,3>(0, 0).setConstant(0.);
	m33.block<1,1>(1, 1)(0, 0) = 1;
	m33 += mat::Array<float, 3, 3>::Constant(-1);
	BOOST_CHECK_EQUAL(m33.abs2().sum(), 8.);
}

BOOST_AUTO_TEST_CASE(TEST_MAT_PREFIX(Resize))
{
	namespace mat = TEST_MAT_NAMESPACE;

	mat::Array<long, mat::Dynamic, mat::Dynamic, 0, 40, 50> dyn4050(4, 5);
	mat::Array<long, mat::Dynamic, mat::Dynamic, 0, 20, 20> dyn2020(2, 2);
	dyn4050.setConstant(4);
	dyn2020.setConstant(9);
	BOOST_CHECK_THROW((dyn4050.block<1, 3>(0,3)), std::runtime_error);
	BOOST_CHECK_THROW(dyn2020 -= dyn4050, std::runtime_error);

	// assignment tests
	BOOST_CHECK_THROW((dyn2020.block<2, 2>(0, 0) = dyn4050), std::runtime_error);
	dyn2020 = dyn4050.block<2, 2>(0, 0);
	BOOST_CHECK(dyn2020.rows() == 2 && dyn2020.cols() == 2);
	dyn2020 = dyn4050;
	BOOST_CHECK(dyn2020.rows() == 4 && dyn2020.cols() == 5);

	// exceeding sizes
	BOOST_CHECK_THROW(dyn2020.resize(80, 1), std::runtime_error);
	BOOST_CHECK_THROW(dyn2020.resize(1, 80), std::runtime_error);

	// operator missmatch and match
	dyn2020.resize(4, 2);
	BOOST_CHECK_THROW(dyn4050 *= dyn2020, std::runtime_error);
	dyn2020.resize(1, 5);
	dyn2020.setConstant(400);
	BOOST_CHECK_THROW(dyn2020 /= dyn4050, std::runtime_error);
	dyn2020.resize(4, 5);
	dyn2020.setConstant(400);
	dyn2020 /= dyn4050;
	const long repl_sum = dyn2020.replicate<3,2>().sum();
	BOOST_CHECK_EQUAL(repl_sum, 100 * (4 * 5) * (3 * 2));

	// partial block operations
	BOOST_CHECK_EQUAL(dyn4050.sum(), 4 * (4 * 5));
	long a123[] = {1, 2, 3};
	mat::Map<mat::Array<long, 1, 3> > a123_map(a123);
	dyn4050.block<1, 2>(0, 3) = a123_map.block<1, 2>(0, 0);
	dyn4050.block(1, 2, 1, 3) = a123_map;
	BOOST_CHECK_EQUAL((dyn4050.block<4, 1>(0,3).sum()), 1 + 2 + 4 + 4);
};

BOOST_AUTO_TEST_SUITE_END()
