#ifndef micromat_h_TEPKbYQAhNrv
#define micromat_h_TEPKbYQAhNrv

#include <stdint.h>
#include <complex>
#include <assert.h>

#ifndef micromat_assert
#define micromat_assert assert
#endif

/**
Implements block and matrix operations. The Eigen API is used if possible.
Limitations: Fixed maximum size, column-major only.
*/

namespace umat
{

typedef int_fast16_t Index;

/// Match the value in Eigen.
const int Dynamic = -1;

template<bool> struct StaticAssert;
template<> struct StaticAssert<true> {};

template<typename T>
struct NumTraits
{};

template<>
struct NumTraits<float>
{
	typedef float Real;
	enum
	{
		IsComplex = 0
	};
};


template<typename T>
struct NumTraits<std::complex<T> >
{
	typedef T Real;
	enum
	{
		IsComplex = 1
	};
};

// combine size with data, because otherwise fixed size dummy
// size objects occupy one byte of RAM.
/// Internal structure
template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime,
		int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
struct DataAndSize
{
	// data first. Alignment might be important
	Scalar data[MaxRowsAtCompileTime * MaxColsAtCompileTime];

	Index rows() const
		{return RowsAtCompileTime;}
	Index cols() const
		{return ColsAtCompileTime;}

	void resize(const Index new_rows, const Index new_cols)
	{
		micromat_assert(new_rows == RowsAtCompileTime);
		micromat_assert(new_cols == ColsAtCompileTime);
	};
};

template<typename Scalar, int RowsAtCompileTime,
		int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
struct DataAndSize<Scalar, RowsAtCompileTime, Dynamic,
		MaxRowsAtCompileTime, MaxColsAtCompileTime>
{
	// data first. Alignment might be important
	Scalar data[MaxRowsAtCompileTime * MaxColsAtCompileTime];
	Index cols_val;

	Index rows() const
		{return RowsAtCompileTime;}
	Index cols() const
		{return cols_val;}
	Index& cols()
		{return cols_val;}

	void resize(const Index new_rows, const Index new_cols)
	{
		micromat_assert(new_rows == RowsAtCompileTime);
		cols_val = new_cols;
	};
};

template<typename Scalar, int ColsAtCompileTime,
		int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
struct DataAndSize<Scalar, Dynamic, ColsAtCompileTime,
		MaxRowsAtCompileTime, MaxColsAtCompileTime>
{
	// data first. Alignment might be important
	Scalar data[MaxRowsAtCompileTime * MaxColsAtCompileTime];
	Index rows_val;

	Index rows() const
		{return rows_val;}
	Index& rows()
		{return rows_val;}
	Index cols() const
		{return ColsAtCompileTime;}

	void resize(const Index new_rows, const Index new_cols)
	{
		rows_val = new_rows;
		micromat_assert(new_cols == ColsAtCompileTime);
	};
};

template<typename Scalar,
		int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
struct DataAndSize<Scalar, Dynamic, Dynamic,
		MaxRowsAtCompileTime, MaxColsAtCompileTime>
{
	// data first. Alignment might be important
	Scalar data[MaxRowsAtCompileTime * MaxColsAtCompileTime];

	Index rows() const
		{return rows_val;}
	Index& rows()
		{return rows_val;}
	Index cols() const
		{return cols_val;}
	Index& cols()
		{return cols_val;}

	Index rows_val;
	Index cols_val;

	void resize(const Index new_rows, const Index new_cols)
	{
		rows_val = new_rows;
		cols_val = new_cols;
	};
};


/// micromat operation functors
/** needed for lazy evaluated expression templates */
namespace umop
{

struct Assign
{
	template<typename T, typename S>
	void operator() (T& left, const S right)
		{left = right;}
};

struct AddAssign
{
	template<typename T, typename S>
	void operator() (T& left, const S right)
		{left += right;}
};

struct SubAssign
{
	template<typename T, typename S>
	void operator() (T& left, const S right)
		{left -= right;}
};

struct MulAssign
{
	template<typename T, typename S>
	void operator() (T& left, const S right)
		{left *= right;}
};

struct DivAssign
{
	template<typename T, typename S>
	void operator() (T& left, const S right)
		{left /= right;}
};

struct AddSubMulDiv
{
	template<typename T, typename S>
	struct Result
		{typedef T type;}; // default compromise. At least the same as with *=

	template<typename T>
	struct Result<std::complex<T>, T>
		{typedef std::complex<T> type;};

	template<typename T>
	struct Result<T, std::complex<T> >
		{typedef std::complex<T> type;};
};

struct Add : public AddSubMulDiv
{
	template<typename T, typename S>
	typename Result<T, S>::type operator() (const T left, const S right) const
		{return left + right;}
};

struct Sub : public AddSubMulDiv
{
	template<typename T, typename S>
	typename Result<T, S>::type operator() (const T left, const S right) const
		{return left - right;}
};

struct Mul : public AddSubMulDiv
{
	template<typename T, typename S>
	typename Result<T, S>::type operator() (const T left, const S right) const
		{return left * right;}
};

struct Div : public AddSubMulDiv
{
	template<typename T, typename S>
	typename Result<T, S>::type operator() (const T left, const S right) const
		{return left / right;}
};

struct Equal
{
	template<typename T, typename S>
	struct Result
		{typedef bool type;};

	template<typename T, typename S>
	bool operator() (const T left, const S right) const
		{return left == right;}
};

struct Abs2
{
	template<typename T>
	struct Result
	{
		typedef T type;

		static type Abs2(const T scalar)
		{
			// Real value specialisation needed, because C++03 has no
			// std::real(float)
			return scalar * scalar;
		}
};

	template<typename T>
	struct Result<std::complex<T> >
	{
		typedef T type;

		static type Abs2(const std::complex<T> cpx)
		{
			return std::real(cpx) * std::real(cpx)
					+ std::imag(cpx) * std::imag(cpx);
		}
	};

	template<typename T>
	typename Result<T>::type operator() (const T input) const
	{
		return Result<T>::Abs2(input);
	}
};

struct Abs
{
	template<typename T>
	struct Result
		{typedef T type;};

	template<typename T>
	struct Result<std::complex<T> >
		{typedef T type;};

	template<typename T>
	typename Result<T>::type operator() (const T input) const
	{
		return std::abs(input);
	}
};

template<typename T>
struct Constant
{
	Constant(const T& val_)
		: val(val_)
	{};

	T operator() (const Index /*row*/, const Index /*col*/) const
	{
		return val;
	};

private:
	const T val;
};

} // namespace op

template<typename Storage>
class View;

template<typename Storage>
class ConstView;

template<typename Op, typename Storage>
class CwiseNullaryOp;

template<typename Op, typename Storage>
class CwiseUnaryOp;

template<typename Op, typename Left, typename Right>
class CwiseBinaryOp;

template<typename Storage>
class Replicate;

template<typename Storage, typename Scalar>
class BlockMixins
{
public:
	template<typename Other>
	Storage& operator = (const Other& other)
		{return CwiseOp<umop::Assign, Other>(other);}

	template<typename Other>
	Storage& operator += (const Other& other)
		{return CwiseOp<umop::AddAssign, Other>(other);}

	template<typename Other>
	Storage& operator -= (const Other& other)
		{return CwiseOp<umop::SubAssign, Other>(other);}

	template<typename Other>
	Storage& operator *= (const Other& other)
		{return CwiseOp<umop::MulAssign, Other>(other);}

	template<typename Other>
	Storage& operator /= (const Other& other)
		{return CwiseOp<umop::DivAssign, Other>(other);}

	template<typename Other>
	CwiseBinaryOp<umop::Add, Storage, Other>
	operator + (const Other& other) const
	{
		return CwiseBinaryOp<umop::Add, Storage, Other>
				(*static_cast<const Storage*>(this), other);
	};

	template<typename Other>
	CwiseBinaryOp<umop::Sub, Storage, Other>
	operator - (const Other& other) const
	{
		return CwiseBinaryOp<umop::Sub, Storage, Other>
				(*static_cast<const Storage*>(this), other);
	};

	template<typename Other>
	CwiseBinaryOp<umop::Mul, Storage, Other>
	operator * (const Other& other) const
	{
		return CwiseBinaryOp<umop::Mul, Storage, Other>
				(*static_cast<const Storage*>(this), other);
	};

	template<typename Other>
	CwiseBinaryOp<umop::Div, Storage, Other>
	operator / (const Other& other) const
	{
		return CwiseBinaryOp<umop::Div, Storage, Other>
				(*static_cast<const Storage*>(this), other);
	};

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
	Constant(const Index c_rows, const Index c_cols, const Scalar val)
	{
		return CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
				(c_rows, c_cols, typename umop::Constant<Scalar>(val));
	}

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
	Constant(const Scalar val)
	{
		return Constant(
			Storage::RowsAtCompileTime,
			Storage::ColsAtCompileTime,
			val
		);
	}

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
	Zero(const Index c_rows, const Index c_cols)
	{
		return Constant(c_rows, c_cols, 0);
	}

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
	Ones(const Index c_rows, const Index c_cols)
	{
		return Constant(c_rows, c_cols, 1);
	}

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
		Zero() {return Constant(0);}

	static CwiseNullaryOp<typename umop::Constant<Scalar>, Storage>
		Ones() {return Constant(1);}

	void setConstant(const Scalar& val)
	{
		Storage& me = *static_cast<Storage*>(this);
		me = Constant(
			me.rows(),
			me.cols(),
			val
		);
	}

	Scalar sum() const
	{
		Scalar accumulate = Scalar();
		const Storage& me = *static_cast<const Storage*>(this);
		const Index rows_end = me.rows();
		const Index cols_end = me.cols();
		for(Index j = 0; j < cols_end; ++j)
		{
			for(Index i = 0; i < rows_end; ++i)
			{
				accumulate += me.coeff(i, j);
			}
		}
		return accumulate;
	}

	void PrintDebug()
	{
		/*Storage& me = *static_cast<Storage*>(this);
		const Index rows_end = me.rows();
		const Index cols_end = me.cols();
		for(Index j = 0; j < rows_end; ++j)
		{
			for(Index i = 0; i < cols_end; ++i)
			{
				std::cout << me(i, j) << "\t";
			}
			std::cout << "\n";
		}*/
	}

	CwiseUnaryOp<umop::Abs, Storage> abs() const
	{
		return CwiseUnaryOp<umop::Abs, Storage>(
			*static_cast<const Storage*>(this)
		);
	}

	CwiseUnaryOp<umop::Abs2, Storage> abs2() const
	{
		return CwiseUnaryOp<umop::Abs2, Storage>(
			*static_cast<const Storage*>(this)
		);
	}

	/// Block of size (p,q), starting at (i,j)
	View<Storage> block(Index i, Index j, Index p, Index q)
	{
		Storage& me = *static_cast<Storage*>(this);
		AssertBounds(i, j);
		micromat_assert(p + i <= me.rows() && q + j <= me.cols());
		return View<Storage>(me, i, j, p, q);
	}

	/// Block of size (p,q), starting at (i,j)
	template<int p, int q>
	View<Storage>
	block(Index i, Index j)
	{
		return block(i, j, p, q);
	}

	/// Block of size (p,q), starting at (i,j)
	ConstView<Storage> block(Index i, Index j, Index p, Index q) const
	{
		const Storage& me = *static_cast<const Storage*>(this);
		AssertBounds(i, j);
		micromat_assert(p + i <= me.rows() && q + j <= me.cols());
		return ConstView<Storage>(me, i, j, p, q);
	}

	/// Block of size (p,q), starting at (i,j)
	template<int p, int q>
	ConstView<Storage>
	block(Index i, Index j) const
	{
		return block(i, j, p, q);
	}

	template<int row_factor, int col_factor>
	Replicate<Storage>
	replicate()
	{
		return Replicate<Storage>(
			*static_cast<Storage*>(this),
			row_factor, col_factor
		);
	}

	Replicate<Storage>
	replicate(Index row_factor, Index col_factor)
	{
		return replicate(row_factor, col_factor);
	}

protected:

	void AssertBounds(const Index row, const Index col) const
	{
		micromat_assert(row >= 0 && col >= 0);
		const Storage& me = *static_cast<const Storage*>(this);
		micromat_assert(me.rows() > row);
		micromat_assert(me.cols() > col);
	}

	template<typename Operation, typename Other>
	Storage& CwiseOp(const Other& other)
	{
		StaticAssert<(
			Index(Storage::RowsAtCompileTime) == Dynamic ||
			Index(Other::RowsAtCompileTime) == Dynamic ||
			Index(Other::RowsAtCompileTime) ==
					Index(Storage::RowsAtCompileTime)
		)> check_rows_missmatch;
		(void)check_rows_missmatch;
		StaticAssert<(
			Index(Storage::ColsAtCompileTime) == Dynamic ||
			Index(Other::ColsAtCompileTime) == Dynamic ||
			Index(Other::ColsAtCompileTime) ==
					Index(Storage::ColsAtCompileTime)
		)> check_cols_missmatch;
		(void)check_cols_missmatch;
		Storage& me = *static_cast<Storage*>(this);
		const Index rows_end = me.rows();
		const Index cols_end = me.cols();
		micromat_assert(rows_end == other.rows());
		micromat_assert(cols_end == other.cols());
		Operation op;
		for(Index j = 0; j < cols_end; ++j)
		{
			for(Index i = 0; i < rows_end; ++i)
			{
				op(
					me.coeffRef(i, j),
					other.coeff(i, j)
				);
			}
		}
		return me;
	}
};

template<typename Storage>
class View : public BlockMixins<View<Storage>,
		typename Storage::Scalar>
{
public:
	typedef View<Storage> OwnType;
	typedef BlockMixins<OwnType, typename Storage::Scalar > BaseType;
	typedef typename Storage::Scalar Scalar;
	enum
	{
		RowsAtCompileTime = Dynamic,
		ColsAtCompileTime = Dynamic
	};

	View(Storage& base_, Index row_offset_, Index col_offset_, Index row_count_, Index col_count_)
		: base(base_)
		, row_offset(row_offset_)
		, col_offset(col_offset_)
		, row_count(row_count_)
		, col_count(col_count_)
	{};

	typename Storage::Scalar& operator() (Index row, Index col)
	{
		BaseType::AssertBounds(row, col);
		return base(row + row_offset, col + col_offset);
	}

	using BaseType::operator= ;

	typename Storage::Scalar& coeffRef(Index row, Index col)
		{return base.coeffRef(row + row_offset, col + col_offset);}

	typename Storage::Scalar coeff(Index row, Index col) const
		{return base.coeff(row + row_offset, col + col_offset);}

	Index rows() const
		{return row_count;};
	Index cols() const
		{return col_count;};

private:
	Storage& base;
	const Index row_offset;
	const Index col_offset;
	const Index row_count;
	const Index col_count;
};

template<typename Storage>
class ConstView : public BlockMixins<ConstView<Storage>,
		typename Storage::Scalar>
{
public:
	typedef ConstView<Storage> OwnType;
	typedef BlockMixins<OwnType, typename Storage::Scalar> BaseType;
	typedef typename Storage::Scalar Scalar;
	enum
	{
		RowsAtCompileTime = Dynamic,
		ColsAtCompileTime = Dynamic
	};

	ConstView(const Storage& base_, Index row_offset_, Index col_offset_, Index row_count_, Index col_count_)
		: base(base_)
		, row_offset(row_offset_)
		, col_offset(col_offset_)
		, row_count(row_count_)
		, col_count(col_count_)
	{};

	const typename Storage::Scalar& operator() (Index row, Index col)
	{
		BaseType::AssertBounds(row, col);
		return base(row + row_offset, col + col_offset);
	}

	const typename Storage::Scalar coeff(Index row, Index col) const
		{return base.coeff(row + row_offset, col + col_offset);}

	Index rows() const
		{return row_count;};
	Index cols() const
		{return col_count;};

private:
	const Storage& base;
	const Index row_offset;
	const Index col_offset;
	const Index row_count;
	const Index col_count;
};

template<typename Op, typename Storage>
class CwiseNullaryOp : public BlockMixins<CwiseNullaryOp<Op, Storage>,
	typename Storage::Scalar>
{
public:
	typedef CwiseNullaryOp<Op, Storage> OwnType;
	typedef typename Storage::Scalar Scalar;
	typedef BlockMixins<OwnType, Scalar> BaseType;

	enum
	{
		RowsAtCompileTime = Storage::RowsAtCompileTime,
		ColsAtCompileTime = Storage::ColsAtCompileTime
	};

	CwiseNullaryOp(const Index rows_, const Index cols_, const Op& op_)
		: op(op_)
		, rows_val(rows_)
		, cols_val(cols_)
	{};

	Scalar operator() (Index row, Index col) const
	{
		BaseType::AssertSize(row, col);
		return coeff(row, col);
	}

	Scalar coeff(Index row, Index col) const
		{return op(row, col);}

	Index rows() const
		{return rows_val;};
	Index cols() const
		{return cols_val;};

private:
	const Op op;
	const Index rows_val;
	const Index cols_val;
};

template<typename Op, typename Storage>
class CwiseUnaryOp : public BlockMixins<CwiseUnaryOp<Op, Storage>,
	typename Op::template Result<typename Storage::Scalar>::type>
{
public:
	typedef CwiseUnaryOp<Op, Storage> OwnType;
	typedef typename Op::template Result<typename Storage::Scalar>::type
			Scalar;
	typedef BlockMixins<OwnType, Scalar> BaseType;

	enum
	{
		RowsAtCompileTime = Storage::RowsAtCompileTime,
		ColsAtCompileTime = Storage::ColsAtCompileTime
	};

	CwiseUnaryOp(const Storage& base_)
		: base(base_)
		, op()
	{};


	Scalar operator() (Index row, Index col) const
	{
		BaseType::AssertBounds(row, col);
		return coeff(row, col);
	}

	Scalar coeff(Index row, Index col) const
		{return op(base.coeff(row, col));}

	Index rows() const
		{return base.rows();};
	Index cols() const
		{return base.cols();};

private:
	const Storage& base;
	const Op op;
};


template<typename Op, typename Left, typename Right>
class CwiseBinaryOp : public BlockMixins<CwiseBinaryOp<Op, Left, Right>,
	typename Op::template Result<typename Left::Scalar,
	typename Right::Scalar>::type>
{
public:
	typedef CwiseBinaryOp<Op, Left, Right> OwnType;
	typedef typename Op::template Result<typename Left::Scalar,
			typename Right::Scalar>::type Scalar;
	typedef BlockMixins<OwnType, Scalar> BaseType;

	enum
	{
		RowsAtCompileTime = Left::RowsAtCompileTime,
		ColsAtCompileTime = Left::ColsAtCompileTime
	};

	CwiseBinaryOp(const Left& left_, const Right& right_)
		: left(left_)
		, right(right_)
		, op()
	{
		StaticAssert<(
			Index(Left::RowsAtCompileTime) == Dynamic ||
			Index(Right::RowsAtCompileTime) == Dynamic ||
			Index(Left::RowsAtCompileTime) ==
					Index(Right::RowsAtCompileTime)
		)> check_rows_missmatch;
		(void)check_rows_missmatch;
		StaticAssert<(
			Index(Left::ColsAtCompileTime) == Dynamic ||
			Index(Right::ColsAtCompileTime) == Dynamic ||
			Index(Left::ColsAtCompileTime) ==
					Index(Right::ColsAtCompileTime)
		)> check_cols_missmatch;
		(void)check_cols_missmatch;
		micromat_assert(left.rows() == right.rows());
		micromat_assert(left.cols() == right.cols());
	};


	Scalar operator() (Index row, Index col) const
	{
		Left::AssertBounds(row, col);
		return coeff(row, col);
	}

	Scalar coeff(Index row, Index col) const
		{return op(left.coeff(row, col), right.coeff(row, col));}

	Index rows() const
		{return left.rows();};
	Index cols() const
		{return left.cols();};

private:
	const Left& left;
	const Right& right;
	const Op op;
};

template<typename Scalar_, int rows_, int cols_,
		int options = 0, int max_rows_ = rows_, int max_cols_ = cols_>
class Array : public BlockMixins<Array<Scalar_, rows_, cols_,
		options, max_rows_, max_cols_>, Scalar_ >
{
public:
	typedef Array<Scalar_, rows_, cols_, options, max_rows_, max_cols_> OwnType;
	typedef BlockMixins<OwnType, Scalar_ > BaseType;
	typedef Scalar_ Scalar;
	enum
	{
		RowsAtCompileTime = rows_,
		ColsAtCompileTime = cols_,
		MaxRowsAtCompileTime = max_rows_,
		MaxColsAtCompileTime = max_cols_
	};

	template<typename Other>
	OwnType& operator = (const Other& other)
	{
		resize(other.rows(), other.cols());
		return BaseType::operator=(other);
	};

	Array()
	{
		StaticSizeCheck();
	};

	Array(const Index new_rows, const Index new_cols)
	{
		StaticSizeCheck();
		resize(new_rows, new_cols);
	};

	template<typename Other>
	Array(const Other& other)
	{
		StaticSizeCheck();
		resize(other.rows(), other.cols());
		BaseType::operator=(other);
	}

	Scalar& operator() (Index row, Index col)
	{
		BaseType::AssertBounds(row, col);
		return coeffRef(row, col);
	}

	const Scalar operator() (Index row, Index col) const
	{
		BaseType::AssertBounds(row, col);
		return coeff(row, col);
	}

	Scalar& coeffRef(Index row, Index col)
		{return data_and_size.data[row + col * MaxRowsAtCompileTime];}

	Scalar coeff(Index row, Index col) const
		{return data_and_size.data[row + col * MaxRowsAtCompileTime];}

	Scalar* data()
		{return data_and_size.data;}

	const Scalar* data() const
		{return data_and_size.data;}

	Index rows() const
		{return data_and_size.rows();}
	Index cols() const
		{return data_and_size.cols();}

	void resize(const Index new_rows, const Index new_cols)
	{
		micromat_assert(new_rows >= 0 && new_rows <= MaxRowsAtCompileTime);
		micromat_assert(new_cols >= 0 && new_cols <= MaxColsAtCompileTime);
		data_and_size.resize(new_rows, new_cols);
	};

private:
	void StaticSizeCheck()
	{
		StaticAssert<(
			Index(RowsAtCompileTime) == Dynamic ||
			Index(RowsAtCompileTime) == MaxRowsAtCompileTime
		)> check_rows_size_error;
		StaticAssert<(
			Index(ColsAtCompileTime) == Dynamic ||
			Index(ColsAtCompileTime) == MaxColsAtCompileTime
		)> check_cols_size_error;
		StaticAssert<(MaxRowsAtCompileTime >= 0)> check_negative_max_rows;
		StaticAssert<(MaxColsAtCompileTime >= 0)> check_negative_max_cols;
		(void)check_cols_size_error;
		(void)check_rows_size_error;
		(void)check_negative_max_rows;
		(void)check_negative_max_cols;
	};

	DataAndSize<Scalar, RowsAtCompileTime, ColsAtCompileTime,
			MaxRowsAtCompileTime, MaxColsAtCompileTime> data_and_size;
};

/// Mimics semantics of PlainObjectType, but with external storage
template<typename PlainObjectType>
class Map : public BlockMixins<Map<PlainObjectType>,
		typename PlainObjectType::Scalar>
{
public:
	typedef Map<PlainObjectType> OwnType;
	typedef BlockMixins<Map<PlainObjectType>, typename PlainObjectType::Scalar>
			BaseType;
	typedef typename PlainObjectType::Scalar Scalar;
	enum
	{
		RowsAtCompileTime = PlainObjectType::RowsAtCompileTime,
		ColsAtCompileTime = PlainObjectType::ColsAtCompileTime
	};

	using BaseType::operator= ;

	Map(Scalar* data_)
		: data(data_)
	{}

	Scalar& operator() (Index row, Index col)
	{
		micromat_assert(
			row >= 0 && row < RowsAtCompileTime &&
			col >= 0 && col < ColsAtCompileTime
		);
		return data[row + col * RowsAtCompileTime];
	}

	Scalar& coeffRef(Index row, Index col)
		{return data[row + col * RowsAtCompileTime];}

	Scalar coeff(Index row, Index col) const
		{return data[row + col * RowsAtCompileTime];}

	Index rows() const
		{return RowsAtCompileTime;};
	Index cols() const
		{return ColsAtCompileTime;};

protected:
	Scalar* data;
};

/// Expression of the multiple replication of a matrix or vector.
template<typename Storage>
class Replicate
		: public BlockMixins<Replicate<Storage>,typename Storage::Scalar>
{
public:

	typedef Replicate<Storage> OwnType;
	typedef BlockMixins<OwnType, typename Storage::Scalar> BaseType;
	typedef typename Storage::Scalar Scalar;
	enum
	{
		RowsAtCompileTime = Dynamic,
		ColsAtCompileTime = Dynamic
	};

	Replicate(Storage& base_, Index row_factor_, Index col_factor_)
		: base(base_)
		, row_factor(row_factor_)
		, col_factor(col_factor_)
	{};

	Scalar& operator() (Index row, Index col)
	{
		BaseType::AssertBounds(row, col);
		return coeffRef(row, col);
	}

	Scalar operator() (Index row, Index col) const
	{
		BaseType::AssertBounds(row, col);
		return coeff(row, col);
	}

	Scalar& coeffRef(Index row, Index col)
	{
		const Index wrapped_row = row % base.rows();
		const Index wrapped_col = col % base.cols();
		return base.coeffRef(wrapped_row, wrapped_col);
	}

	Scalar coeff(Index row, Index col) const
	{
		const Index wrapped_row = row % base.rows();
		const Index wrapped_col = col % base.cols();
		return base.coeff(wrapped_row, wrapped_col);
	}

	Index rows() const
		{return row_factor * base.rows();};
	Index cols() const
		{return col_factor * base.cols();};

	using BaseType::operator= ;

protected:
	Storage& base;
	const Index row_factor;
	const Index col_factor;
};

} // namespace umat

#endif // micromat_h_TEPKbYQAhNrv
