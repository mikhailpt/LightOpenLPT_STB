#pragma once

#include <iterator>
#include "math\tnt.h"
#include "math\jama_lu.h"

namespace Kriging
{
	/*First, Last - iterators to known elements in container
	FirstEstimate, LastEstimate - iterators to unknown elements in container
	Dist - function to calc distance between two elements e.g double(const element&, const element&)
	DataAccess - function to access element value e.g float&(element&) */
	template<class _It, class _ItEst, class _Dist, class _DataAccess>
	void OrdinaryKriging(_It First, _It Last, _ItEst FirstEstimate, _ItEst LastEstimate, _Dist Dist, _DataAccess DataAccess)
	{
		static_assert(std::is_reference<decltype(DataAccess(*FirstEstimate))>::value, "DataAccess must return lvalue");
		using DistType = std::remove_reference<decltype(Dist(*First, *First))>::type;
		using ValueType = std::remove_reference<decltype(DataAccess(*FirstEstimate))>::type;

		size_t size = std::distance(First, Last);
		Array2D<DistType> A(size + 1, size + 1);	//distance matrix
		for (size_t i = 0; i < size; ++i)
		{
			auto& ei = *(First + i);
			for (size_t j = 0; j <= i; ++j)	//use symmetry
			{
				auto& ej = *(First + j);
				A[i][j] = A[j][i] = Dist(ei, ej);
			}

			A[i][size] = A[size][i] = DistType(1);	//last row and column fill with 1
		}
		A[size][size] = DistType(0);	//last diagonal element = 0
		A = __matinvert(A);

		Array2D<ValueType> v(1, size + 1);	//vector of known values
		for (size_t j = 0; j < size; ++j)
		{
			auto ej = *(First + j);	//make a copy of element here to not conflict with const_iterator
			v[0][j] = DataAccess(ej);
		}
		ValueType zero;
		zero -= zero;
		v[0][size] = zero;
		auto w = __matmult(v, A);

		for (auto& it = FirstEstimate; it != LastEstimate; ++it)
		{
			Array2D<DistType> b(size + 1, 1);
			for (size_t j = 0; j < size; ++j)
			{
				auto ej = *(First + j);
				b[j][0] = Dist(*it, ej);
			}
			b[size][0] = DistType(1);

			auto est = __matmult(w, b)[0][0];	//est = v * A^-1 * b
			DataAccess(*it) = (ValueType)est;
		}
	}

	template<class TypeA, class TypeB>
	auto __matmult(const Array2D<TypeA>& A, const Array2D<TypeB>& B)
	{
		using TypeC = decltype(A[0][0] * B[0][0] + A[0][0] * B[0][0]);
		if (A.dim2() != B.dim1())
			throw std::exception("Can't multiply matrix");

		TypeC zero;
		zero -= zero;
		int M = A.dim1();
		int N = A.dim2();
		int K = B.dim2();
				
		Array2D<TypeC> C(M, K);
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < K; j++)
			{
				TypeC sum = zero;
				for (int k = 0; k < N; k++)
				{
					sum += A[i][k] * B[k][j];
				}
				C[i][j] = sum;
			}
		}
		return C;
	}

	template <class Type>
	Array2D<Type> __matinvert(const Array2D<Type>& A)
	{
		if (A.dim2() != A.dim1())
			throw std::exception("Can't invert non square matrix");

		Array2D<Type> E(A.dim1(), A.dim1(), Type(0));
		for (int i = 0; i < A.dim1(); ++i)
		{
			E[i][i] = Type(1);
		}
		JAMA::LU<Type> S(A);
		if (!S.isNonsingular())
			throw std::exception("Can't invert matrix with determinant=0");
		return S.solve(E);
	}
}
