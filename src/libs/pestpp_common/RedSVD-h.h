/*
 * A header-only version of RedSVD
 *
 * Copyright (c) 2014 Nicolas Tessore
 *
 * based on RedSVD
 *
 * Copyright (c) 2010 Daisuke Okanohara
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above Copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above Copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the authors nor the names of its contributors
 *    may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 */

#ifndef REDSVD_MODULE_H
#define REDSVD_MODULE_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cstdlib>
#include <cmath>
#include <random>

namespace RedSVD
{
	template<typename Scalar>
	inline void sample_gaussian_old(Scalar& x, Scalar& y)
	{
		using std::sqrt;
		using std::log;
		using std::cos;
		using std::sin;
		
		const Scalar PI(3.1415926535897932384626433832795028841971693993751);
		Scalar v1 = (Scalar)(std::rand() + Scalar(1)) / ((Scalar)RAND_MAX+Scalar(2));
		Scalar v2 = (Scalar)(std::rand() + Scalar(1)) / ((Scalar)RAND_MAX+Scalar(2));
		/*Scalar v1 = (Scalar)(rand_gen() + Scalar(1)) / ((Scalar)rand_gen.max() + Scalar(2));
		Scalar v1 = (Scalar)(rand_gen() + Scalar(1)) / ((Scalar)rand_gen.max() + Scalar(2));*/

		
		
		Scalar len = sqrt(Scalar(-2) * log(v1));
		x = len * cos(Scalar(2) * PI * v2);
		y = len * sin(Scalar(2) * PI * v2);
	}

	template<typename MatrixType>
	inline void sample_gaussian_old(MatrixType& mat)
	{
		typedef typename MatrixType::Index Index;

		for(Index i = 0; i < mat.rows(); ++i)
		{
			for(Index j = 0; j+1 < mat.cols(); j += 2)
				sample_gaussian(mat(i, j), mat(i, j+1));
			if(mat.cols() % 2)
				sample_gaussian(mat(i, mat.cols()-1), mat(i, mat.cols()-1));
		}
	}

	// std::normal_distribution is implementation-defined, and a Box-Muller transform relies on
	// std::log/std::sin, which are not correctly-rounded and differ by ~1 ULP across libm
	// implementations (glibc / Apple libm / MSVC CRT) - so the randomized projection, and therefore
	// every RedSVD/RedSymEigen result, diverges across platforms even from an identical engine state.
	// Draw standard normals instead via Acklam's rational approximation to the inverse normal CDF
	// over the raw (fully portable) mt19937 integer stream, using only + - * / - no transcendentals -
	// so the same bits are produced on every platform. The uniform is clamped to Acklam's central
	// rational region (~+/-1.96 sigma), which keeps the draw entirely transcendental-free (only the
	// tail branch would need a log). The downstream range finder orthonormalizes the sketch, so the
	// bounded support and any overall scale are immaterial to the recovered subspace.
	//
	// NOTE: a Rademacher (+/-1) projection was tried for full bit-level determinism but it is too
	// crude a sketch for a low-rank range finder; this inverse-CDF draw is a smooth, near-Gaussian,
	// isotropic sub-Gaussian sketch instead. The few genuinely platform-sensitive tests are still
	// pinned to the deterministic SVD_EIGEN package.
	template<typename MatrixType>
	inline void sample_gaussian(MatrixType& mat)
	{
		typedef typename MatrixType::Index Index;
		std::mt19937 generator;
		const double span = static_cast<double>(generator.max()) - static_cast<double>(generator.min());
		// Acklam (2003) inverse normal CDF, central-region rational coefficients
		const double a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02, a3 = -2.759285104469687e+02,
		             a4 = 1.383577518672690e+02, a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00;
		const double b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02, b3 = -1.556989798598866e+02,
		             b4 = 6.680131188771972e+01, b5 = -1.328068155288572e+01;
		const double plow = 0.02425, phigh = 1.0 - plow;  // central rational validity window
		for (Index i = 0; i < mat.rows(); ++i)
		{
			for (Index j = 0; j < mat.cols(); ++j)
			{
				// uniform in (0,1) from the portable integer stream, clamped to the central region
				double p = (static_cast<double>(generator() - generator.min()) + 1.0) / (span + 2.0);
				if (p < plow) p = plow;
				else if (p > phigh) p = phigh;
				double q = p - 0.5, r = q * q;
				mat(i, j) = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6) * q /
				            (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r + 1.0);
			}
		}
	}


	template<typename MatrixType>
	inline void gram_schmidt(MatrixType& mat)
	{
		typedef typename MatrixType::Scalar Scalar;
		typedef typename MatrixType::Index Index;

		static const Scalar EPS(1E-4);

		for(Index i = 0; i < mat.cols(); ++i)
		{
			for(Index j = 0; j < i; ++j)
			{
				Scalar r = mat.col(i).dot(mat.col(j));
				mat.col(i) -= r * mat.col(j);
			}

			Scalar norm = mat.col(i).norm();

			if(norm < EPS)
			{
				for(Index k = i; k < mat.cols(); ++k)
					mat.col(k).setZero();
				return;
			}
			mat.col(i) /= norm;
		}
	}

	template<typename _MatrixType>
	class RedSVD
	{
	public:
		typedef _MatrixType MatrixType;
		typedef typename MatrixType::Scalar Scalar;
		typedef typename MatrixType::Index Index;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

		RedSVD() {}

		RedSVD(const MatrixType& A)
		{
			int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
			compute(A, r);
		}

		RedSVD(const MatrixType& A, const Index rank)
		{
			compute(A, rank);
		}

		void compute(const MatrixType& A, const Index rank)
		{
			if(A.cols() == 0 || A.rows() == 0)
				return;

			Index r = (rank < A.cols()) ? rank : A.cols();

			r = (r < A.rows()) ? r : A.rows();

			// oversampling: build the random sketch with l = r + p columns (capped at the matrix
			// dimension) so the dominant rank-r subspace is captured accurately, then truncate the
			// factorization back to r at the end. this widens the two A-products (A^T*O and A*Y) by p
			// columns but adds no extra A-passes; when the caller requests full rank (r already == the
			// min dimension) it is a no-op.
			const Index oversample = 10;
			Index maxl = (A.rows() < A.cols()) ? A.rows() : A.cols();
			Index l = ((r + oversample) < maxl) ? (r + oversample) : maxl;

			// Gaussian Random Matrix for A^T
			DenseMatrix O(A.rows(), l);
			sample_gaussian(O);

			// Compute Sample Matrix of A^T
			DenseMatrix Y = A.transpose() * O;

			// Orthonormalize Y
			gram_schmidt(Y);

			// Range(B) = Range(A^T)
			DenseMatrix B = A * Y;

			// Gaussian Random Matrix
			DenseMatrix P(B.cols(), l);
			sample_gaussian(P);

			// Compute Sample Matrix of B
			DenseMatrix Z = B * P;

			// Orthonormalize Z
			gram_schmidt(Z);

			// Range(C) = Range(B)
			DenseMatrix C = Z.transpose() * B;

			Eigen::JacobiSVD<DenseMatrix> svdOfC(C, Eigen::ComputeThinU | Eigen::ComputeThinV);

			// C = USV^T  ->  A = Z * U * S * V^T * Y^T ; keep only the leading r triplets (truncate
			// off the p oversampled directions)
			m_matrixU = Z * svdOfC.matrixU().leftCols(r);
			m_vectorS = svdOfC.singularValues().head(r);
			m_matrixV = Y * svdOfC.matrixV().leftCols(r);
		}

		DenseMatrix matrixU() const
		{
			return m_matrixU;
		}

		ScalarVector singularValues() const
		{
			return m_vectorS;
		}

		DenseMatrix matrixV() const
		{
			return m_matrixV;
		}

	private:
		DenseMatrix m_matrixU;
		ScalarVector m_vectorS;
		DenseMatrix m_matrixV;
	};

	template<typename _MatrixType>
	class RedSymEigen
	{
	public:
		typedef _MatrixType MatrixType;
		typedef typename MatrixType::Scalar Scalar;
		typedef typename MatrixType::Index Index;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

		RedSymEigen() {}

		RedSymEigen(const MatrixType& A)
		{
			int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
			compute(A, r);
		}

		RedSymEigen(const MatrixType& A, const Index rank)
		{
			compute(A, rank);
		}

		void compute(const MatrixType& A, const Index rank)
		{
			if(A.cols() == 0 || A.rows() == 0)
				return;

			Index r = (rank < A.cols()) ? rank : A.cols();

			r = (r < A.rows()) ? r : A.rows();

			// Gaussian Random Matrix
			DenseMatrix O(A.rows(), r);
			sample_gaussian(O);

			// Compute Sample Matrix of A
			DenseMatrix Y = A.transpose() * O;

			// Orthonormalize Y
			gram_schmidt(Y);

			DenseMatrix B = Y.transpose() * A * Y;
			Eigen::SelfAdjointEigenSolver<DenseMatrix> eigenOfB(B);

			m_eigenvalues = eigenOfB.eigenvalues();
			m_eigenvectors = Y * eigenOfB.eigenvectors();
		}

		ScalarVector eigenvalues() const
		{
			return m_eigenvalues;
		}

		DenseMatrix eigenvectors() const
		{
			return m_eigenvectors;
		}

	private:
		ScalarVector m_eigenvalues;
		DenseMatrix m_eigenvectors;
	};

	template<typename _MatrixType>
	class RedPCA
	{
	public:
		typedef _MatrixType MatrixType;
		typedef typename MatrixType::Scalar Scalar;
		typedef typename MatrixType::Index Index;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

		RedPCA() {}

		RedPCA(const MatrixType& A)
		{
			int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
			compute(A, r);
		}

		RedPCA(const MatrixType& A, const Index rank)
		{
			compute(A, rank);
		}

		void compute(const DenseMatrix& A, const Index rank)
		{
			RedSVD<MatrixType> redsvd(A, rank);

			ScalarVector S = redsvd.singularValues();

			m_components = redsvd.matrixV();
			m_scores = redsvd.matrixU() * S.asDiagonal();
		}

		DenseMatrix components() const
		{
			return m_components;
		}

		DenseMatrix scores() const
		{
			return m_scores;
		}

	private:
		DenseMatrix m_components;
		DenseMatrix m_scores;
	};
}

#endif
