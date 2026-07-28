#ifndef ENSEMBLEVIEW_H_
#define ENSEMBLEVIEW_H_

#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "Ensemble.h"

/**
 * @brief Thrown when an EnsembleView is used after the ensemble it borrows from has been
 *        modified, replaced or destroyed.
 */
class EnsembleViewInvalidated : public std::runtime_error
{
public:
	EnsembleViewInvalidated() : std::runtime_error(
		"EnsembleView used after the ensemble was modified, replaced, or destroyed") {}
};

/**
 * @brief A borrowed, column-major window into an Ensemble's numeric storage.
 *
 * The point is to let a caller look at (and edit in place) a live ensemble without copying
 * the matrix - for a large ensemble that copy is the expensive part. `data()` is the
 * ensemble's own buffer, so reads see current values and writes land directly in it.
 *
 * A view is only meaningful while the ensemble's storage stays put. It goes invalid when:
 *
 *   - the ensemble is **destroyed** or **assigned over** (the guard token expires or is
 *     replaced - see EnsembleViewGuard);
 *   - `reals` is **reallocated** by any mutator - drop_rows, keep_rows, reorder,
 *     append_other_rows, resize, ... - which is detected by re-checking the buffer address
 *     and the dimensions rather than by bookkeeping at each of those ~28 call sites;
 *   - the whole matrix is swapped for a **same-shaped** one by set_eigen()/from_eigen_mat().
 *     Eigen reuses the buffer in that case, so the address/dimension check cannot see it and
 *     those two call invalidate() explicitly.
 *
 * Every accessor re-checks first and throws EnsembleViewInvalidated, so a stale view raises
 * cleanly instead of reading freed memory. In-place edits that do not resize (writing
 * through `data()`, `replace_col`, `update_real_ip`) keep a view valid, which is exactly the
 * window in which editing is useful.
 *
 * Not copyable: two handles onto one borrowed buffer is a mistake worth making loud.
 *
 * Threading: none. The tools are single-threaded at the iteration boundaries where a view
 * is meant to be taken.
 */
class EnsembleView
{
public:
	explicit EnsembleView(Ensemble& ens)
		: ens_(&ens), guard_(ens.get_view_guard()),
		  base_(ens.get_eigen_ptr()->data()),
		  rows_(ens.get_eigen_ptr()->rows()),
		  cols_(ens.get_eigen_ptr()->cols()) {}

	EnsembleView(const EnsembleView&) = delete;
	EnsembleView& operator=(const EnsembleView&) = delete;

	/// Is the borrowed buffer still the ensemble's current storage?
	bool valid() const
	{
		// order matters: if the ensemble is gone the token has expired, and we must not
		// touch ens_ at all
		if (guard_.expired())
			return false;
		const Eigen::MatrixXd* m = ens_->get_eigen_ptr();
		return (m->data() == base_) && (m->rows() == rows_) && (m->cols() == cols_);
	}

	void check() const { if (!valid()) throw EnsembleViewInvalidated(); }

	/// The ensemble's own buffer. Writes through this pointer are seen by the algorithm.
	double* data() const { check(); return const_cast<double*>(base_); }

	int64_t rows() const { return rows_; }
	int64_t cols() const { return cols_; }

	/// Eigen::MatrixXd is COLUMN-major: element (i,j) lives at data()[i + j*rows()].
	std::array<int64_t, 2> strides_bytes() const
	{
		return { (int64_t)sizeof(double), (int64_t)sizeof(double) * rows_ };
	}

	/// Read one element, bounds- and validity-checked.
	double at(int64_t i, int64_t j) const
	{
		check();
		if ((i < 0) || (i >= rows_) || (j < 0) || (j >= cols_))
			throw std::out_of_range("EnsembleView::at() index out of range");
		return base_[i + (j * rows_)];
	}

	/// Write one element in place, bounds- and validity-checked.
	void set(int64_t i, int64_t j, double value) const
	{
		check();
		if ((i < 0) || (i >= rows_) || (j < 0) || (j >= cols_))
			throw std::out_of_range("EnsembleView::set() index out of range");
		const_cast<double*>(base_)[i + (j * rows_)] = value;
	}

	// Labels are read from the ensemble rather than cached, so a rename cannot leave the
	// view reporting stale names alongside live numbers.
	std::vector<std::string> row_names() const { check(); return ens_->get_real_names(); }
	std::vector<std::string> col_names() const { check(); return ens_->get_var_names(); }

private:
	Ensemble* ens_;
	std::weak_ptr<int> guard_;
	const double* base_;
	int64_t rows_, cols_;
};

#endif // ENSEMBLEVIEW_H_
