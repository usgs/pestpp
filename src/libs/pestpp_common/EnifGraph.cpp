#include <set>
#include <thread>
#include <atomic>
#include <map>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "EnifGraph.h"
#include "utilities.h"

using namespace std;


void EnifGraph::from_file(const string& filename, const vector<string>& par_names,
	ofstream& frec, const string& order)
{
	stringstream ss;
	Mat m;
	//Mat::from_file dispatches on the extension: jcb/jco binary (the extended
	//coo variant is genuine sparse triplet storage), mat/cov ascii, csv
	m.from_file(filename);

	vector<string> mrow = m.get_row_names();
	set<string> have(mrow.begin(), mrow.end());
	vector<string> missing;
	for (const auto& n : par_names)
		if (have.find(n) == have.end())
			missing.push_back(n);
	if (!missing.empty())
	{
		ss << "EnifGraph::from_file(): graph file '" << filename << "' is missing "
			<< missing.size() << " adjustable parameters, e.g. '" << missing[0] << "'";
		frec << ss.str() << endl;
		frec << "   a partially covered graph would silently drop parameters from the "
			<< "prior precision, so this is treated as an error" << endl;
		throw runtime_error(ss.str());
	}

	//reorder to the solve ordering
	Mat sub = m.get(par_names, par_names);
	Eigen::SparseMatrix<double> a = sub.get_matrix();

	//structural pattern only: any non-zero is an edge, the value is ignored.
	//symmetrise (G union G^T) and force the diagonal on, so a thresholded
	//correlation matrix, a localiser product or a hand-built adjacency all mean
	//the same thing.
	int p = (int)par_names.size();
	vector<Eigen::Triplet<double>> trips;
	trips.reserve(a.nonZeros() * 2 + p);
	for (int k = 0; k < a.outerSize(); k++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(a, k); it; ++it)
		{
			if (it.value() == 0.0)
				continue;
			trips.push_back(Eigen::Triplet<double>((int)it.row(), (int)it.col(), 1.0));
			trips.push_back(Eigen::Triplet<double>((int)it.col(), (int)it.row(), 1.0));
		}
	for (int i = 0; i < p; i++)
		trips.push_back(Eigen::Triplet<double>(i, i, 1.0));

	adj.resize(p, p);
	//duplicate entries are summed by setFromTriplets, so clamp back to a pattern
	adj.setFromTriplets(trips.begin(), trips.end());
	for (int k = 0; k < adj.outerSize(); k++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(adj, k); it; ++it)
			it.valueRef() = 1.0;
	adj.makeCompressed();

	names = par_names;
	nnz_offdiag = (int)adj.nonZeros() - p;

	//--- ordering and symbolic factorisation ---------------------------------
	//the precision is estimated through a cholesky-like factor, so what each
	//node is regressed on is the support of that FACTOR, not the graph.  the
	//factor carries fill: eliminating a node ties its remaining neighbours to
	//each other, so L has entries where the precision has none.  how much fill
	//depends entirely on the elimination order - natural order on a grid is
	//close to the worst case - so the order is chosen to reduce it and the
	//pattern is found by a symbolic factorisation.  this is what graphite-maps
	//does (metis + a cholmod symbolic factorisation); amd is eigen's built-in
	//equivalent.
	order_method = order;
	for (auto& c : order_method)
		c = (char)tolower(c);
	if (order_method.size() == 0)
		order_method = "amd";
	if ((order_method != "amd") && (order_method != "natural"))
		throw runtime_error("EnifGraph::from_file(): unknown ies_enif_order '" +
			order_method + "', should be 'amd' or 'natural'");

	//a positive definite matrix carrying the graph pattern.  gershgorin: with
	//the diagonal at max degree + 1 the matrix is strictly diagonally dominant,
	//so the factorisation cannot fail for numerical reasons and the pattern it
	//returns is the symbolic one
	int maxdeg = 0;
	for (int i = 0; i < p; i++)
		maxdeg = max(maxdeg, (int)adj.col(i).nonZeros() - 1);
	Eigen::SparseMatrix<double> pattern = adj;
	for (int k = 0; k < pattern.outerSize(); k++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, k); it; ++it)
			it.valueRef() = (it.row() == it.col()) ? (double)(maxdeg + 1) : -1.0;

	Eigen::SparseMatrix<double> Lp;
	vector<int> orig_of_pos(p);
	if (order_method == "natural")
	{
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
			Eigen::NaturalOrdering<int>> sym;
		sym.compute(pattern);
		if (sym.info() != Eigen::Success)
			throw runtime_error("EnifGraph::from_file(): symbolic factorisation failed");
		Lp = sym.matrixL();
		for (int i = 0; i < p; i++)
			orig_of_pos[i] = i;
	}
	else
	{
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
			Eigen::AMDOrdering<int>> sym;
		sym.compute(pattern);
		if (sym.info() != Eigen::Success)
			throw runtime_error("EnifGraph::from_file(): symbolic factorisation failed");
		Lp = sym.matrixL();
		//P maps original -> elimination position; invert it so the solve can walk
		//positions and recover which parameter each one is
		const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>& P = sym.permutationP();
		for (int i = 0; i < p; i++)
			orig_of_pos[P.indices()(i)] = i;
	}

	//L is lower triangular in elimination space: a non-zero at (row, col) with
	//row > col says the node at position 'row' is regressed on the node at
	//position 'col'.  carry that back to parameter indices.
	solve_order = orig_of_pos;
	pred_sets.assign(p, vector<int>());
	int lnnz = 0;
	for (int c = 0; c < Lp.outerSize(); c++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Lp, c); it; ++it)
		{
			int r = (int)it.row(), cc = (int)it.col();
			lnnz++;
			if (r > cc)
				pred_sets[orig_of_pos[r]].push_back(orig_of_pos[cc]);
		}
	fill_edges = (lnnz - p) - (nnz_offdiag / 2);

	//the same sets without the fill: direct graph neighbours that come earlier in the
	//ordering.  pos_of_orig inverts solve_order
	vector<int> pos_of_orig(p);
	for (int pos = 0; pos < p; pos++)
		pos_of_orig[orig_of_pos[pos]] = pos;
	direct_pred_sets.assign(p, vector<int>());
	int max_direct = 0;
	for (int i = 0; i < p; i++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(adj, i); it; ++it)
		{
			int j = (int)it.row();
			if ((j != i) && (pos_of_orig[j] < pos_of_orig[i]))
				direct_pred_sets[i].push_back(j);
		}
		sort(direct_pred_sets[i].begin(), direct_pred_sets[i].end());
		max_direct = max(max_direct, (int)direct_pred_sets[i].size());
	}

	initialized = true;
	report(frec);
	frec << "...solve ordering: " << order_method << "; cholesky factor carries "
		<< (lnnz - p) << " off-diagonal entries, " << fill_edges << " of them fill "
		<< "beyond the " << (nnz_offdiag / 2) << " graph edges" << endl;
	frec << "...largest predecessor set with fill: " << [&]() { size_t m = 0; for (auto& v : pred_sets) m = max(m, v.size()); return m; }()
		<< ", direct neighbours only: " << max_direct << endl;
}


void EnifGraph::report(ofstream& frec) const
{
	if (!initialized)
		return;
	int p = (int)names.size();
	//degree statistics tell the user whether the graph they supplied is sane
	int mind = p, maxd = 0;
	double sumd = 0.0;
	for (int i = 0; i < p; i++)
	{
		int d = (int)(adj.col(i).nonZeros()) - 1;
		mind = min(mind, d);
		maxd = max(maxd, d);
		sumd += d;
	}
	frec << endl << "  ---  enif conditional-independence graph  ---  " << endl;
	frec << "...nodes (adjustable parameters): " << p << endl;
	frec << "...edges: " << num_edges() << endl;
	frec << "...neighbours per node: min " << mind << ", mean "
		<< setprecision(3) << (sumd / (double)p) << ", max " << maxd << endl;
	frec << "...density: " << setprecision(4)
		<< (100.0 * (double)adj.nonZeros() / ((double)p * (double)p)) << " percent" << endl;
	if (maxd == 0)
		frec << "...WARNING: the graph has no edges - the prior precision will be "
		<< "diagonal, which discards all parameter correlation" << endl;
}


void EnifGraph::estimate_precision(const Eigen::MatrixXd& anomalies, double shrink,
	ofstream& frec, bool direct_only, bool into_current)
{
	stringstream ss;
	if (!initialized)
		throw runtime_error("EnifGraph::estimate_precision(): no graph");
	int p = (int)names.size();
	if (anomalies.rows() != p)
		throw runtime_error("EnifGraph::estimate_precision(): anomaly rows != graph nodes");
	int nreal = (int)anomalies.cols();

	//Estimate the precision through its cholesky-like factor, which guarantees a
	//symmetric positive definite result.  For node i, regress its anomalies on the
	//anomalies of the graph neighbours that PRECEDE it in the solve ordering:
	//
	//    a_i = sum_{j in pred(i)} b_j a_j + e_i,   var(e_i) = d_i
	//
	//then the row of the upper-triangular factor is L(i,i) = 1/sqrt(d_i) and
	//L(i,j) = -b_j/sqrt(d_i), and Lam = L^T L.  Restricting to predecessors is
	//what makes the factor triangular; the resulting precision carries the full
	//(symmetric) graph sparsity.
	vector<Eigen::Triplet<double>> ltrips;
	ltrips.reserve(adj.nonZeros());
	int n_shrunk = 0, max_pred = 0, n_floored = 0;
	//the ridge each node's regression actually got, for the rec
	vector<double> ridges;
	//shrink < 0 is the default: pick each node's ridge by k-fold cross-validation over
	//the realizations.  the ridge exists to stop the regression over-fitting the
	//ensemble, and held-out prediction error measures exactly that, so it is chosen by
	//the thing it is for rather than by a formula.  the candidates are relative to the
	//mean diagonal of the local gram; folds are contiguous blocks of realizations.  too
	//few realizations to hold any out falls back to k/(N-1) per node
	const bool cv_ridge = (shrink < 0.0);
	const int cv_min_reals = 20;
	const int cv_folds = 5;
	const vector<double> cv_cands = {1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1, 3.0e-1, 1.0};
	int n_cv_fallback = 0;

	for (int pos = 0; pos < p; pos++)
	{
		//walk the elimination ordering, not the parameter ordering
		int i = solve_order.empty() ? pos : solve_order[pos];
		vector<int> pred = direct_only ? (direct_pred_sets.empty() ? vector<int>() : direct_pred_sets[i])
			: (pred_sets.empty() ? vector<int>() : pred_sets[i]);

		//A neighbourhood anywhere near the ensemble size overfits: the regression
		//interpolates in-sample, the residual variance collapses, and since
		//Lam_ii = 1/d_i the prior precision explodes.  An exploded prior precision
		//is not a loud failure - it is infinite prior confidence, so the update
		//silently goes to zero and the ensemble never moves.  Cap the
		//neighbourhood at half the ensemble to keep the regression honest.
		int max_allowed = max(1, (nreal - 1) / 2);
		if ((int)pred.size() > max_allowed)
		{
			//keep the strongest predecessors, not the lowest-numbered ones.
			//truncating by index is a geometric bias on a grid - it keeps
			//whichever neighbours happen to come first in the file
			vector<pair<double, int>> strength;
			strength.reserve(pred.size());
			double ni = anomalies.row(i).norm();
			for (int j : pred)
			{
				double nj = anomalies.row(j).norm();
				double den = ni * nj;
				double c = (den > 0.0) ? abs(anomalies.row(i).dot(anomalies.row(j))) / den : 0.0;
				strength.push_back(make_pair(c, j));
			}
			sort(strength.begin(), strength.end(),
				[](const pair<double, int>& a, const pair<double, int>& b)
				{ return a.first > b.first; });
			pred.clear();
			for (int k = 0; k < max_allowed; k++)
				pred.push_back(strength[k].second);
			sort(pred.begin(), pred.end());
			n_shrunk++;
		}
		max_pred = max(max_pred, (int)pred.size());

		double var_i = anomalies.row(i).squaredNorm();
		double d_i = var_i;
		Eigen::VectorXd b;

		if (!pred.empty())
		{
			Eigen::MatrixXd Ap(pred.size(), nreal);
			for (size_t k = 0; k < pred.size(); k++)
				Ap.row(k) = anomalies.row(pred[k]);
			Eigen::MatrixXd G = Ap * Ap.transpose();
			//stein-type pull toward the diagonal, strengthened as the
			//neighbourhood grows relative to the ensemble.  a fixed shrinkage is
			//not enough: what matters is k/N, and the penalty has to bite hardest
			//exactly where the local gram is worst conditioned.
			//
			//shrink < 0 (the default) means automatic: the ridge is k/(N-1) itself,
			//the number of predictors over the samples, per node.  on the hkpp
			//problem with the precision from the current ensemble a fixed 1e-3 let
			//the regressions over-fit a contracted ensemble, the hessian claimed
			//more certainty than the ensemble had, and only heavily damped steps
			//still reduced phi; a sweep put the useful value near k/(N-1) (0.1 for
			//the rook graph at N=50) with 0.5 already too much.  a positive shrink
			//is the old behaviour: a user floor, with kfrac^2 underneath it.
			double kfrac = (double)pred.size() / (double)max(1, nreal - 1);
			double tr = G.trace() / (double)pred.size();
			double eff_shrink;
			if (!cv_ridge)
				eff_shrink = max(shrink, kfrac * kfrac);
			else if (nreal < cv_min_reals)
			{
				eff_shrink = kfrac;
				n_cv_fallback++;
			}
			else
			{
				//k-fold: for each candidate, fit on the other folds, score the held-out
				//columns, keep the candidate with the smallest total held-out error
				Eigen::VectorXd ai = anomalies.row(i).transpose();
				double best_err = -1.0;
				eff_shrink = kfrac;
				for (double cand : cv_cands)
				{
					double err = 0.0;
					for (int f = 0; f < cv_folds; f++)
					{
						int lo = (f * nreal) / cv_folds, hi = ((f + 1) * nreal) / cv_folds;
						int ntr = nreal - (hi - lo);
						if ((ntr < 2) || (hi <= lo))
							continue;
						Eigen::MatrixXd Atr(pred.size(), ntr), Ate(pred.size(), hi - lo);
						Eigen::VectorXd btr(ntr), bte(hi - lo);
						int c = 0;
						for (int j = 0; j < nreal; j++)
						{
							if ((j >= lo) && (j < hi))
							{
								Ate.col(j - lo) = Ap.col(j);
								bte[j - lo] = ai[j];
							}
							else
							{
								Atr.col(c) = Ap.col(j);
								btr[c] = ai[j];
								c++;
							}
						}
						Eigen::MatrixXd Gtr = Atr * Atr.transpose();
						double trtr = Gtr.trace() / (double)pred.size();
						Gtr.diagonal().array() += max(cand * trtr, 1.0e-12 * max(trtr, 1.0));
						Eigen::VectorXd btrial = Gtr.ldlt().solve(Atr * btr);
						err += (bte - Ate.transpose() * btrial).squaredNorm();
					}
					if ((best_err < 0.0) || (err < best_err))
					{
						best_err = err;
						eff_shrink = cand;
					}
				}
			}
			ridges.push_back(eff_shrink);
			G.diagonal().array() += max(eff_shrink * tr, 1.0e-12 * max(tr, 1.0));
			Eigen::VectorXd rhs = Ap * anomalies.row(i).transpose();
			b = G.ldlt().solve(rhs);
			d_i = (anomalies.row(i).transpose() - Ap.transpose() * b).squaredNorm();
		}

		//floor the residual variance against the node's OWN variance, not an
		//absolute constant: this bounds Lam_ii at RESID_FLOOR^-1 times the
		//diagonal precision, so no node can claim near-infinite certainty
		const double RESID_FLOOR = 1.0e-3;
		double floor_d = RESID_FLOOR * max(var_i, 1.0e-30);
		if (!(d_i > floor_d))
		{
			d_i = floor_d;
			n_floored++;
		}

		double s = 1.0 / sqrt(d_i);
		ltrips.push_back(Eigen::Triplet<double>(i, i, s));
		for (size_t k = 0; k < pred.size(); k++)
			if (b[k] != 0.0)
				ltrips.push_back(Eigen::Triplet<double>(i, pred[k], -b[k] * s));
	}

	Eigen::SparseMatrix<double> L(p, p);
	L.setFromTriplets(ltrips.begin(), ltrips.end());
	if (into_current)
	{
		//the damping metric for the hybrid step, re-estimated each iteration.  it is
		//never factored on its own, only as part of the posterior precision
		prec_cur = (Eigen::SparseMatrix<double>)(L.transpose() * L);
		prec_cur.makeCompressed();
		cur_ready = true;
		frec << "...estimated current-ensemble precision on the graph for the hybrid damping: "
			<< prec_cur.nonZeros() << " non-zeros" << endl;
		return;
	}
	prec = (Eigen::SparseMatrix<double>)(L.transpose() * L);
	prec.makeCompressed();

	solver.compute(prec);
	if (solver.info() != Eigen::Success)
		throw runtime_error("EnifGraph::estimate_precision(): failed to factor the "
			"estimated prior precision - try a sparser graph or more shrinkage");
	prec_ready = true;

	frec << "...estimated prior precision on the graph: " << prec.nonZeros()
		<< " non-zeros, density " << setprecision(4)
		<< (100.0 * (double)prec.nonZeros() / ((double)p * (double)p)) << " percent" << endl;
	frec << "...largest neighbourhood used: " << max_pred << " of " << nreal
		<< " realizations" << (direct_only ? "  (direct graph neighbours only, fill ignored: the factor is "
			"not exact, the regressions are sized to the ensemble)" : "") << endl;
	if (!ridges.empty())
	{
		sort(ridges.begin(), ridges.end());
		string how = cv_ridge ? (n_cv_fallback > 0 ? "k/(N-1), too few realizations to cross-validate"
			: "cross-validated, " + to_string(cv_folds) + " folds") : "user floor " + to_string(shrink);
		frec << "...regression ridge (" << how << ") min / median / max over nodes: " << setprecision(3)
			<< ridges.front() << " / " << ridges[ridges.size() / 2] << " / " << ridges.back() << endl;
	}
	if (n_shrunk > 0)
		frec << "...WARNING: " << n_shrunk << " nodes had more neighbours than the "
		<< "ensemble can support and were truncated to " << ((nreal - 1) / 2)
		<< " - consider a sparser graph or more realizations" << endl;
	if (n_floored > 0)
		frec << "...WARNING: " << n_floored << " nodes hit the residual-variance "
		<< "floor, meaning their neighbours explain them almost perfectly in "
		<< "sample.  the prior precision there is capped; a sparser graph would "
		<< "be a better description of this ensemble" << endl;
}


Eigen::MatrixXd EnifGraph::apply_cov(const Eigen::MatrixXd& M) const
{
	if (!prec_ready)
		throw runtime_error("EnifGraph::apply_cov(): precision not estimated");
	//C * M == solve(Lam, M): the covariance implied by the graph is never formed
	Eigen::MatrixXd out = solver.solve(M);
	if (solver.info() != Eigen::Success)
		throw runtime_error("EnifGraph::apply_cov(): sparse solve failed");
	return out;
}


void EnifGraph::load_precision(const string& filename, ofstream& frec)
{
	if (!initialized)
		throw runtime_error("EnifGraph::load_precision(): call from_file() with the same file first");
	int p = (int)names.size();
	Mat m;
	m.from_file(filename);
	Mat sub = m.get(names, names);
	Eigen::SparseMatrix<double> a = sub.get_matrix();
	if ((a.rows() != p) || (a.cols() != p))
		throw runtime_error("EnifGraph::load_precision(): precision file '" + filename + "' is not square over the adjustable parameters");
	Eigen::SparseMatrix<double> at = a.transpose();
	Eigen::SparseMatrix<double> diff = (a - at).eval();
	double asym = 0.0, amax = 0.0;
	for (int k = 0; k < diff.outerSize(); k++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(diff, k); it; ++it)
			asym = max(asym, abs(it.value()));
	for (int k = 0; k < a.outerSize(); k++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(a, k); it; ++it)
			amax = max(amax, abs(it.value()));
	if (asym > 1.0e-8 * max(amax, 1.0e-300))
		throw runtime_error("EnifGraph::load_precision(): precision file '" + filename + "' is not symmetric");
	prec = (0.5 * (a + at)).eval();
	prec.makeCompressed();
	for (int i = 0; i < p; i++)
		if (!(prec.coeff(i, i) > 0.0))
			throw runtime_error("EnifGraph::load_precision(): precision file '" + filename + "' has a non-positive diagonal at '" + names[i] + "'");
	solver.compute(prec);
	if (solver.info() != Eigen::Success)
		throw runtime_error("EnifGraph::load_precision(): precision file '" + filename + "' is not positive definite");
	prec_ready = true;
	frec << "...prior precision read from '" << filename << "': " << prec.nonZeros()
		<< " non-zeros, density " << setprecision(4)
		<< (100.0 * (double)prec.nonZeros() / ((double)p * (double)p)) << " percent" << endl;
}


Eigen::MatrixXd EnifGraph::information_step(const Eigen::SparseMatrix<double>& H,
	const Eigen::VectorXd& rinv, const Eigen::MatrixXd& e,
	const Eigen::MatrixXd& resid, double lam, ostream& frec) const
{
	if (!prec_ready)
		throw runtime_error("EnifGraph::information_step(): precision not estimated");
	int p = (int)names.size();

	//H^T Rinv H stays sparse because H is sparse - that is the whole reason the
	//graph is worth having.  with a dense H this product would be a full p x p
	//matrix and the sparse cholesky below would be pointless.
	Eigen::SparseMatrix<double> Rinv((int)H.rows(), (int)H.rows());
	{
		vector<Eigen::Triplet<double>> t;
		t.reserve(H.rows());
		for (int i = 0; i < (int)H.rows(); i++)
			t.push_back(Eigen::Triplet<double>(i, i, rinv[i]));
		Rinv.setFromTriplets(t.begin(), t.end());
	}
	Eigen::SparseMatrix<double> HtRH = (H.transpose() * Rinv * H).eval();

	Eigen::SparseMatrix<double> lam_post;
	if (hybrid_damp)
	{
		//lambda on the current-ensemble precision: the curvature keeps the prior, the
		//trust region follows the ensemble
		if (!cur_ready)
			throw runtime_error("EnifGraph::information_step(): hybrid damping asked for but no current-ensemble precision estimated");
		lam_post = (prec + (lam * prec_cur) + HtRH).eval();
	}
	else
		lam_post = ((1.0 + lam) * prec + HtRH).eval();
	lam_post.makeCompressed();

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> post_solver;
	post_solver.compute(lam_post);
	if (post_solver.info() != Eigen::Success)
		throw runtime_error("EnifGraph::information_step(): failed to factor the "
			"posterior precision");

	//gradient of the rml objective at the current iterate
	Eigen::MatrixXd wresid = rinv.asDiagonal() * resid;
	Eigen::MatrixXd g = prec * e + H.transpose() * wresid;
	Eigen::MatrixXd delta = -post_solver.solve(g);
	if (post_solver.info() != Eigen::Success)
		throw runtime_error("EnifGraph::information_step(): posterior solve failed");

	frec << "...posterior precision: " << lam_post.nonZeros() << " non-zeros, density "
		<< setprecision(4)
		<< (100.0 * (double)lam_post.nonZeros() / ((double)p * (double)p))
		<< " percent (H contributed " << HtRH.nonZeros() << ")" << endl;
	return delta;
}


//one coordinate descent sweep over the columns listed in cols for
//    min 0.5||b - X x||^2 + alpha ||x||_1
//X is realizations x parameters, so each parameter is a contiguous column; cnorm holds
//the squared column norms.  x and the residual r = b - X x are updated in place.
//returns the largest coefficient change.
static double lasso_sweep(const Eigen::MatrixXd& X, const Eigen::VectorXd& cnorm,
	double alpha, Eigen::VectorXd& x, Eigen::VectorXd& r, const vector<int>& cols)
{
	double max_chg = 0.0;
	for (int j : cols)
	{
		if (cnorm[j] <= 0.0)
			continue;
		double xj = x[j];
		double rho = X.col(j).dot(r) + cnorm[j] * xj;
		double nx = 0.0;
		if (rho > alpha)
			nx = (rho - alpha) / cnorm[j];
		else if (rho < -alpha)
			nx = (rho + alpha) / cnorm[j];
		if (nx != xj)
		{
			r.noalias() -= X.col(j) * (nx - xj);
			x[j] = nx;
			max_chg = max(max_chg, abs(nx - xj));
		}
	}
	return max_chg;
}


//coordinate descent to convergence with an active set: a full sweep finds which
//coefficients are nonzero, sweeps over just those run until they settle, and the next
//full sweep checks nothing outside them wants in.  warm starts from x and r.
static void lasso_cd(const Eigen::MatrixXd& X, const Eigen::VectorXd& cnorm,
	double alpha, Eigen::VectorXd& x, Eigen::VectorXd& r, const vector<int>& all_cols,
	int max_sweeps, double tol)
{
	vector<int> active;
	for (int outer = 0; outer < max_sweeps; outer++)
	{
		double chg = lasso_sweep(X, cnorm, alpha, x, r, all_cols);
		if (chg < tol)
			break;
		active.clear();
		for (int j : all_cols)
			if (x[j] != 0.0)
				active.push_back(j);
		for (int inner = 0; inner < max_sweeps; inner++)
			if (lasso_sweep(X, cnorm, alpha, x, r, active) < tol)
				break;
	}
}


//the reference implementation's H: standardized rows, penalty per observation by
//k-fold cross-validation, then a refit on every realization
static Eigen::SparseMatrix<double> estimate_sparse_H_cv(const Eigen::MatrixXd& A,
	const Eigen::MatrixXd& B, int cv_folds, int num_threads,
	Eigen::VectorXd& unexplained, ostream& frec, double unexp_divisor)
{
	int p = (int)A.rows();
	int nreal = (int)A.cols();
	int nobs = (int)B.rows();
	unexplained.setZero(nobs);
	if (unexp_divisor <= 0.0)
		unexp_divisor = (double)nreal;
	int nfold = max(2, min(cv_folds, nreal));
	const int n_alphas = 50;
	const double alpha_eps = 1.0e-3;
	const int max_sweeps = 1000;
	//coefficients are in standardized units, so this is a relative tolerance
	const double tol = 1.0e-6;

	//unit-length parameter anomalies, stored realizations x parameters so each
	//parameter is a contiguous column.  a parameter with no spread stays zero
	Eigen::VectorXd asd = A.rowwise().norm();
	Eigen::MatrixXd Xs(nreal, p);
	for (int j = 0; j < p; j++)
		Xs.col(j) = (asd[j] > 0.0) ? Eigen::VectorXd(A.row(j).transpose() / asd[j]) : Eigen::VectorXd::Zero(nreal);
	vector<int> all_cols(p);
	for (int j = 0; j < p; j++)
		all_cols[j] = j;

	//contiguous folds over the realizations, as a plain k-fold split does
	vector<vector<int>> test_idx(nfold), train_idx(nfold);
	for (int k = 0; k < nreal; k++)
		test_idx[(int)((long)k * nfold / nreal)].push_back(k);
	vector<Eigen::MatrixXd> Xtr(nfold), Xte(nfold);
	vector<Eigen::VectorXd> xtr_norm(nfold);
	for (int f = 0; f < nfold; f++)
	{
		set<int> te(test_idx[f].begin(), test_idx[f].end());
		for (int k = 0; k < nreal; k++)
			if (te.find(k) == te.end())
				train_idx[f].push_back(k);
		Xtr[f].resize(train_idx[f].size(), p);
		Xte[f].resize(test_idx[f].size(), p);
		for (size_t c = 0; c < train_idx[f].size(); c++)
			Xtr[f].row(c) = Xs.row(train_idx[f][c]);
		for (size_t c = 0; c < test_idx[f].size(); c++)
			Xte[f].row(c) = Xs.row(test_idx[f][c]);
		xtr_norm[f] = Xtr[f].colwise().squaredNorm().transpose();
	}
	Eigen::VectorXd xs_norm = Xs.colwise().squaredNorm().transpose();

	vector<vector<Eigen::Triplet<double>>> row_trips(nobs);
	vector<double> chosen_frac(nobs, 0.0);

	auto solve_row = [&](int i)
	{
		double bsd = B.row(i).norm();
		if (bsd <= 0.0)
			return;
		Eigen::VectorXd bs = B.row(i).transpose() / bsd;

		//penalty path from the value that zeros the whole row down to alpha_eps of it
		double amax = (Xs.transpose() * bs).cwiseAbs().maxCoeff();
		if (amax <= 0.0)
			return;
		vector<double> alphas(n_alphas);
		for (int a = 0; a < n_alphas; a++)
			alphas[a] = amax * pow(alpha_eps, (double)a / (double)(n_alphas - 1));

		//held-out squared error summed over folds, per penalty
		vector<double> cv_err(n_alphas, 0.0);
		for (int f = 0; f < nfold; f++)
		{
			Eigen::VectorXd btr(train_idx[f].size()), bte(test_idx[f].size());
			for (size_t c = 0; c < train_idx[f].size(); c++)
				btr[c] = bs[train_idx[f][c]];
			for (size_t c = 0; c < test_idx[f].size(); c++)
				bte[c] = bs[test_idx[f][c]];
			//the penalty is on the sum of squares, so scale it to the training size
			//to keep the path comparable to the full fit
			double nscale = (double)train_idx[f].size() / (double)nreal;
			Eigen::VectorXd x = Eigen::VectorXd::Zero(p);
			Eigen::VectorXd r = btr;
			for (int a = 0; a < n_alphas; a++)
			{
				lasso_cd(Xtr[f], xtr_norm[f], alphas[a] * nscale, x, r, all_cols, max_sweeps, tol);
				cv_err[a] += (bte - Xte[f] * x).squaredNorm();
			}
		}
		int best = (int)(min_element(cv_err.begin(), cv_err.end()) - cv_err.begin());
		chosen_frac[i] = alphas[best] / amax;

		//refit on every realization, warm starting down the path to the chosen penalty
		Eigen::VectorXd x = Eigen::VectorXd::Zero(p);
		Eigen::VectorXd r = bs;
		for (int a = 0; a <= best; a++)
			lasso_cd(Xs, xs_norm, alphas[a], x, r, all_cols, max_sweeps, tol);

		//back to original units: h_ij = x_j * ||b_i|| / ||a_j||, and the residual
		//variance in original units is what the observation error is inflated by
		Eigen::VectorXd h = Eigen::VectorXd::Zero(p);
		for (int j = 0; j < p; j++)
			if ((x[j] != 0.0) && (asd[j] > 0.0))
				h[j] = x[j] * bsd / asd[j];
		Eigen::VectorXd res = B.row(i).transpose() - A.transpose() * h;
		unexplained[i] = res.squaredNorm() * (double)(nreal - 1) / unexp_divisor;
		for (int j = 0; j < p; j++)
			if (h[j] != 0.0)
				row_trips[i].push_back(Eigen::Triplet<double>(i, j, h[j]));
	};

	//ies_num_threads defaults to -1; take every core then, the rows are independent
	if (num_threads < 1)
		num_threads = max(1, (int)std::thread::hardware_concurrency());
	if (num_threads < 2)
	{
		for (int i = 0; i < nobs; i++)
			solve_row(i);
	}
	else
	{
		vector<thread> threads;
		std::atomic<int> next(0);
		int nt = min(num_threads, nobs);
		for (int t = 0; t < nt; t++)
			threads.push_back(thread([&]() {
				int i;
				while ((i = next.fetch_add(1)) < nobs)
					solve_row(i);
				}));
		for (auto& th : threads)
			th.join();
	}

	vector<Eigen::Triplet<double>> all;
	for (auto& rt : row_trips)
		all.insert(all.end(), rt.begin(), rt.end());
	Eigen::SparseMatrix<double> H(nobs, p);
	H.setFromTriplets(all.begin(), all.end());
	H.makeCompressed();

	vector<double> fr = chosen_frac;
	sort(fr.begin(), fr.end());
	frec << "...sparse H by " << nfold << "-fold cross-validated lasso on standardized anomalies: "
		<< H.nonZeros() << " of " << ((long)nobs * (long)p) << " entries (" << setprecision(3)
		<< (100.0 * (double)H.nonZeros() / ((double)nobs * (double)p)) << " percent)" << endl;
	frec << "...chosen penalty as a fraction of the row-zeroing value, min / median / max: "
		<< fr.front() << " / " << fr[fr.size() / 2] << " / " << fr.back() << endl;
	frec << "...mean unexplained variance per observation: " << unexplained.mean() << endl;
	return H;
}


Eigen::SparseMatrix<double> estimate_sparse_H(const Eigen::MatrixXd& A,
	const Eigen::MatrixXd& B, double lasso_frac, int num_threads,
	Eigen::VectorXd& unexplained, ostream& frec, int cv_folds, double unexp_divisor)
{
	if (cv_folds > 0)
		return estimate_sparse_H_cv(A, B, cv_folds, num_threads, unexplained, frec, unexp_divisor);

	int p = (int)A.rows();
	int nreal = (int)A.cols();
	int nobs = (int)B.rows();
	unexplained.setZero(nobs);
	if (unexp_divisor <= 0.0)
		unexp_divisor = (double)nreal;

	//squared norms of each parameter's anomaly vector - the coordinate descent
	//denominators, computed once
	Eigen::VectorXd anorm(p);
	for (int j = 0; j < p; j++)
		anorm[j] = A.row(j).squaredNorm();

	vector<vector<Eigen::Triplet<double>>> row_trips(nobs);

	//each observation row is an independent lasso problem
	auto solve_row = [&](int i)
	{
		Eigen::VectorXd b = B.row(i).transpose();
		//penalty scaled by the smallest value that would zero the whole row, so
		//lasso_frac is dimensionless
		double amax = 0.0;
		for (int j = 0; j < p; j++)
			amax = max(amax, abs(A.row(j).dot(b)));
		double alpha = lasso_frac * amax;

		Eigen::VectorXd x = Eigen::VectorXd::Zero(p);
		Eigen::VectorXd r = b;
		for (int sweep = 0; sweep < 100; sweep++)
		{
			double max_chg = 0.0;
			for (int j = 0; j < p; j++)
			{
				if (anorm[j] <= 0.0)
					continue;
				double xj = x[j];
				double rho = A.row(j).dot(r) + anorm[j] * xj;
				double nx = 0.0;
				if (rho > alpha)
					nx = (rho - alpha) / anorm[j];
				else if (rho < -alpha)
					nx = (rho + alpha) / anorm[j];
				if (nx != xj)
				{
					r -= A.row(j).transpose() * (nx - xj);
					x[j] = nx;
					max_chg = max(max_chg, abs(nx - xj));
				}
			}
			if (max_chg < 1.0e-10)
				break;
		}
		//residual variance is what the observation error gets inflated by.  the
		//anomalies carry a 1/sqrt(nreal-1), so r^2 (nreal-1) is the raw sum of squares
		unexplained[i] = r.squaredNorm() * (double)(nreal - 1) / unexp_divisor;
		for (int j = 0; j < p; j++)
			if (x[j] != 0.0)
				row_trips[i].push_back(Eigen::Triplet<double>(i, j, x[j]));
	};

	if (num_threads < 2)
	{
		for (int i = 0; i < nobs; i++)
			solve_row(i);
	}
	else
	{
		vector<thread> threads;
		std::atomic<int> next(0);
		int nt = min(num_threads, nobs);
		for (int t = 0; t < nt; t++)
			threads.push_back(thread([&]() {
				int i;
				while ((i = next.fetch_add(1)) < nobs)
					solve_row(i);
				}));
		for (auto& th : threads)
			th.join();
	}

	vector<Eigen::Triplet<double>> all;
	for (auto& rt : row_trips)
		all.insert(all.end(), rt.begin(), rt.end());
	Eigen::SparseMatrix<double> H(nobs, p);
	H.setFromTriplets(all.begin(), all.end());
	H.makeCompressed();

	frec << "...sparse H: " << H.nonZeros() << " of " << ((long)nobs * (long)p)
		<< " entries (" << setprecision(3)
		<< (100.0 * (double)H.nonZeros() / ((double)nobs * (double)p))
		<< " percent), lasso fraction " << lasso_frac << endl;
	frec << "...mean unexplained variance per observation: "
		<< unexplained.mean() << endl;
	return H;
}


Eigen::MatrixXd enif_woodbury_step(const Eigen::SparseMatrix<double>& H,
	const Eigen::VectorXd& noise_var, const Eigen::SparseMatrix<double>& C,
	const Eigen::MatrixXd& e, const Eigen::MatrixXd& resid, double lam)
{
	int n = (int)H.rows();
	int p = (int)H.cols();
	if ((C.rows() != p) || (C.cols() != p))
		throw runtime_error("enif_woodbury_step(): prior covariance is not p x p");
	if ((noise_var.size() != n) || (resid.rows() != n) || (e.rows() != p))
		throw runtime_error("enif_woodbury_step(): H, noise, e and resid sizes disagree");
	double s = 1.0 / (1.0 + lam);
	//C_lam H^T: only the non-zero rows of H touch C
	Eigen::SparseMatrix<double> CHt_s = (C * H.transpose()).eval();
	Eigen::MatrixXd CHt = s * Eigen::MatrixXd(CHt_s);
	Eigen::MatrixXd G = Eigen::MatrixXd(H * CHt);
	G.diagonal() += noise_var;
	Eigen::LDLT<Eigen::MatrixXd> g_fact(G);
	if (g_fact.info() != Eigen::Success)
		throw runtime_error("enif_woodbury_step(): failed to factor the innovation covariance");
	Eigen::MatrixXd He = H * e;
	return -((s * (e - (CHt * g_fact.solve(He)))) + (CHt * g_fact.solve(resid)));
}


Eigen::MatrixXd enif_lowrank_step(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Q,
	const Eigen::MatrixXd& L, const Eigen::VectorXd& noise_var,
	const Eigen::MatrixXd& e, const Eigen::MatrixXd& He, const Eigen::MatrixXd& resid, double lam)
{
	int n = (int)U.rows();
	int N = (int)U.cols();
	int p = (int)L.rows();
	if ((Q.rows() != N) || (Q.cols() != N) || (L.cols() != N))
		throw runtime_error("enif_lowrank_step(): U, Q and L do not share the ensemble dimension");
	if ((noise_var.size() != n) || (resid.rows() != n) || (He.rows() != n) || (e.rows() != p))
		throw runtime_error("enif_lowrank_step(): U, L, noise, e, He and resid sizes disagree");
	if ((e.cols() != resid.cols()) || (He.cols() != resid.cols()))
		throw runtime_error("enif_lowrank_step(): e, He and resid must have the same number of columns");
	for (int i = 0; i < n; i++)
		if (!(noise_var[i] > 0.0))
			throw runtime_error("enif_lowrank_step(): noise variance must be positive");
	double s = 1.0 / (1.0 + lam);
	Eigen::VectorXd rinv = noise_var.cwiseInverse();
	//Rinv U once; then S = U^T Rinv U and the N x N system (I + S sQ)
	Eigen::MatrixXd RU = rinv.asDiagonal() * U;
	Eigen::MatrixXd S = U.transpose() * RU;
	Eigen::MatrixXd Qs = s * Q;
	Eigen::MatrixXd Kmat = Eigen::MatrixXd::Identity(N, N) + S * Qs;
	Eigen::PartialPivLU<Eigen::MatrixXd> k_fact(Kmat);
	//W = resid - s He: the prior term carries s once for the (1+lam) on it and once
	//more inside C_lam, so the He part of the gain is s^2 against s for the residual
	Eigen::MatrixXd W = resid - (s * He);
	Eigen::MatrixXd UtRW = RU.transpose() * W;
	Eigen::MatrixXd GinvW = (rinv.asDiagonal() * W) - (RU * (Qs * k_fact.solve(UtRW)));
	if (!GinvW.allFinite())
		throw runtime_error("enif_lowrank_step(): the innovation solve is not finite");
	return -s * (e + (L * (U.transpose() * GinvW)));
}


map<string, double> enif_mm_inflation_csv(const vector<string>& real_names,
	const vector<string>& obs_names, const vector<string>& groups,
	const vector<Eigen::VectorXd>& weights, const vector<Eigen::VectorXd>& unexplained,
	const vector<const Eigen::SparseMatrix<double>*>& H, bool applied,
	const string& csv_filename)
{
	int nreal = (int)real_names.size();
	int nobs = (int)obs_names.size();
	if (((int)weights.size() != nreal) || ((int)unexplained.size() != nreal) || ((int)H.size() != nreal))
		throw runtime_error("enif_mm_inflation_csv(): per-realization vectors are not aligned with real_names");
	if ((int)groups.size() != nobs)
		throw runtime_error("enif_mm_inflation_csv(): groups not aligned with obs_names");
	ofstream csv(csv_filename);
	if (!csv.good())
		throw runtime_error("enif_mm_inflation_csv(): error opening " + csv_filename);
	csv << "real_name,obs_name,group,weight,noise_var,unexplained_var,inflated_var,inflate_ratio,effective_weight,h_row_nnz" << endl;
	csv << setprecision(10);
	map<string, double> mean_ratio;
	for (int k = 0; k < nreal; k++)
	{
		const Eigen::VectorXd& w = weights[k];
		const Eigen::VectorXd& u = unexplained[k];
		if ((w.size() != nobs) || (u.size() != nobs) || (H[k] == nullptr) || (H[k]->rows() != nobs))
			throw runtime_error("enif_mm_inflation_csv(): realization '" + real_names[k] + "' vectors not aligned with obs_names");
		//row counts of a compressed row-major or column-major sparse matrix: walk the non-zeros once
		vector<int> row_nnz(nobs, 0);
		for (int oc = 0; oc < H[k]->outerSize(); oc++)
			for (Eigen::SparseMatrix<double>::InnerIterator it(*H[k], oc); it; ++it)
				if (it.value() != 0.0)
					row_nnz[it.row()]++;
		double rsum = 0.0;
		for (int i = 0; i < nobs; i++)
		{
			double noise = (w[i] > 0.0) ? 1.0 / (w[i] * w[i]) : 0.0;
			double uv = max(0.0, u[i]);
			double inflated = noise + uv;
			double ratio = (noise > 0.0) ? inflated / noise : 0.0;
			double eff = applied ? ((inflated > 0.0) ? 1.0 / sqrt(inflated) : 0.0) : w[i];
			csv << real_names[k] << "," << obs_names[i] << "," << groups[i] << "," << w[i] << "," << noise << ","
				<< uv << "," << inflated << "," << ratio << "," << eff << "," << row_nnz[i] << endl;
			rsum += ratio;
		}
		mean_ratio[real_names[k]] = (nobs > 0) ? rsum / (double)nobs : 0.0;
	}
	csv.close();
	return mean_ratio;
}


map<string, EnifInflateGroupStats> enif_inflation_report(
	const vector<string>& obs_names, const vector<string>& groups,
	const Eigen::VectorXd& weights, const Eigen::VectorXd& unexplained,
	bool applied, int iter, const string& csv_filename, ofstream& frec)
{
	int nobs = (int)obs_names.size();
	if (((int)groups.size() != nobs) || (weights.size() != nobs) || (unexplained.size() != nobs))
		throw runtime_error("enif_inflation_report(): obs names, groups, weights and unexplained variance are different lengths");

	//one row per observation.  inflate_ratio is always the ratio the inflation
	//WOULD apply, so a run with ies_enif_resid_inflate false still shows what it
	//is missing.  effective_weight is what the update actually used
	ofstream csv(csv_filename);
	if (!csv.good())
		throw runtime_error("enif_inflation_report(): error opening " + csv_filename);
	csv << "obs_name,group,weight,noise_var,unexplained_var,inflated_var,inflate_ratio,effective_weight" << endl;
	csv << setprecision(10);
	map<string, EnifInflateGroupStats> stats;
	for (int i = 0; i < nobs; i++)
	{
		double w = weights[i];
		double nv = 1.0 / (w * w);
		double uv = max(0.0, unexplained[i]);
		double iv = nv + uv;
		double ratio = iv / nv;
		double used = applied ? iv : nv;
		double ew = 1.0 / sqrt(used);
		csv << obs_names[i] << "," << groups[i] << "," << w << "," << nv << "," << uv << ","
			<< used << "," << ratio << "," << ew << endl;
		EnifInflateGroupStats& s = stats[groups[i]];
		if (s.count == 0)
		{
			s.ratio_min = ratio;
			s.ratio_max = ratio;
		}
		s.count++;
		s.noise_var += nv;
		s.unexp_var += uv;
		s.ratio_mean += ratio;
		s.ratio_min = min(s.ratio_min, ratio);
		s.ratio_max = max(s.ratio_max, ratio);
		s.weight_mean += w;
		s.eff_weight_mean += ew;
	}
	csv.close();
	for (auto& kv : stats)
	{
		EnifInflateGroupStats& s = kv.second;
		s.noise_var /= s.count;
		s.unexp_var /= s.count;
		s.ratio_mean /= s.count;
		s.weight_mean /= s.count;
		s.eff_weight_mean /= s.count;
	}

	//the rec section, in the shape of the group phi summary
	vector<pair<string, double>> order;
	int len = 5;
	for (auto& kv : stats)
	{
		order.push_back(make_pair(kv.first, kv.second.ratio_mean));
		len = max(len, (int)kv.first.size());
	}
	len++;
	sort(order.begin(), order.end(),
		[](const pair<string, double>& a, const pair<string, double>& b) { return a.second > b.second; });
	frec << endl << "  ---  enif observation noise inflation summary, iteration " << iter << "  ---  " << endl;
	frec << "       (noise variance is 1/weight^2, inflated by the variance H fails to explain)" << endl;
	frec << "       (inflation applied to the update: " << (applied ? "yes" : "no, ies_enif_resid_inflate is false") << ")" << endl;
	frec << "           (sorted by mean inflation ratio)" << endl;
	frec << left << setw(len) << "group" << right << setw(7) << "count" << setw(12) << "noise_var"
		<< setw(12) << "unexp_var" << setw(11) << "ratio" << setw(11) << "ratio_min" << setw(11) << "ratio_max"
		<< setw(11) << "weight" << setw(11) << "eff_wght" << endl;
	for (auto& o : order)
	{
		const EnifInflateGroupStats& s = stats[o.first];
		frec << left << setw(len) << pest_utils::lower_cp(o.first) << " ";
		frec << right << setw(6) << s.count << " ";
		frec << right << setw(11) << setprecision(3) << s.noise_var << " ";
		frec << setw(11) << setprecision(3) << s.unexp_var << " ";
		frec << setw(10) << setprecision(3) << s.ratio_mean << " ";
		frec << setw(10) << setprecision(3) << s.ratio_min << " ";
		frec << setw(10) << setprecision(3) << s.ratio_max << " ";
		frec << setw(10) << setprecision(3) << s.weight_mean << " ";
		frec << setw(10) << setprecision(3) << s.eff_weight_mean << endl;
	}
	frec << "    Note: 'ratio' is inflated over original noise variance, 'eff_wght' is the weight the update used." << endl;
	frec << "...saved per-observation inflation to " << csv_filename << endl << endl;
	return stats;
}
