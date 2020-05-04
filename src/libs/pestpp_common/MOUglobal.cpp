#include "MOUglobal.h"
//#include "rand.h"

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace mou;
using namespace std;

extern random_gen rgen; // global common random generator

individual::individual() throw () :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    config(0) {}

individual::individual(const individual_config& c
                       /*const unsigned int nreal,
                       const unsigned int nbin,
                       const unsigned int ncon,
                       const std::vector<int>& nbits,
                       const unsigned int nobj*/) throw (mou::MOEAexception) :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    evaluated(false),
    config(&c) {

    xreal.resize(config->nreal,0);
    xbin.resize(config->nbin,0);
    gene.resize(config->nbin);
    if (config->nbits.size() != config->nbin)
        throw mou::MOEAexception("nbits size != nbin");
    for (int j = 0; j < config->nbin; ++j) {
        gene[j].resize(config->nbits[j],0);
    }
    obj.resize(config->nobj,0);
    constr.resize(config->ncon,0);
}

individual::~individual() {
}

void individual::initialize() throw (mou::MOEAexception) {
    if (!config)
        throw mou::MOEAexception("Individual not configured");

    for (int i = 0; i < config->nreal; ++i) {
        xreal[i] = rgen.real(config->limits_realvar[i].first,
			     config->limits_realvar[i].second);
    }

    for (int i = 0; i < config->nbin; ++i) {
        for (int j = 0; j < config->nbits[i]; ++j) {
            gene[i][j] = rgen.realu() <= 0.5 ? 0 : 1;
        }
    }
}

void individual::decode() {
    int sum;
    for (int i = 0; i < config->nbin; ++i) {
        sum = 0;
        for (int j = 0; j < config->nbits[i]; ++j) {
            sum += (1 << (config->nbits[i]-1-j));  // TODO: check
        }

        xbin[i] = config->limits_binvar[i].first +
            (double)sum*( config->limits_binvar[i].second - config->limits_binvar[i].first) / (double)((1 << (config->nbits[i]))-1); // TODO: check
    }
}

void individual::evaluate() {

    // workaround to respect the signature of test_problem and its (int**)
    vector<int*> tmp(gene.size());
    for (unsigned i=0; i < gene.size(); ++i) {
        tmp[i] = &(gene[i][0]);
    }

    (*config->function) (&xreal[0], &xbin[0], &tmp[0], &obj[0], &constr[0]);

    if (config->ncon) {
      constr_violation = 0.0;
      for (int i = 0; i < config->ncon; ++i)
        if (constr[i] < 0.0)
          constr_violation += constr[i];
    } else {
      constr_violation = 0.0;
    }

    evaluated = true;
}

// returns:  1 if this < b (this dominates b),
//          -1 if this > b (this is dominated by b),
//           0 if they are nondominated
int individual::check_dominance(const individual& b) const {

    if (constr_violation < 0 && b.constr_violation < 0) {
        // both have constraint violations

        if (constr_violation > b.constr_violation)
            return 1; // this violates less
        else if (constr_violation < b.constr_violation)
            return -1; // b violates less
        else
            return 0; // they both violate equally

    } else if (constr_violation < 0 && b.constr_violation == 0) {
        // this violates and b doesn't => b dominates

        return -1;

    } else if (constr_violation == 0 && b.constr_violation < 0) {
        // this doesn't violate and b does => this dominates

        return 1;

    } else {
        // no constraint violations

        int flag1 = 0, // to check if this has a smaller objective
            flag2 = 0; // to check if b    has a smaller objective

        for (int i=0; i<config->nobj; ++i) {
	    if (config->nobj > 1) { // Normal multi objective comparison
		if (obj[i] < b.obj[i]) {
		    flag1 = 1;
		} else if (obj[i] > b.obj[i]) {
		    flag2 = 1;
		}
	    } else { // mono objective comparison with an epsilon
		if (obj[i] < b.obj[i] && fabs(obj[i]-b.obj[i]) > config->epsilon_c) {
		    flag1 = 1;
		} else if (obj[i] > b.obj[i] && fabs(obj[i]-b.obj[i]) > config->epsilon_c) {
		    flag2 = 1;
		}
	    }
        }

        if (flag1==1 && flag2==0) {
            // there is at least one smaller objective for this and none for b

            return 1;

        } else if (flag1==0 && flag2==1) {
            // there is at least one smaller objective for b and none for this

            return -1;

        } else {
            // no smaller objective or both have one smaller

            return 0;
        }

    }
}

// returns num_mut_real, num_mut_bin
pair<int,int> individual::mutate() {
    pair<int,int> num_mut = make_pair(0,0);
    if (config->nreal)
        num_mut.first  += real_mutate();
    if (config->nbin)
        num_mut.second += bin_mutate();
    return num_mut;
}

int individual::real_mutate() {
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    int num_mut = 0;
    for (j=0; j<config->nreal; j++) {
        if (rgen.realu() <= config->pmut_real) {
            y = xreal[j];
            yl = config->limits_realvar[j].first;
            yu = config->limits_realvar[j].second;
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rgen.realu();
            mut_pow = 1.0/(config->eta_m+1.0);
            if (rnd <= 0.5) {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(config->eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            } else {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(config->eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            xreal[j] = y;
            num_mut+=1;
        }
    }
    return num_mut;
}

int individual::bin_mutate() {
    int j, k;
    double prob;
    int num_mut = 0;
    for (j=0; j<config->nbin; j++) {
        for (k=0; k<config->nbits[j]; k++) {
            prob = rgen.realu();
            if (prob <=config->pmut_bin) {
                if (gene[j][k] == 0) {
                    gene[j][k] = 1;
                } else {
                    gene[j][k] = 0;
                }
                num_mut+=1;
            }
        }
    }
    return num_mut;
}


ostream& mou::operator<< (ostream& os, const individual& ind) {

    os << "{Individual rank=" << ind.rank
       << "\nconstr_violation=" << ind.constr_violation;

    os << "\nxreal=[";
    vector<double>::const_iterator it;
    for (it = ind.xreal.begin(); it != ind.xreal.end(); ++it) {
        os << *it;
        if (it+1 != ind.xreal.end())
            os << ",";
    }

    os << "]\ngene=";
    vector< vector<int> >::const_iterator it1;
    for (it1 = ind.gene.begin(); it1 != ind.gene.end(); ++it1) {
        const vector<int>& tmp = *it1;
        vector<int>::const_iterator it2;
        if (it1 != ind.gene.begin())
            os << "     "; // tab space
        for (it2 = tmp.begin(); it2 != tmp.end(); ++it2) {
            os << *it2;
        }
        //       gene=
        os << '\n';
    }

    os << "xbin=";
    for (it = ind.xbin.begin(); it != ind.xbin.end(); ++it) {
        os << *it;
        if (it+1 != ind.xbin.end())
            os << ",";
    }

    os << "\nobj=";
    for (it = ind.obj.begin(); it != ind.obj.end(); ++it) {
        os << *it;
        if (it+1 != ind.obj.end())
            os << ",";
    }

    os << "\nconstr=";
    for (it = ind.constr.begin(); it != ind.constr.end(); ++it) {
        os << *it;
        if (it+1 != ind.constr.end())
            os << ",";
    }

    os << "\ncrowd_dist=" << ind.crowd_dist;

    os << " }";

    return os;
}

population::population(const int size,
                       const int nreal,
                       const int nbin,
                       const int ncon,
                       const vector<int>& nbits,
                       const vector< pair<double,double> >& limreal,
                       const vector< pair<double,double> >& limbin,
                       const int nobj,
                       const double pmut_real,
                       const double pmut_bin,
                       const double eta_m,
		       const double epsilon_c,
                       const individual_config::funcType func)
	  throw (mou::MOEAexception) :
	  crowd_obj(true),
	  ind_config(),
	  eval_pop_function(NULL) {

    generation = 1;
    ind_config.nreal          = nreal;
    ind_config.nbin           = nbin;
    ind_config.nobj           = nobj;
    ind_config.ncon           = ncon;
    ind_config.nbits          = nbits;
    ind_config.limits_realvar = limreal;
    ind_config.limits_binvar  = limbin;
    ind_config.pmut_real      = pmut_real;
    ind_config.pmut_bin       = pmut_bin;
    ind_config.eta_m          = eta_m;
    ind_config.function       = func;
    ind_config.epsilon_c      = epsilon_c;

    for (int i = 0; i < size; ++i) {
        ind.push_back(individual(ind_config));
    }

}

population::~population() {
}

void population::initialize() throw (mou::MOEAexception) {
    vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->initialize();
    }
}

void population::decode() {
    vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->decode();
    }
}

void population::evaluate() {
    normal_evaluate_openmp();
}

void population::custom_evaluate() {
    if (eval_pop_function != NULL)
	(*eval_pop_function)(*this);
    else 
	normal_evaluate_openmp();
}

void population::normal_evaluate() {
    vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->evaluate();
    }
}

void population::normal_evaluate_openmp() {
#ifdef USE_OPENMP
#pragma omp parallel for
    for (int i = 0; i < ind.size(); ++i) {
        ind[i].evaluate();
    }
#else
    normal_evaluate();
#endif
}

void population::set_popfunction(individual_config::popFuncType f) {
    eval_pop_function = f;
}

void population::fast_nds() {
    front.resize(1);
    front[0].clear();
    //std::vector< std::vector<int> >  F(1);
#pragma omp parallel for
    for (int i = 0; i < ind.size(); ++i) {
	
        vector<int> dom;
        int dcount = 0;
	
        individual& p = ind[i];
        // p.dcounter  = 0;
        // p.dominated.clear();
	
        for (int j = 0; j < ind.size(); ++j) {
	    
            individual& q = ind[j];
	    
            int compare = p.check_dominance(q);
            if (compare == 1) { // p dominates q
                //p.dominated.push_back(j);
                dom.push_back(j);
            } else if (compare == -1) { // q dominates p
                //p.dcounter += 1;
                dcount += 1;
            }
        }
	
#pragma omp critical
        {
            p.dcounter  = dcount;
            p.dominated.clear();
            p.dominated = dom;
	    
	    
            if (p.dcounter == 0) {
                p.rank = 1;
                front[0].push_back(i);
            }
        }
	
    }
    
    // using OpenMP can have different orders in the front[0]
    // so let's sort it so that the algorithm is deterministic
    // given a seed
    sort(front[0].begin(), front[0].end());    

    int fi = 1;
    while (front[fi-1].size() > 0) {

        vector<int>& fronti = front[fi-1];
        vector<int> Q;
        for (int i = 0; i < fronti.size(); ++i) {

            individual& p = ind[fronti[i]];

            for (int j = 0; j < p.dominated.size() ; ++j) {


                individual& q = ind[p.dominated[j]];
                q.dcounter -= 1;

                if (q.dcounter == 0) {
                    q.rank = fi+1;
                    Q.push_back(p.dominated[j]);
                }
            }
        }


        fi += 1;
        front.push_back(Q);
    }

}

struct comparator_obj {
    comparator_obj(const population& population, int index) :
        pop(population), m(index) {};
    const population& pop;
    int m;
    bool operator() (int i, int j) {
        return pop.crowd_obj?
              pop.ind[i].obj[m] < pop.ind[j].obj[m]
            : pop.ind[i].xreal[m] < pop.ind[j].xreal[m];
    };
};

void population::crowding_distance_all() {
    for (int i = 0; i < front.size(); ++i)
        crowding_distance(i);
}

void population::crowding_distance(int fronti) {

    vector<int> F = front[fronti];
    if (F.size() == 0 ) return;

    const int l = F.size();

    for (int i = 0; i < l; ++i)
        ind[F[i]].crowd_dist = 0;

    // for (int m = 0; m < ind_config.nobj; ++m) {
    //     std::sort(F.begin(), F.end(), comparator_obj(*this,m));
    //     ind[F[0]].crowd_dist = INF;
    // }

    const int limit = crowd_obj?ind_config.nobj:ind_config.nreal;
    for (int m = 0; m < limit; ++m) {

        sort(F.begin(), F.end(), comparator_obj(*this,m));

        // in the paper dist=INF for the first and last, in the code
        // this is only done to the first one or to the two first when size=2
        // ind[F[0]].crowd_dist = INF;
        // if (l == 2)
        //      ind[F[0]].crowd_dist = ind[F[1]].crowd_dist = INF;
        ind[F[0]].crowd_dist = INF;
        if (l > 1)
            ind[F[l-1]].crowd_dist = INF;
	cout << "min " << ind[F[0]].xreal[0];
	cout << "\tmax " << ind[F[l-1]].xreal[0] << endl;

        for (int i = 1; i < l-1; ++i) {
            if (ind[F[i]].crowd_dist != INF) {
                if (crowd_obj && ind[F[l-1]].obj[m] != ind[F[0]].obj[m]) {
                // ind[F[l-1]].xreal[m] != ind[F[0]].xreal[m])
                    ind[F[i]].crowd_dist +=
                        (ind[F[i+1]].obj[m] - ind[F[i-1]].obj[m]) // crowd over obj
                        / (ind[F[l-1]].obj[m] - ind[F[0]].obj[m]);
                 } else if (!crowd_obj && ind[F[l-1]].xreal[m] != ind[F[0]].xreal[m]) {
                    ind[F[i]].crowd_dist +=
                        (ind[F[i+1]].xreal[m] - ind[F[i-1]].xreal[m]) // crowd over vars
                        / (ind[F[l-1]].xreal[m] - ind[F[0]].xreal[m]);
                }
            }
        }
    }

    // for (int i=0; i < l; ++i) { // this is deduced from code, not mentioned in paper
    //     if (ind[F[i]].crowd_dist != INF)
    //         ind[F[i]].crowd_dist /= ind_config.nobj;
    // }
}

void population::merge(const population& pop1, const population& pop2)
    throw (mou::MOEAexception) {

    if (size() < pop1.size() + pop2.size())
        throw mou::MOEAexception("Merge: target population not big enough");

    copy(pop1.ind.begin(), pop1.ind.end(), ind.begin());
    copy(pop2.ind.begin(), pop2.ind.end(), ind.begin() + pop1.size());

}

void population::report(ostream& os) const {

    vector<individual>::const_iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        for (int j = 0; j < ind_config.nobj; ++j)
            os << it->obj[j] << '\t';
        for (int j = 0; j < ind_config.ncon; ++j)
            os << it->constr[j] << '\t';
        for (int j = 0; j < ind_config.nreal; ++j)
            os << it->xreal[j] << '\t';
        for (int j = 0; j < ind_config.nbin; ++j)
            for (int k = 0; k < ind_config.nbits[j]; ++k)
                os << it->gene[j][k] << '\t';

        os << it->constr_violation << '\t'
           << it->rank << '\t'
           << it->crowd_dist << '\n';

    }
}

void population::dump(ostream& os) const {

    vector<individual>::const_iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        if (ind_config.nobj > 0)
            os.write(reinterpret_cast<const char*>(&(it->obj[0])),
                     sizeof(double)*ind_config.nobj);

        if (ind_config.ncon > 0)
            os.write(reinterpret_cast<const char*>(&(it->constr[0])),
                     sizeof(double)*ind_config.ncon);

        if (ind_config.nreal > 0)
            os.write(reinterpret_cast<const char*>(&(it->xreal[0])),
                     sizeof(double)*ind_config.nreal);

        for (int j = 0; j < ind_config.nbin; ++j)
            os.write(reinterpret_cast<const char*>(&(it->gene[j][0])),
                     sizeof(int)*ind_config.nbits[j]);

        os.write(reinterpret_cast<const char*>(&(it->constr_violation)),
                 sizeof(double));
        os.write(reinterpret_cast<const char*>(&(it->rank)),
                 sizeof(int));
        os.write(reinterpret_cast<const char*>(&(it->crowd_dist)),
                 sizeof(double));

    }
}

void population::load(istream& os) {

    vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        if (ind_config.nobj > 0)
            os.read(reinterpret_cast<char*>(&(it->obj[0])),
                     sizeof(double)*ind_config.nobj);

        if (ind_config.ncon > 0)
            os.read(reinterpret_cast<char*>(&(it->constr[0])),
                     sizeof(double)*ind_config.ncon);

        if (ind_config.nreal > 0)
            os.read(reinterpret_cast<char*>(&(it->xreal[0])),
                     sizeof(double)*ind_config.nreal);

        for (int j = 0; j < ind_config.nbin; ++j)
            os.read(reinterpret_cast<char*>(&(it->gene[j][0])),
                     sizeof(int)*ind_config.nbits[j]);

        os.read(reinterpret_cast<char*>(&(it->constr_violation)),
                 sizeof(double));
        os.read(reinterpret_cast<char*>(&(it->rank)),
                 sizeof(int));
        os.read(reinterpret_cast<char*>(&(it->crowd_dist)),
                 sizeof(double));

    }
}

pair<int,int> population::mutate() {
    pair<int,int>
        num_mut = make_pair(0,0),
        tmp = make_pair(0,0);
    vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        tmp = it->mutate();
        num_mut.first  += tmp.first;
        num_mut.second += tmp.second;
    }
    return num_mut;
}

ostream& mou::operator<< (ostream& os, const population& pop) {
    os << "Population: {\n";
    vector<individual>::const_iterator it;
    for (it = pop.ind.begin(); it != pop.ind.end(); ++it) {
        os << *it;
    }
    os << '}';
    return os;
}
