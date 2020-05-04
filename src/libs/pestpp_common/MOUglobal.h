#ifndef MOUGLOBAL_H_
#define MOUGLOBAL_H_

#include <vector>
#include <ostream>
#include <string>
#include <utility>

#define EPS 1e-14
#define INF 1e+14

namespace mou {
    class MOEAexception : public std::exception {
    public:
        explicit MOEAexception(const std::string& msg) :
            std::exception(),
            message(msg) {};
        virtual ~MOEAexception() throw () {};

        const char* what() const throw() {
            return message.c_str();
        };

    private:
        std::string message;
    };
    class NSGA2;
    struct population;

    struct individual_config {
    
        typedef void (*funcType)(double *,
                                 double *,
                                 int **,
                                 double *,
                                 double *);
        typedef void (*popFuncType)(population&);
    
        int nreal;
        int nbin;
        int nobj;
        int ncon;
        // int popsize;
        // int ngen;
        // double pcross_real;
        // double pcross_bin;
        double pmut_real;
        double pmut_bin;
        // double eta_c;
        double eta_m;
        std::vector<int> nbits;
        std::vector< std::pair<double,double> > limits_realvar;
        std::vector< std::pair<double,double> > limits_binvar;
        funcType function;
        double epsilon_c;
    };

struct individual {

    individual() throw (); // needed for std::vector<individual> allocator

    individual(const individual_config& c) throw (mou::MOEAexception);
    virtual ~individual();

    // individual& operator=(const individual& ind);

    void initialize() throw (mou::MOEAexception);

    void decode();
    void evaluate();
    std::pair<int,int> mutate();
    int real_mutate();
    int bin_mutate();

    int check_dominance(const individual& b) const;

    int rank;
    double constr_violation;
    std::vector<double> xreal;
    std::vector< std::vector<int> > gene;
    std::vector<double> xbin;
    std::vector<double> obj;
    std::vector<double> constr;
    double crowd_dist;

    int dcounter; // domination counter n_p
    std::vector<int> dominated;
    bool evaluated;

private:
    const individual_config* config;
    friend std::ostream& operator<< (std::ostream& os, const individual& ind);
};

std::ostream& operator<< (std::ostream& os, const individual& ind);

struct population {

    population(const int size,
               const int nreal,
               const int nbin,
               const int ncon,
               const std::vector<int>& nbits,
               const std::vector< std::pair<double,double> >& limreal,
               const std::vector< std::pair<double,double> >& limbin,
               const int nobj,
               const double pmut_real,
               const double pmut_bin,
               const double eta_m,
	       const double epsilon_c,
               const individual_config::funcType func)
        throw (mou::MOEAexception);
    virtual ~population();

    void initialize() throw (mou::MOEAexception);

    void decode();
    void evaluate();
    void custom_evaluate(); // this one takes into account custom evaluation by user
    void set_popfunction(individual_config::popFuncType f);
    void fast_nds();
    void crowding_distance_all();
    void crowding_distance(int fronti);

    std::pair<int,int> mutate();
    void merge(const population& pop1, const population& pop2)
    throw (mou::MOEAexception);

    void report(std::ostream& os) const;
    void dump(std::ostream& os) const;
    void load(std::istream& os);

    int size() const {
        return ind.size();
    };

    std::vector<individual> ind;
    std::vector< std::vector<int > > front;

    bool crowd_obj; // true: crowding over objective (default) false: crowding over real vars
    int generation;

private:
    individual_config ind_config;
    individual_config::popFuncType eval_pop_function;

    void normal_evaluate();
    void normal_evaluate_openmp();


    friend std::ostream& operator<< (std::ostream& os, const population& pop);
};

std::ostream& operator<< (std::ostream& os, const population& pop);
}

#endif /* MOUGLOBAL_H_ */
