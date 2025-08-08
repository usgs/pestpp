/*


	This file is part of PEST++.

	PEST++ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	PEST++ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_


/** @file
 @brief Transformation Classes

 This file defines the classes used to implement transformations.
*/
#include <string>
#include <map>
#include <set>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Transformable.h"
#include "SVD_PROPACK.h"
//#include "covariance.h"

class Jacobian;
class QSqrtMatrix;
class ParameterGroupInfo;
class Parameters;
class PerformanceLog;

using namespace std;

/**
 @brief Transformation Base Class

 This is the pure virtual base class for all transformations.  The class can not be
 instantiated and all transformations classes must defive from it.
*/
class Transformation {
public:
	Transformation(const string &_name="unnamed Transformation") : name(_name) {}
	Transformation(const Transformation &rhs) : name(rhs.name) {}
	virtual ~Transformation(){};
	 /** Perform a forward transformation on the dataset stored in data.  The
	 dataset contained in data is changed in place.
	 */
	virtual void forward(Transformable &data) = 0;
	 /** Perform a reverse transformation on the dataset stored in data.  The
	 dataset contained in data is changed in place.
	 */
	virtual void reverse(Transformable &data) = 0;

	virtual void jacobian_forward(Jacobian &jac) = 0;
	virtual void jacobian_reverse(Jacobian &jac) = 0;

	virtual void d1_to_d2(Transformable &del_data, Transformable &data) = 0;
	virtual void d2_to_d1(Transformable &del_data, Transformable &data) = 0;

	/** Returns the name of the transformation.
	 */
	string get_name()const {return name;}
	/** Prints the information associated with the transformation.
	This is used primarily for debugging purposes
	 */
	virtual void print(ostream &os) const = 0;
	/** If very member of a transformed dataset maps to single and independent member
	under both the forward and reverse transformation this method will return true.  Otherwise
	it returns false.
	 */
	virtual bool is_one_to_one() const {return false;}
	/** Returns a pointer to a copy of current instance of this class
	 */
	virtual Transformation* clone() const= 0;
protected:
	string name;
};

/**
 @brief Transformation Base Class for Transformations based on an STL map

 This is the pure virtual base class for Transformations based on an STL map.  The class can not be
 instantiated but provides a common framework to build transformations based using an STL map
*/
class TranMapBase: public Transformation {
public:
	TranMapBase(const string &_name="unnamed TranMapBase"): Transformation(_name){};
	TranMapBase(const TranMapBase &rhs) : Transformation(rhs), items(rhs.items) {}
	virtual ~TranMapBase(){};
	virtual void forward(Transformable &data)=0;
	virtual void reverse(Transformable &data)=0;
	/** Adds a transformation item to the items member map.  _name is the name of
	   the parameter to be transformed and _value is the value to be applied in the forward
	   transformtion.  How _value is applied in the transformation is determined in the
	   child classes
	 */
	void insert(const string &item_name, double item_value);
	void clear() {items.clear();}
	void insert (const Parameters &pars);
	void reset(const Parameters &pars);
	/** Returns the transformation value associated with the name of an transformable item.
	The results are returned as a pair<bool, double> where the first element in the pair will be
	false if the items is not part of the transformation and true if it is
	 */
	pair<bool, double> get_value(const string &name) const;
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return false;}
	virtual TranMapBase* clone() const = 0;

protected:
	map<string, double> items;
};

/**
 @brief Transformation Base Class for Transformations based on an STL set

 This is the pure virtual base class for Transformations based on an STL set.  The class can not be
 instantiated but provides a common framework to build transformations based using an STL set
*/
class TranSetBase: public Transformation {
public:
	TranSetBase(const string &_name="unnamed TranSetBase"): Transformation(_name){};
	TranSetBase(const TranSetBase &rhs) : Transformation(rhs), items(rhs.items) {}
	virtual void forward(Transformable &data) = 0;
	virtual void reverse(Transformable &data) = 0;
	void insert(const string &data_name);
	virtual bool has_value(const string &name) const;
	virtual ~TranSetBase(){}
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return false;}
	virtual TranSetBase* clone() const = 0;
	virtual const set<string>& get_items() const {return items;}
protected:
	set<string> items;
};


/**
 @brief Offset Transformation Class

 This class provides an offset transformation in which an offset is added in the forward transformation
 and subtracted in the reverse transformation.
*/
class TranOffset: public TranMapBase {
public:
	TranOffset(const string &_name="unnamed TranOffset"): TranMapBase(_name){};
	TranOffset(const TranOffset &rhs) : TranMapBase(rhs) {}
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	virtual ~TranOffset(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranOffset* clone() const {return new TranOffset(*this);}
private:
};

/**
 @brief Scale Transformation Class

 This class provides an scaled transformation in which value is multiple by scaling value in the forward transformation
 and dvided by it in the reverse transformation.
*/
class TranScale: public TranMapBase {
public:
	TranScale(const string &_name="unnamed TranScale"): TranMapBase(_name){};
	TranScale(const TranScale &rhs) : TranMapBase(rhs) {}
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	virtual ~TranScale(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranScale* clone() const {return new TranScale(*this);}
private:
};


/**
 @brief Log10 Transformation Class

 This class provides an log10 transformation in which value is log10 transformed is applied in the forward direction
 and inverse log10 transform is applied in the reverse direction.
*/
class TranLog10: public TranSetBase {
public:
	TranLog10(const string &_name="unnamed TranLog10"): TranSetBase(_name){};
	TranLog10(const TranLog10 &rhs) : TranSetBase(rhs) {}
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	virtual ~TranLog10(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranLog10* clone() const {return new TranLog10(*this);}
};

/**
 @brief Fixed Value Transformation Class

 This class provides transformation that adds an additional transformable items with a fixed
 value in the forward transformation and removes the item in the reverse transformation.
*/
class TranFixed: public TranMapBase {
public:
	TranFixed(const string &_name="unknown TranFixed"): TranMapBase(_name){};
	TranFixed(const TranFixed &rhs) : TranMapBase(rhs) {}
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	virtual ~TranFixed(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranFixed* clone() const {return new TranFixed(*this);}
private:
};

class TranFrozen: public TranFixed {
public:
	TranFrozen(const string &_name="unknown TranFrozen"): TranFixed(_name){};
	TranFrozen(const TranFrozen &rhs) : TranFixed(rhs) {}
	virtual ~TranFrozen(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranFrozen* clone() const {return new TranFrozen(*this);}
private:
};


/**
 @brief Tied Transformation Class

 This class provides transformation that ties a transformable item to another transformable
 item.  In the forard direction this transformation will add a new Tied transformable item.
 This item is assigned the scaled value of the items that it is tied to.  In the reverse direction
 the tied parameters are simple removed from the transformable set.
*/
class TranTied: public Transformation {
public:
	typedef pair<string, double> pair_string_double;
	TranTied(const string &_name="unknown TranTied"): Transformation(_name){};
	TranTied(const TranTied &rhs) : Transformation(rhs), items(rhs.items) {}
	void insert(const string &item_name, const pair<string, double> &item_value);
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	virtual ~TranTied(){};
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return true;}
	virtual TranTied* clone() const {return new TranTied(*this);}
	const map<string, pair_string_double> get_items() { return items; }
protected:
	map<string, pair_string_double> items;
};

/**
 @brief Super Parameter or SVD Assist (SVDA) Transformation

 This class provides a transformation for super parameters.  In the forward direction base parameters
 are transformed to super parameters and in the reverse direction super parameters are transformed to
 base parameters
*/
class TranSVD: public Transformation {
public:
	TranSVD(int _max_sing, double _eign_thresh, const string &_name = "unnamed TranSVD");
	TranSVD(const TranSVD &rhs);
	void set_SVD_pack();
	void update_reset_frozen_pars(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const Parameters &base_numeric_pars,
		int maxsing, double eigthresh, const vector<string> &par_names, const vector<string> &obs_names, 
		Eigen::SparseMatrix<double>& parcov_inv, map<string,double> dss,const Parameters &_frozen_derivative_pars=Parameters());
	void update_add_frozen_pars(const Parameters &_frozen_derivative_pars);
	Parameters& get_frozen_derivative_pars() {return frozen_derivative_parameters;}
	const vector<string>& get_super_parameter_names(){return super_parameter_names;}
	virtual void forward(Transformable &data);
	virtual void reverse(Transformable &data);
	virtual void jacobian_forward(Jacobian &jac);
	virtual void jacobian_reverse(Jacobian &jac);
	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
	void insert(const string &item_name, double item_value)
	    {/*can't insert into SVD.  This is intentionally a no-op*/};
	virtual ~TranSVD();
	virtual void print(ostream &os) const;
	virtual bool is_one_to_one() const {return false;}
	virtual TranSVD* clone() const {return new TranSVD(*this);}
	ParameterGroupInfo build_par_group_info(const ParameterGroupInfo &base_pg_info);
	Parameters map_basepar_to_super(const Parameters &base_pars);
	const Eigen::SparseMatrix<double>& get_vt() const;
	void save(ostream &fout) const;
	void read(istream &fin);
	void set_performance_log(PerformanceLog *_performance_log);
protected:
	SVDPackage *tran_svd_pack;

	vector<string> base_parameter_names;
	vector<string> super_parameter_names;
	vector<string> obs_names;
	//Eigen::SparseMatrix<double> SqrtQ_J;
	Eigen::SparseMatrix<double> jtqj;
	Eigen::VectorXd Sigma;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> Vt;
	Parameters init_base_numeric_parameters;
	Parameters frozen_derivative_parameters;
	void calc_svd();
	map<string, double> dss;
};

//class TranNormalize: public Transformation {
//public:
//	class NormData {
//	public:
//		double offset;
//		double scale;
//		NormData(double _offset=0.0, double _scale=1.0) :offset(_offset), scale(_scale){}
//		~NormData(){};
//	};
//	TranNormalize(const string &_name="unknown TranNormalize"): Transformation(_name){};
//	TranNormalize(const TranNormalize &rhs) : Transformation(rhs), items(rhs.items) {}
//	void insert(const string &item_name, double _offset, double _scale);
//	virtual void forward(Transformable &data);
//	virtual void reverse(Transformable &data);
//	virtual void jacobian_forward(Jacobian &jac);
//	virtual void jacobian_reverse(Jacobian &jac);
//	virtual void d1_to_d2(Transformable &del_data, Transformable &data);
//	virtual void d2_to_d1(Transformable &del_data, Transformable &data);
//	virtual ~TranNormalize(){};
//	virtual void print(ostream &os) const;
//	virtual bool is_one_to_one() const {return true;}
//	virtual TranNormalize* clone() const {return new TranNormalize(*this);}
//protected:
//	map<string, NormData> items;
//};
#endif /* TRANSFORMATION_H_ */
