#include <string>
#include <unordered_map>
#include <set>
//#include <utility>
#include <iostream>
//#include "hashtable_wrapper.h"
#include "utilities.h"
//#include "pest_error.h"

using namespace std;
using namespace pest_utils;


typedef class RunManagerAbstract RunManagerAbstract;

set<unordered_map<string, int>* > hashtable_set;


extern "C"
{
  void hash_add_table_(char *names, int  *names_str_len, int *names_array_len,
		    void **ptr)

  {
	vector<string> names_vec =  fortran_str_array_2_vec(names, *names_str_len, *names_array_len);

	unordered_map<string, int> *new_table =  new  unordered_map<string, int>;
	::hashtable_set.insert(new_table);
	int i = 0;
	int n_names = names_vec.size();
	for(int i=0; i<n_names; ++i)
	{
	  (*new_table).insert(std::make_pair(names_vec[i], i+1));
	}
        *ptr = new_table;
  }

  void  hash_which1_cpp(void **ptr, string &name, int *index, int *ifail)
  {
    *ifail = 1;
    *index = 0;
    unordered_map<string, int> *table = static_cast< unordered_map<string, int> *>(*ptr);
    set<unordered_map<string, int> *>::iterator it = hashtable_set.find(table);
    if (it == hashtable_set.end())
    {
      *ifail = 1;
    }
    else
    {
      unordered_map<string, int>::iterator it_name = (*table).find(name);
      if (it_name != (*table).end())
      {
        *index = (*it_name).second;
        *ifail = 0;
      }
      else
      {
        *ifail = 1;
      }
    }
  }

  void  hash_which1_(void **ptr, char *f_name, int *f_name_len, int *index, int *ifail)
  {
    string name =  fortran_str_2_string(f_name, *f_name_len);
    hash_which1_cpp(ptr, name, index, ifail);
  }

  void  hash_which1_lower_(void **ptr, char *f_name, int *f_name_len, int *index, int *ifail)
  {
    string name =  fortran_str_2_string(f_name, *f_name_len);
    lower_ip(name);
    hash_which1_cpp(ptr, name, index, ifail);
  }


  void hash_free_(void **ptr)
  {
    unordered_map<string, int> *table = static_cast< unordered_map<string, int> *>(*ptr);
    set<unordered_map<string, int> *>::iterator it = hashtable_set.find(table);
    if (it == hashtable_set.end())
    {
      cerr << "hash_free send illegal pointer" << endl;
    }
    else{
      hashtable_set.erase(*it);
      delete (table);
    }

  }
}
