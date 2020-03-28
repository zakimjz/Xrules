#ifndef __treeminer_H
#define __treeminer_H
 
#include <vector>
#include <functional>
using namespace std;

#define BranchIt -1 //-1 indicates a branch in the tree
#define Max_Num_Classes 25 //max classes by default
//enums
enum sort_vals {nosort, incr, decr};
enum alg_vals {treeminer, maxtreeminer};
enum prune_vals {noprune, prune};
  
//externs
extern vector<double> MINSUP_PER;
extern vector<int> MINSUPPORT;
extern int DBASE_MAXITEM;
extern int DBASE_NUM_TRANS;
extern vector<int> DBASE_NUMCL_TRANS;

extern bool output;
extern bool output_idlist;
extern bool use_class_ratio;
extern bool print_ratio;

extern sort_vals sort_type;  
extern alg_vals alg_type;
extern prune_vals prune_type;
extern int Num_Classes;

template<class T>
struct delnode: public unary_function<T, void>{
   void operator() (T x){ delete x; }
};


#endif
