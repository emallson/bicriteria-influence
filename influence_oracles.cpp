#include <igraph.h>
#include <vector>
#include <unordered_set>

using namespace std;
typedef unsigned long myint;
typedef unordered_set< myint >  uset;


class influence_oracles {
public:
  myint ell;
  myint k;

  vector< igraph_t* > v_instances;

  vector< vector< myint > > v_sketches;

  void compute_oracles();
  
  void estimate_influence() {}

  void merge_sketches(vector< myint >& sketch_1, vector< myint >& sketch_2, 
		      vector< myint >& result ) {}p

  influence_oracles( vector< igraph_t* > in_instances, myint in_ell, myint in_k ) :
    v_instances( in_instances ), ell( in_ell ), k( in_k ) { }

};

void influence_oracles::compute_oracles() {
  // create the permutation of size n*ell

  // group ranks by instance

  // compute combined bottom-k rank sketches

  // for each instance i

  // for each vertex u in i

  // Run BFS in instance i from u until sketch is of
  // size k, or BFS ends

  // merge the reachability sketch in instance i
  // into the global sketch. 

}
