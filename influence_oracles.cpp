#include <igraph.h>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>

using namespace std;
typedef unsigned long myint;
typedef unordered_set< myint >  uset;

struct mypair {
  myint first;
  myint second;
};

class influence_oracles {
public:
  myint ell;
  myint k;
  myint n;

  vector< igraph_t* > v_instances;

  vector< vector< myint > > v_sketches;

  void compute_oracles();
  
  void estimate_influence() {}

  void merge_sketches(vector< myint >& sketch_1, vector< myint >& sketch_2, 
		      vector< myint >& result ) {}

  influence_oracles( vector< igraph_t* > in_instances, myint in_ell, myint in_k ) :
    v_instances( in_instances ), ell( in_ell ), k( in_k ) { }

};

void influence_oracles::compute_oracles() {
  // create the permutation of size n*ell
  myint K = n * ell;
  vector< mypair > perm;
  perm.reserve( K );
  for (myint u = 0; u < n; ++u) {
    for (myint i = 0; i < ell; ++i) {

      mypair tmp;
      tmp.first = u;
      tmp.second = i;
      perm.push_back( tmp );
    }
  }

  random_shuffle( perm.begin(), perm.end() );

  // group ranks by instance
  vector< vector < mypair > > instanceRanks;
  for (myint r = 1; r <= K; ++r ) {
    myint i = perm[ r - 1 ].second;
    mypair tmp;
    tmp.first = r;
    tmp.second = perm[ r - 1 ].first;
    instanceRanks[i].push_back( tmp );
  }

  // compute combined bottom-k rank sketches
  
  // for each instance i
  for (myint i = 0; i < ell; ++i) {
    // for each vertex u in i by increasing rank
    for (myint j = 0; j < n; ++j) {
      myint vertex = instanceRanks[i].second;
      // Run BFS in instance i from u until sketch is of
  // size k, or BFS ends

  // merge the reachability sketch in instance i
  // into the global sketch. 

}
