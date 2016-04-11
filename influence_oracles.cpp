#include <igraph.h>
#include <vector>

using namespace std;

class influence_oracles {
public:
  unsigned ell;
  unsigned k;

  vector< igraph_t* > v_instances;

  vector< vector< unsigned > > v_sketches;

  void compute_oracles() {}
  
  void estimate_influence() {}

  void merge_sketches(vector< unsigned >& sketch_1, vector< unsigned >& sketch_2, 
		      vector< unsigned >& result ) {}

  influence_oracles( vector< igraph_t* > in_instances, unsigned in_ell, unsigned in_k ) :
    v_instances( in_instances ), ell( in_ell ), k( in_k ) { }

};
