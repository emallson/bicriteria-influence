#include <igraph.h>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <queue>
#include <iostream>

using namespace std;
typedef igraph_integer_t myint;
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

  vector< vector< myint > > global_sketches;

  void compute_oracles();
  void alt_compute_oracles();

  double estimate_reachability_sketch( vector< myint >& asketch ) {
    double uniform_rank;
    double estimate;
    if (asketch.size() == k) {
      myint T = asketch.back();
      uniform_rank = ((double) T - 1.0)/(ell * n - 1.0 );
      estimate = ((double) k - 1)/ uniform_rank;
    }
    else {
      estimate = ((double) asketch.size());
    }
      
    
    //double estimate = 1.0 + 
    //      ((double)(k - 1)*(n*ell - 1)) / (T - 1);

    return estimate / ell;

  }

  double estimate_reachability( myint vertex ) {
    double uniform_rank;
    double estimate;
    if (global_sketches[ vertex ].size() == k) {
      myint T = global_sketches[ vertex ].back();
      uniform_rank = ((double) T - 1.0)/(ell * n - 1.0 );
      estimate = ((double) k - 1)/ uniform_rank;
    }
    else {
      estimate = ((double) global_sketches[vertex].size());
    }
      




    
    //double estimate = 1.0 + 
    //      ((double)(k - 1)*(n*ell - 1)) / (T - 1);

    return estimate / ell;
  }

  double average_reachability( myint vertex ) {
    myint tot_reach = 0;
    for (myint i = 0; i < ell; ++i) {
      tot_reach += forwardBFS( v_instances[i], vertex );

    }

    return ((double) tot_reach) / ell;
  }

  void merge_sketches(vector< myint >& sketch_1, vector< myint >& sketch_2, 
		      vector< myint >& result ) {}

  influence_oracles( vector< igraph_t* >& in_instances,
                     myint in_ell, 
                     myint in_k,
                     myint in_n) :
    v_instances( in_instances ), 
    ell( in_ell ), 
    k( in_k ),
    n( in_n )
  { }

  influence_oracles( vector< igraph_t >& in_instances,
                     myint in_ell, 
                     myint in_k,
                     myint in_n) :
    ell( in_ell ), 
    k( in_k ),
    n( in_n )
  { 
    igraph_t* nullpoint;
    v_instances.assign( ell, nullpoint);

    for (myint i = 0; i < ell; ++i) {
      v_instances[i] = &(in_instances[i]);

    }
  }

  myint forwardBFS( igraph_t* G_i, myint vertex ) {

    queue <myint> Q;
    vector < int > dist;
    dist.reserve( n );
    for (myint i = 0; i < n; ++i) {
      dist.push_back( -1 ); //infinity
    }

    dist[ vertex ] = 0;
    Q.push( vertex );
    while (!(Q.empty())) {

      myint current = Q.front();
      Q.pop();
      //get forwards neighbors of current
      igraph_vector_t neis;
      igraph_vector_init( &neis, 0 );
      igraph_neighbors( G_i, &neis, current, IGRAPH_OUT );
      for (myint i = 0; i < igraph_vector_size( &neis ); ++i) {
	myint aneigh = VECTOR( neis )[i];
	if (dist[ aneigh ] == -1 ) {
	  dist[ aneigh ] = dist[ current ] + 1;
	  Q.push( aneigh );

	}
      }

      igraph_vector_destroy( &neis );
    } //BFS finished

    myint count = 0;
    for (myint i = 0; i < n; ++i) {
      if (dist[i] != -1) {
	//i is reachable from vertex in G_i
	++count;
      }
    }

    return count;

  }

  void rankorder_BFS( igraph_t* G_i, 
                      igraph_vector_t& vRO ) {

    //get neighborhoods for all nodes
    igraph_vs_t vs;
    igraph_vs_vector( &vs, &vRO );

    igraph_vector_ptr_t res;
    igraph_neighborhood(
			G_i,
			&res,
			vs,
			n,
			IGRAPH_IN );

    

  }

  void update_local_sketches( igraph_t* G_i, 
			      mypair root, 
			      vector< vector< myint > >& local_sketches ) {
    //root is a pair, first element is rank,
    //second is vertex

    //want to insert the rank of root to sketches of backwards reachable nodes
    //BFS backwards
    queue <myint> Q;
    vector < int > dist;
    dist.reserve( n );
    for (myint i = 0; i < n; ++i) {
      dist[i] = -1; //infinity

    }

    dist[ root.second ] = 0;
    Q.push( root.second );
    while (!(Q.empty())) {

      myint current = Q.front();
      Q.pop();
      //get backwards neighbors of current
      igraph_vector_t neis;
      igraph_vector_init( &neis, 0 );
      igraph_neighbors( G_i, &neis, current, IGRAPH_IN );
      for (myint i = 0; i < igraph_vector_size( &neis ); ++i) {
	myint aneigh = VECTOR( neis )[i];
	if (dist[ aneigh ] == -1 ) {
	  dist[ aneigh ] = dist[ current ] + 1;
	  Q.push( aneigh );
	  //update the local sketch of 'aneigh'
	  //with the rank of the root
	  //since we look at ranks in increasing order,
	  //pushing onto the end will maintain the
	  //property that the sketches are sorted
	  if ( local_sketches[ aneigh ].size() < k ) {
	    local_sketches[ aneigh ].push_back( root.first );
	  }

	}
      }
      igraph_vector_destroy( &neis );
    } //BFS finished

  }

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

  cerr << "Permutation shuffled" << endl;

  // group ranks by instance
  // the pairs are in order of rank
  vector< mypair > emptyInstance;
  vector< vector < mypair > > instanceRanks( ell, emptyInstance );

  for (myint r = 1; r <= K; ++r ) {
    myint i = perm[ r - 1 ].second;
    mypair tmp;
    tmp.first = r;
    tmp.second = perm[ r - 1 ].first;
    instanceRanks[i].push_back( tmp );
  }

  // compute combined bottom-k rank sketches
  
  // for each instance i, compute local bottom-k sketches
  vector < vector< myint > > local_sketches;
  vector < myint > empty_sketch;
  global_sketches.assign( n, empty_sketch );

  for (myint i = 0; i < ell; ++i) {
    cerr << i << endl;
    local_sketches.clear();
    local_sketches.assign( n, empty_sketch );
      
    // for each vertex in instance i by increasing rank
    cerr << "Updating local ranks..." << endl;
    for (myint j = 0; j < n; ++j) {
      myint vertex = instanceRanks[i][j].second;

      // Run reverse BFS in instance i from 'vertex', 
      // updating sketches of discovered vertices
      // (and not including vertex itself)
      update_local_sketches( v_instances[ i ],
			     instanceRanks[i][j],
			     local_sketches );
    }
    
    cerr << "Performing merges..." << endl;
    // merge local_sketches in instance i
    // into the global sketch for each node
    vector< myint > new_sketch;
    for (myint u = 0; u < n; ++u) {
      new_sketch.assign( local_sketches[u].size()
                         + global_sketches[u].size(),
                         0
                        );

      merge( local_sketches[u].begin(),
             local_sketches[u].end(),
             global_sketches[u].begin(),
             global_sketches[u].end(),
             new_sketch.begin() );

      if (new_sketch.size() > k)
        new_sketch.resize( k );

      global_sketches[u].swap( new_sketch );
      
    }


  }

}

// void influence_oracles::alt_compute_oracles() {
//   // create the permutation of size n*ell
//   myint K = n * ell;
//   vector< mypair > perm;
//   perm.reserve( K );
//   for (myint u = 0; u < n; ++u) {
//     for (myint i = 0; i < ell; ++i) {

//       mypair tmp;
//       tmp.first = u;
//       tmp.second = i;
//       perm.push_back( tmp );
//     }
//   }

//   random_shuffle( perm.begin(), perm.end() );

//   cerr << "Permutation shuffled" << endl;

//   // group ranks by instance
//   // the pairs are in order of rank
//   vector< mypair > emptyInstance;
//   vector< vector < mypair > > instanceRanks( ell, emptyInstance );

//   for (myint r = 1; r <= K; ++r ) {
//     myint i = perm[ r - 1 ].second;
//     mypair tmp;
//     tmp.first = r;
//     tmp.second = perm[ r - 1 ].first;
//     instanceRanks[i].push_back( tmp );
//   }

//   // compute combined bottom-k rank sketches
  
//   // for each instance i, compute local bottom-k sketches
//   vector < vector< myint > > local_sketches;
//   vector < myint > empty_sketch;
//   global_sketches.assign( n, empty_sketch );
  
//   for (myint i = 0; i < ell; ++i) {
//     cerr << i << endl;
//     local_sketches.clear();
//     local_sketches.assign( n, empty_sketch );
      
//     // for each vertex in instance i by increasing rank
//     cerr << "Updating local ranks..." << endl;

//     igraph_vector_t vRankOrder;
//     igraph_vector_init( &vRankOrder, 0 );

//     for (myint j = 0; j < n; ++j) {
//       myint vertex = instanceRanks[i][j].second;
//       igraph_vector_push_back( &vRankOrder, vertex );
//     }


//     // Run reverse BFS in instance i from 'vertex', 
//     // updating sketches of discovered vertices
//     // (and not including vertex itself)
    
//     rankorder_BFS( v_instances[i],
//                    vRankOrder );
    
//       update_local_sketches( v_instances[ i ],
// 			     instanceRanks[i][j],
// 			     local_sketches );
//     }
    
//     cerr << "Performing merges..." << endl;
//     // merge local_sketches in instance i
//     // into the global sketch for each node
//     vector< myint > new_sketch;
//     for (myint u = 0; u < n; ++u) {
//       new_sketch.assign( local_sketches[u].size()
//                          + global_sketches[u].size(),
//                          0
//                         );

//       merge( local_sketches[u].begin(),
//              local_sketches[u].end(),
//              global_sketches[u].begin(),
//              global_sketches[u].end(),
//              new_sketch.begin() );

//       if (new_sketch.size() > k)
//         new_sketch.resize( k );

//       global_sketches[u].swap( new_sketch );
      
//     }


//   }

// }





