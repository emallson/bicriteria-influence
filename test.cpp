#include <igraph.h>
#include "influence_oracles.cpp"
#include <iostream>
#include <ctime>

int main() {
  myint n = 100;

  igraph_t base_graph;

  // igraph_erdos_renyi_game(&base_graph, 
  //                         IGRAPH_ERDOS_RENYI_GNP, 
  //                         n, 2.0 / n,
  //                         IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  igraph_empty( &base_graph, 100, false );
  igraph_add_edge( &base_graph, 0, 1 );
  igraph_add_edge( &base_graph, 0, 2 );
  igraph_add_edge( &base_graph, 0, 3 );

  vector< igraph_t* > v1;
  v1.push_back( &base_graph );
  influence_oracles my_oracle( v1, 1, 1, n );

  cout << "Testing igraph_neighborhood..." << endl;

  igraph_vs_t vs;
  //  igraph_vs_all( &vs );
  igraph_vector_t myvec;
  igraph_vector_init( &myvec, 0 );
  igraph_vector_push_back( &myvec, 10 );
  igraph_vector_push_back( &myvec, 0 );

  igraph_vs_vector( &vs, &myvec );

  clock_t t_start = clock();
  igraph_vector_ptr_t res;
  igraph_vector_ptr_init( &res, 0 );
  igraph_neighborhood( &base_graph,
                       &res,
                       vs,
                       4,
                       IGRAPH_IN );

  clock_t t_end = clock();
  double t_elapsed = (double) (t_end - t_start) / CLOCKS_PER_SEC;

  cout << "Done in " << t_elapsed << " seconds" << endl;

  cout << "Testing my BFS..." << endl;
  t_start = clock();
  vector< vector< myint > > local_sketches;
  vector< myint > empty_sketch;
  
  local_sketches.assign( n, empty_sketch );
  for (myint i = 0; i < n; ++i ) {
    mypair tmp;
    tmp.first = n - i;
    tmp.second = i;
    my_oracle.update_local_sketches( &base_graph,
                           tmp,
                           local_sketches );
                           
  }
  
  t_end = clock();
  
  t_elapsed = (double) (t_end - t_start) / CLOCKS_PER_SEC;

  cout << "Done in " << t_elapsed << " seconds" << endl;

  for (myint i = 0; i < igraph_vector_ptr_size( &res );
       ++i ) {
    igraph_vector_t* v;
    v = (igraph_vector_t* ) igraph_vector_ptr_e( &res, i );
    cout << "node " << i << ": ";
    for (myint j = 0; j < igraph_vector_size( v ); ++j) {
      cout << igraph_vector_e( v, j ) << ' ';
    }

    cout << endl;
  }

}
