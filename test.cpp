#include <igraph.h>
#include "influence_oracles.cpp"
#include <iostream>
#include <ctime>

int main() {
  myint n = 10000;

  igraph_t base_graph;

  igraph_erdos_renyi_game(&base_graph, 
                          IGRAPH_ERDOS_RENYI_GNP, 
                          n, 2.0 / n,
                          IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  vector< igraph_t* > v1;
  v1.push_back( &base_graph );
  influence_oracles my_oracle( v1, 1, 1, n );

  cout << "Testing igraph_neighborhood..." << endl;
  clock_t t_start = clock();

  igraph_vs_t vs;
  igraph_vs_all( &vs );
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

}
