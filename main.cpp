#include <igraph.h>
#include "influence_oracles.cpp"
#include <vector>
#include <iostream>

using namespace std;

int main(void)
{
  igraph_integer_t diameter;
  igraph_t graph;
  igraph_rng_seed(igraph_rng_default(), 42);
  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 1000, 2.0/1000,
			  IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);


  vector< igraph_t* > in_inst(10, &graph);

  influence_oracles my_oracles( in_inst, (myint) 10, (myint) 10, igraph_vcount( &graph ) );

  cerr << "computing oracles...\n";

  my_oracles.compute_oracles();

  cerr << "done" << endl;

  cout << "avg reachability and estimates: ";

  for (myint i = 0; i < igraph_vcount( &graph ); ++i) {
    cout << my_oracles.average_reachability( 1 ) << ' ' << my_oracles.estimate_reachability( 1 ) << endl;
  }

  igraph_destroy(&graph);
  return 0;
}
