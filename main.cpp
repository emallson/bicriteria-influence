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
  myint ell = 10;
  myint n = 1000;
  vector < igraph_t > vgraphs( ell, graph );

  for (myint i = 0; i < ell; ++i) {
    igraph_erdos_renyi_game(&(vgraphs[i]), IGRAPH_ERDOS_RENYI_GNP, n, 1.0/1000,
			  IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  }

  influence_oracles my_oracles( vgraphs, ell, (myint) 10, n );

  cerr << "computing oracles...\n";

  my_oracles.compute_oracles();

  cerr << "done" << endl;

  cout << "avg reachability and estimates: ";
  for (myint i = 0; i < n; ++i) {
    cout << i << ' ' << my_oracles.average_reachability( i ) << ' ' << my_oracles.estimate_reachability( i ) << endl;
  }

  return 0;
}
