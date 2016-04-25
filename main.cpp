#include <igraph.h>
#include "influence_oracles.cpp"
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdio>

using namespace std;

//global random number generator
std::random_device rd;
std::mt19937 gen(rd());


int main(int argc, char** argv) {
  if (argc < 2) {

    cerr << "Usage: " << argv[0] << " <base graph filename>";
    return 1;
  }
  
  //  string bg_filename( argv[1] );

  //  ifstream ifile( bg_filename.c_str());

  FILE* fp;
  fp = fopen( argv[1], "r" );

  igraph_t base_graph;
  igraph_read_graph_edgelist( &base_graph, fp, 0, true ); 

  vector< double > IC_weights;

  construct_independent_cascade( base_graph, IC_weights );

  vector< double> node_probs;

  //this is a simple model of external influence
  construct_external_influence( base_graph, node_probs );

  
  
  fclose( fp );
  return 0;
}

void construct_reachability_instance( igraph_t& G, 
				      vector< igraph_t >& vgraphs,
				      myint ell ) {
  //create a gen. reachability instance
  igraph_t G_i;
  myint n = igraph_vcount( &G );
  for (myint i = 0; i < ell; ++i) {
    //    igraph_copy( &G_i, &G );
    sample_independent_cascade( G, IC, G_i );
    //G_i now has an IC instance. Let's select the
    //externally activated nodes
    igraph_vector_t ext_act;
    //need to initialize ext_act
    sample_external_influence( G, NP, ext_act );
    

    //remove the reachable set in G_i, of the externally activated set.
    igraph_vs_t vids; //vertex selector for the seeds
    igraph_vs_vector( &vids, &ext_act );
    //this function only does neighborhood of the individual vertices,
    //not the set. //better to just implement a BFS
    igraph_neighborhood( &G, &res, vids, n, IGRAPH_OUT );

  }

}



int test_estimation()
{
  igraph_integer_t diameter;
  igraph_t graph;
  igraph_rng_seed(igraph_rng_default(), 42);
  myint ell = 10;
  myint n = 100000;
  vector < igraph_t > vgraphs( ell, graph );

  for (myint i = 0; i < ell; ++i) {
    igraph_erdos_renyi_game(&(vgraphs[i]), IGRAPH_ERDOS_RENYI_GNP, n, 1.0/ n,
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



void construct_independent_cascade( igraph_t& G, vector< double >& edge_weights ) {
  edge_weights.clear();

  std::random_device rd;


  myint m = igraph_ecount( &G );

  for (myint i = 0; i < m; ++i) {
    edge_weights.push_back( dis(gen) );
  }

}

void construct_external_influence( igraph_t& G, vector< double >& node_probs ) {
  node_probs.clear();

  std::uniform_real_distribution<> dis(0, 1);

  myint n = igraph_vcount( &G );

  for (myint i = 0; i < n; ++i) {
    node_probs.push_back( dis(gen) );
  }
}

void sample_independent_cascade( igraph_t& G, vector< double >& edge_weights, igraph_t& sample_graph ) {
  
  igraph_copy( &sample_graph, &G );

  std::uniform_real_distribution<> dis(0, 1);

  myint m = igraph_ecount( &G );

  double cvalue;
  igraph_vector_t edges_to_delete;
  igraph_vector_init( &edges_to_delete, 0L );
  igraph_vector_reserve( &edges_to_delete, m );

  for (myint i = 0; i < m; ++i ) {
    //toss a coin
    cvalue = dis(gen);
    if (cvalue < edge_weights[i]) {
      //keep edge i
    } else {
      //delete edge i
      igraph_vector_push_back( &edges_to_delete, i );
    }
  }

  igraph_es_t es; //edge selector of edges to destroy
  igraph_es_vector( &es, &edges_to_delete );

  //actually delete the edges from the sample graph
  igraph_delete_edges( &sample_graph, es );

  igraph_es_destroy( &es );
  igraph_vector_destroy( &edges_to_delete );
  
}

void sample_external_influence( igraph_t& G, vector< double >& node_probs, igraph_vector_t& nodes_sampled ) {
  

  std::uniform_real_distribution<> dis(0, 1);

  myint n = igraph_vcount( &G );

  double cvalue;

  igraph_vector_reserve( &nodes_sampled, n );

  for (myint i = 0; i < n; ++i ) {
    //toss a coin
    cvalue = dis(gen);
    if (cvalue < node_probs[i]) {
      //sample this node i
      igraph_vector_push_back( &nodes_sampled, i );
    } 
  }

}

