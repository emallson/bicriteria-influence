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

//function prototypes
void construct_independent_cascade( igraph_t& G, vector< double >& edge_weights );
void construct_external_influence( igraph_t& G, vector< double >& node_probs );
void construct_reachability_instance( igraph_t& G, 
				      vector< igraph_t* >& vgraphs,
				      vector< double >& IC, //IC weights
				      vector< double >& NP, //Node probabilities
				      myint ell );
void sample_independent_cascade( igraph_t& G, 
				 vector< double >& edge_weights, 
				 igraph_t& sample_graph );

void sample_external_influence( igraph_t& G, 
				vector< double >& node_probs, 
				vector< myint >& nodes_sampled );

myint forwardBFS( igraph_t* G_i, vector< myint >& vseeds, vector< myint >& v_neighborhood );



int main(int argc, char** argv) {
  if (argc < 2) {

    cerr << "Usage: " << argv[0] << " <base graph filename>\n";
    return 1;
  }
  
  //  string bg_filename( argv[1] );

  //  ifstream ifile( bg_filename.c_str());

  //  FILE* fp;
  //  fp = fopen( argv[1], "r" );

  igraph_t base_graph;
  //  igraph_read_graph_edgelist( &base_graph, fp, 0, true ); 
  myint n = 50000;
  igraph_erdos_renyi_game(&base_graph, IGRAPH_ERDOS_RENYI_GNP, n, 1.0/ n,
			  IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  //  fclose( fp );

  cerr << "Constructing the IC model...\n";
  vector< double > IC_weights;
  construct_independent_cascade( base_graph, IC_weights );

  vector< double> node_probs;

  //this is a simple model of external influence
  cerr << "Constructing external influence...\n";
  construct_external_influence( base_graph, node_probs );

  //create the set of graphs for the alg.
  cerr << "Constructing the gen. reach. instance...\n";

  myint ell = 10;
  myint k_cohen = ell;

  vector< igraph_t* > v_graphs; // the ell graphs
  construct_reachability_instance( base_graph, 
				   v_graphs, 
				   IC_weights,
				   node_probs,
				   ell );

  //create the oracles
  influence_oracles my_oracles( v_graphs, ell, k_cohen, n );
  
  cerr << "computing oracles...\n";

  my_oracles.compute_oracles();

  cerr << "done" << endl;

  cout << "avg reachability and estimates: ";
  for (myint i = 0; i < n; ++i) {
    cout << i << ' ' << my_oracles.average_reachability( i ) << ' ' << my_oracles.estimate_reachability( i ) << endl;
  }

  return 0;
}

void construct_reachability_instance( igraph_t& G, 
				      vector< igraph_t* >& vgraphs,
				      vector< double >& IC, //IC weights
				      vector< double >& NP, //Node probabilities
				      myint ell ) {
  //create a gen. reachability instance
  vgraphs.clear();

  myint n = igraph_vcount( &G );
  for (myint i = 0; i < ell; ++i) {
    if (i % 1 == 0) {
      cerr << "\r                                     \r"
	   << ((double) i )/ ell * 100 
	//	   << i
	   << "\% done";
    }

    igraph_t G_i; //want a new address in memory each iteration

    sample_independent_cascade( G, IC, G_i );
    //G_i now has an IC instance. Let's select the
    //externally activated nodes
    vector< myint > ext_act;
    //need to initialize ext_act
    sample_external_influence( G, NP, ext_act );
    //remove the reachable set in G_i, of the externally activated set.
    //actually, what we want is to remove all edges incident to reachable
    //set.
    vector< myint > v_reach;

    igraph_vector_t v_edges_to_remove;
    igraph_vector_init( &v_edges_to_remove, 0 );
    igraph_vector_t v_tmp;
    igraph_vector_init( &v_tmp, 0 );

    cerr << "\n" << forwardBFS( &G_i, ext_act, v_reach ) << ' ';

    for (myint i = 0; i < v_reach.size(); ++i) {

      igraph_incident( &G_i, &v_tmp, v_reach[i], IGRAPH_ALL );

      igraph_vector_append( &v_edges_to_remove, &v_tmp );

    }

    cerr << igraph_vector_size( &v_edges_to_remove ) << "\n";

    igraph_vector_destroy(&v_tmp);

    igraph_es_t es; //vertex selector for the reachable set to remove
    igraph_es_vector( &es, &v_edges_to_remove );
    igraph_delete_edges( &G_i, es );
    igraph_es_destroy( &es );
    igraph_vector_destroy( &v_edges_to_remove );

    cerr << "BFS done\n";
    //All edges incident to the reachable set have been removed
    //That is, H_i has been created
    vgraphs.push_back( &G_i );
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

  std::uniform_real_distribution<> dis(0, 1);

  myint m = igraph_ecount( &G );

  for (myint i = 0; i < m; ++i) {
    edge_weights.push_back( dis(gen) );
  }

}

void construct_external_influence( igraph_t& G, vector< double >& node_probs ) {
  node_probs.clear();

  std::uniform_real_distribution<> dis(0, 0.1);

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

void sample_external_influence( igraph_t& G, 
				vector< double >& node_probs, 
				vector< myint >& nodes_sampled ) {
  

  std::uniform_real_distribution<> dis(0, 1);

  myint n = igraph_vcount( &G );

  double cvalue;

  nodes_sampled.clear();

  for (myint i = 0; i < n; ++i ) {
    //toss a coin
    cvalue = dis(gen);
    if (cvalue < node_probs[i]) {
      //sample this node i
      //igraph_vector_push_back( &nodes_sampled, i );
      nodes_sampled.push_back( i );
    } 
  }

}

myint forwardBFS( igraph_t* G_i, vector< myint >& vseeds, vector< myint >& v_neighborhood ) {
  queue <myint> Q;
  vector < int > dist;
  myint n = igraph_vcount( G_i );
  dist.reserve( n );
  for (myint i = 0; i < n; ++i) {
    dist.push_back( -1 ); //infinity
  }

  //  igraph_vector_t v_edges_to_remove;
  //  igraph_vector_init( &v_edges_to_remove, 0);

  for (myint i = 0; i < vseeds.size(); ++i) {
    dist[ vseeds[i] ] = 0;
    Q.push( vseeds[i] );
    //    igraph_vector_push_back( &v_neighborhood, vseeds[i] );
    v_neighborhood.push_back( vseeds[i] );
  }

  while (!(Q.empty())) {

    myint current = Q.front();
    Q.pop();
    //get forwards neighbors of current
    igraph_vector_t neis;
    igraph_vector_init( &neis, 0 );
    igraph_neighbors( G_i, &neis, current, IGRAPH_OUT );
    for (myint i = 0; i < igraph_vector_size( &neis ); ++i) {
      myint aneigh = VECTOR( neis )[i];
      if (dist[ aneigh ] == -1 ) { //if aneigh hasn't been discovered yet
	dist[ aneigh ] = dist[ current ] + 1;
	Q.push( aneigh );
	//	igraph_vector_push_back( &v_neighborhood, aneigh );
	v_neighborhood.push_back( aneigh );

	//flag this edge for removal from the graph
	//	myint eid;
	//	igraph_get_eid( G_i, &eid, current, 
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
