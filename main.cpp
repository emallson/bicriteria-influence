#include <igraph.h>
#include "influence_oracles.cpp"
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

//global random number generator
std::random_device rd;
std::mt19937 gen(rd());

//function prototypes
void construct_independent_cascade( igraph_t& G, vector< double >& edge_weights );
void construct_external_influence( igraph_t& G, vector< double >& node_probs, double max_prob );
double construct_reachability_instance( igraph_t& G, 
				      vector< igraph_t* >& vgraphs,
				      vector< double >& IC, //IC weights
				      vector< double >& NP, //Node probabilities
				      myint ell );
void sample_independent_cascade( igraph_t& G, 
				 vector< double >& edgne_weights, 
				 igraph_t& sample_graph );

void sample_external_influence( igraph_t& G, 
				vector< double >& node_probs, 
				vector< myint >& nodes_sampled );

myint forwardBFS( igraph_t* G_i, vector< myint >& vseeds, vector< myint >& v_neighborhood );
void my_merge( vector< myint >& sk1, vector< myint >& sk2,
	       vector< myint >& res_sk, 
               myint k );

void bicriteria( influence_oracles& oracles, 
                 myint n, 
                 double beta,
                 double offset,
		 //out parameters
		 vector< myint >& seeds);

void my_merge2( vector< myint >& sk1, vector< myint >& sk2,
	       vector< myint >& res_sk, 
		myint k );


int main(int argc, char** argv) {
  igraph_t base_graph;
  myint n;
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <input filename>\n";
    return 1;
  }

  string str_ifile( argv[1] );
  ifstream ifile( str_ifile.c_str() );

  string bg_graph_type;
  getline( ifile, bg_graph_type );
  //  ifile >> bg_graph_type;

  if (bg_graph_type == "ER") {
    ifile >> n;
    igraph_erdos_renyi_game(&base_graph, 
                          IGRAPH_ERDOS_RENYI_GNP, 
                          n, 2.0 / n,
                          IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);


  } else {
    //    ifstream ifile( bg_filename.c_str());
    FILE* fp;
    fp = fopen( bg_graph_type.c_str(), "r" );
    //if the graph is directed
    igraph_read_graph_edgelist( &base_graph, fp, 0, true ); 

    fclose( fp );
    n = igraph_vcount( &base_graph );

    cout << "Base graph read from " << bg_graph_type << "..." << endl;
  }
  
  cout << "n = " << n << endl;
  cout << "m = " << igraph_ecount( &base_graph ) << endl;
  double beta;
  ifile >> beta;
  myint C = 2;
  double alpha;
  ifile >> alpha;
  myint K = 1.0 / (C * alpha);
  double delta = 0.5;
  myint ell = log( 2 / delta ) / (alpha * alpha) / 2;
  myint k_cohen = 3 * (myint)( log ( ((double) n) ) );

  double ext_maxprob;
  ifile >> ext_maxprob;
  cout << "beta = " << beta << endl;
  cout << "alpha = " << alpha << endl;
  cout << "ell = " << ell << endl;
  cout << "K = " << K << endl;
  cout << "ext_maxprob = " << ext_maxprob << endl;

  system("sleep 2");

  cout << "Constructing the IC model..." << endl;
  vector< double > IC_weights;
  construct_independent_cascade( base_graph, IC_weights );

  vector< double> node_probs;
  //this is a simple model of external influence
  cout << "Constructing external influence..." << endl;
  construct_external_influence( base_graph, node_probs, ext_maxprob );

  //create the set of graphs for the alg.
  cerr << "Constructing the gen. reach. instance...\n";

  vector< igraph_t* > v_graphs; // the ell graphs
  double offset = construct_reachability_instance( base_graph, 
				   v_graphs, 
				   IC_weights,
				   node_probs,
				   ell );
  cerr << "offset=" << offset << endl;

  system("sleep 1");


  //create the oracles
  influence_oracles my_oracles( v_graphs, ell, k_cohen, n );
  
  cerr << "computing oracles...\n";

  my_oracles.compute_oracles();

  cerr << "done" << endl;

  // cout << "avg reachability and estimates: ";
  // for (myint i = 0; i < n; ++i) {
  //   cout << i << ' ' << my_oracles.average_reachability( i ) << ' ' << my_oracles.estimate_reachability( i ) << endl;
  // }

  //run the bicriteria alg.
  vector< myint > seed_set;
  bicriteria( my_oracles, n, beta, offset,
	      seed_set );

  cout << "Size of seed set: " << seed_set.size() << endl;

  return 0;
}

void print_sketch( vector< myint >& sk1 ) {
  for (myint i = 0; i < sk1.size(); ++i) {
    cout << sk1[i] << ' ';
  }

  cout << endl;
}

void bicriteria( influence_oracles& oracles, 
                 myint n, 
                 double beta,
                 double offset,
		 //out parameters
		 vector< myint >& seeds) {
  cerr << "offset = " << offset << endl;

  vector< myint > sketch;
  vector< vector < myint > > track_sketches;


  double est_infl = offset;
    double max_marg = 0.0;
  myint next_node = 0;
  double curr_tau = 0.0;

  vector< myint > tmp_sketch;
  while (est_infl < beta * n ) {
    track_sketches.push_back( sketch );
    //select the node with the max. marginal gain
    //for each u

    max_marg = 0.0;
    curr_tau = oracles.estimate_reachability_sketch( sketch );
    for (myint u = 0; u < n; ++u) {
      //tmp_sketch = merge( sketch, sketch_u )
      my_merge( sketch, oracles.global_sketches[ u ], tmp_sketch, oracles.k );
      
      double tmp_marg = oracles.estimate_reachability_sketch( tmp_sketch ) - curr_tau;
      if (tmp_marg > max_marg) {
	max_marg = tmp_marg;
	next_node = u;
      }
    }

    //pick the max. one, update sketch
    // cout << "Merging sketches:\n";
    // print_sketch( sketch );
    // print_sketch( oracles.global_sketches[ next_node ] );
    my_merge( sketch, oracles.global_sketches[ next_node ], tmp_sketch, oracles.k );
    sketch.swap( tmp_sketch );
    // cout << "to get sketch:\n";
    // print_sketch( sketch );

    //    print_sketch( oracles.global_sketches[ next_node ] );
    // for (myint jj = 0; jj < track_sketches.size(); ++jj) {
    //   my_merge( track_sketches[jj], oracles.global_sketches[ next_node ], tmp_sketch, oracles.k );
      // cout << "tau_init: " << oracles.estimate_reachability_sketch( track_sketches[jj] )
      // 	   << " tau_addu: " << oracles.estimate_reachability_sketch( tmp_sketch )
      // 	   << " margin: " << oracles.estimate_reachability_sketch( tmp_sketch ) - 
      // 	oracles.estimate_reachability_sketch( track_sketches[jj] )
      // 	   << endl;

      //      print_sketch( track_sketches[jj] );
      //      print_sketch( tmp_sketch );
    //   }

    seeds.push_back( next_node );

    est_infl = offset + curr_tau + max_marg + seeds.size();

    cerr << est_infl << ' ' << curr_tau << ' ' <<  max_marg << ' ' << next_node << endl;
    cerr << "perm rank " << sketch.back() << endl;
  }


}

void my_merge2( vector< myint >& sk1, vector< myint >& sk2,
	       vector< myint >& res_sk, 
               myint k ) {
  //this is not really a merge, seems to be what cohen is saying
  //in the 1997 paper. Rather is coordinate-wise min.
  
  vector< myint >::iterator it1 = sk1.begin();
  vector< myint >::iterator it2 = sk2.begin();

  myint l = sk1.size();
  if (sk2.size() > l)
    l = sk2.size();

  res_sk.assign( l, 0 );

  vector< myint >::iterator res_it = res_sk.begin();

  myint s = 0; //size

  while (s < l) {
    //add the smallest of coordinates
    if (it1 != sk1.end()) {
      if (it2 != sk2.end()) {
	if ( (*it1) < (*it2) ) {
	  (*res_it) = (*it1);
	  ++it1;
	  ++it2;
	} else {
	  if ( (*it1) > (*it2) ) {
	    (*res_it) = (*it2);
	    ++it2;
	    ++it1;
	  } else {
	    //the values are equal
	    //but we can only have
	    //an element appear once
	    //so insert it, and increment both
	    (*res_it) = (*it1);
	    ++it1;
	    ++it2;
	  }
	}
      } else { //we're done with sketch2.
	//just add from sk1
	(*res_it) = (*it1);
	++it1;
      }
    } else {
      //done with sk1,
      //add from sk2
      (*res_it) = (*it2);
      ++it2;
    }

    ++res_it;
    ++s;
  }

}

void my_merge( vector< myint >& sk1, vector< myint >& sk2,
	       vector< myint >& res_sk, 
               myint k ) {
  //	       , vector< myint >& seeds ) {

  vector< myint >::iterator it1 = sk1.begin();
  vector< myint >::iterator it2 = sk2.begin();

  myint l = k;
  if (l > (sk1.size() + sk2.size())) {
    l = sk1.size() + sk2.size();
  }

  res_sk.assign( l, 0 );
  vector< myint >::iterator res_it = res_sk.begin();

  myint s = 0; //size

  while (s < l) {
    //add the smallest next one
    if (it1 != sk1.end()) {
      if (it2 != sk2.end()) {
	if ( (*it1) < (*it2) ) {
	  (*res_it) = (*it1);
	  ++it1;
	} else {
	  if ( (*it1) > (*it2) ) {
	    (*res_it) = (*it2);
	    ++it2;
	  } else {
	    //the values are equal
	    //but we can only have
	    //an element appear once
	    //so insert it, and increment both
	    (*res_it) = (*it1);
	    ++it1;
	    ++it2;
	  }
	}
      } else { //we're done with sketch2.
	//just add from sk1
	(*res_it) = (*it1);
	++it1;
      }
    } else {
      //done with sk1,
      //add from sk2
      (*res_it) = (*it2);
      ++it2;
    }

    ++res_it;
    ++s;
  }


}


double construct_reachability_instance( igraph_t& G, 
				      vector< igraph_t* >& vgraphs,
				      vector< double >& IC, //IC weights
				      vector< double >& NP, //Node probabilities
				      myint ell ) {
  //create a gen. reachability instance
  vgraphs.clear();

  myint n = igraph_vcount( &G );
  double offset = 0.0;
  for (myint i = 0; i < ell; ++i) {
    if (i % 1 == 0) {
      cerr << "\r                                     \r"
	   << ((double) i )/ ell * 100 
	//	   << i
	   << "\% done";
    }

    igraph_t* G_i = new igraph_t; //want a new address in memory each iteration    

    sample_independent_cascade( G, IC, *G_i );
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

    //    cerr << "\n" << ' ';
    forwardBFS( G_i, ext_act, v_reach );

    offset += v_reach.size();

    for (myint i = 0; i < v_reach.size(); ++i) {

      igraph_incident( G_i, &v_tmp, v_reach[i], IGRAPH_ALL );

      igraph_vector_append( &v_edges_to_remove, &v_tmp );

    }

    //cerr << igraph_vector_size( &v_edges_to_remove ) << "\n";

    igraph_vector_destroy(&v_tmp);

    igraph_es_t es; //vertex selector for the reachable set to remove
    igraph_es_vector( &es, &v_edges_to_remove );
    igraph_delete_edges( G_i, es );
    igraph_es_destroy( &es );
    igraph_vector_destroy( &v_edges_to_remove );

    //    cerr << "BFS done\n";
    //All edges incident to the reachable set have been removed
    //That is, H_i has been created
    vgraphs.push_back( G_i );
  }

  cout << "\r                               \r100% done" << endl;

  return offset / ell;
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

void construct_external_influence( igraph_t& G, vector< double >& node_probs, double max_prob ) {
  node_probs.clear();
  myint n = igraph_vcount( &G );
  std::uniform_real_distribution<> dis(0, max_prob );

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
