#include <igraph.h>
#include "influence_oracles.cpp"
#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <ctime>
#include <pthread.h>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

//global random number generator
std::random_device rd;
std::mt19937 gen(rd());

//global thread variables

igraph_t base_graph;
string graph_filename;
myint n, m;
myreal beta;
myreal alpha;
myreal int_maxprob;
myreal ext_maxprob;
string output_filename;
vector< myreal > IC_weights;
vector< myreal> node_probs;
myint ell;
influence_oracles my_oracles( 0, 0, 0 );

//function prototypes
void construct_independent_cascade( igraph_t& G, vector< myreal >& edge_weights,
                                    myreal int_maxprob);
void construct_external_influence( igraph_t& G, vector< myreal >& node_probs, myreal max_prob );
myreal construct_reachability_instance( igraph_t& G,
                                        vector< igraph_t* >& vgraphs,
                                        vector< myreal >& IC, //IC weights
                                        vector< myreal >& NP, //Node probabilities
                                        myint ell );
void sample_independent_cascade( igraph_t& G,
                                 vector< myreal >& edgne_weights,
                                 igraph_t& sample_graph );

void sample_external_influence( igraph_t& G,
                                vector< myreal >& node_probs,
                                vector< myint >& nodes_sampled );

myint forwardBFS( igraph_t* G_i, vector< myint >& vseeds, vector< myint >& v_neighborhood );
myint forwardBFS( igraph_t* G_i,
                  vector< myint >& vseeds,
                  vector< myint >& vseeds2,
                  vector< myint >& v_neighborhood,
                  int max_distance );
void my_merge( vector< myint >& sk1, vector< myint >& sk2,
               vector< myint >& res_sk,
               myint k );
void my_merge( vector< myreal >& sk1, vector< myreal >& sk2,
               vector< myreal >& res_sk,
               myint k );

void* compute_oracles_online( void * );

void print_sketch( vector< myreal >& sk1 );

void bicriteria( influence_oracles& oracles,
                 myint n,
                 myint T,
                 myreal offset,
                 //out parameters
                 vector< myint >& seeds);

void my_merge2( vector< myint >& sk1, vector< myint >& sk2,
                vector< myint >& res_sk,
                myint k );
double kempe_greedy_max( igraph_t& G,
                         vector< myreal >& IC,
                         vector< myint >& seed_set,
                         unsigned L, //number of samples
                         unsigned kk //number of seed nodes
  );

void read_params(
  bool& bdir,
  myint& n,
  unsigned& nthreads,
  string& graph_filename,
  myreal& beta,
  myreal& alpha,
  myreal& int_maxprob,
  myreal& ext_maxprob,
  string& output_filename,
  istream& is
  ) {
  is >> graph_filename;
  if (graph_filename == "ER") {
    is >> n;
  }

  is >> beta;
  is >> alpha;
  is >> int_maxprob;
  is >> ext_maxprob;
  is >> output_filename;
  is >> nthreads;
  string sdir;
  is >> sdir;
  if (sdir == "true") {
    bdir = true;
  } else {
    bdir = false;
  }
}

myreal actual_influence(
  vector< myint >& seed_set,
  igraph_t& base_graph,
  vector< myreal >& IC,
  vector< myreal >& NP,
  unsigned L ) {
  myreal activated = 0.0;

  for (unsigned i = 0; i < L; ++i) {
    if (i % 1 == 0) {
      cout << "\r                                     \r"
           << ((myreal) i )/ L * 100
        //	   << i
           << "\% done";
      cout.flush();
    }

    igraph_t G_i;
    sample_independent_cascade( base_graph, IC, G_i );
    vector< myint > ext_act;
    sample_external_influence( base_graph, NP, ext_act );

    vector< myint > v_reach;
    forwardBFS( &G_i,
                ext_act,
                seed_set,
                v_reach,
                igraph_vcount( &G_i ) );

    activated += v_reach.size();
    //    cout << ext_act.size() + seed_set.size() << ' ' << v_reach.size() << endl;

    igraph_destroy( &G_i );
  }

  cout << "\r                              \r100% done" << endl;

  return activated / L;
}



int main(int argc, char** argv) {
  unsigned int L, k;
  po::options_description opts("Kempe Simulation for CIKM'16");
  opts.add_options()
    ("help", "Show this help")
    ("graph,g", po::value<string>(&graph_filename)->required()->value_name("GRAPH.txt"),
      "The graph topology to load. This should be in the 2-3 data format.")
    ("samples,L", po::value<unsigned int>(&L)->required(),
      "The number of samples to generate in Monte-Carlo.")
    ("seeds,k", po::value<unsigned int>(&k)->required(),
      "The number of nodes in the seed set.")
    ;

  po::positional_options_description popts;
  popts.add("graph", 1);
  popts.add("samples", 1);
  popts.add("seeds", 1);

  po::variables_map args;
  po::store(po::command_line_parser(argc, argv)
            .options(opts).positional(popts).run(), args);

  if(args.count("help")) {
    cerr << opts << endl;
    return 1;
  }

  try {
    po::notify(args);
  } catch(po::required_option e) {
    cerr << "Error: " << e.get_option_name() << " required but not provided." << endl;
    cerr << opts << endl;
    return -1;
  }

  unsigned nthreads;
  bool bdir;

  ifstream file(graph_filename, fstream::in);
  file >> n >> m;

  IC_weights.reserve(m);
  igraph_vector_t edgelist;
  igraph_vector_init(&edgelist, m * 2);
  long u, v;
  for(long i = 0; i < m; i++) {
    file >> u >> v >> IC_weights[i];
    VECTOR(edgelist)[2 * i] = u;
    VECTOR(edgelist)[2 * i + 1] = v;
  }
  igraph_t base_graph;
  igraph_create(&base_graph, &edgelist, n, 1);

  n = igraph_vcount( &base_graph );

  cout << "Base graph read from " << graph_filename
       << "...";

  if (bdir)
    cout << "is directed.";
  else
    cout << "is undirected.";

  cout << endl;


  cout << "n = " << n << endl;
  cout << "m = " << igraph_ecount( &base_graph ) << endl;

  vector< myint > seed_set;

  auto start_time = std::chrono::system_clock::now();
  kempe_greedy_max(base_graph, IC_weights, seed_set, L, k);
  auto end_time = std::chrono::system_clock::now();

  cout << "Size of seed set: " << seed_set.size() << endl;
  cout << "Finished in: " << (end_time - start_time).count() << " seconds" << endl;
  //compute "actual" influence of seed set
  cout << "Estimated Influence: ";
  node_probs = vector< myreal >(n, 0.0); // no external influence
  myreal act_infl = actual_influence( seed_set, base_graph, IC_weights, node_probs, L );
  cout << act_infl << endl;

  igraph_destroy(&base_graph);
  igraph_vector_destroy(&edgelist);

  return 0;
}

void print_sketch( vector< myint >& sk1 ) {
  for (myint i = 0; i < sk1.size(); ++i) {
    cout << sk1[i] << ' ';
  }

  cout << endl;
}

void print_sketch( vector< myreal >& sk1 ) {
  for (myint i = 0; i < sk1.size(); ++i) {
    cout << sk1[i] << ' ';
  }

  cout << endl;
}

void bicriteria( influence_oracles& oracles,
                 myint n,
                 myint T,
                 myreal offset,
                 //out parameters
                 vector< myint >& seeds) {
  cerr << "offset = " << offset << endl;

  //  vector< myint > sketch;
  vector< myreal > sketch;

  myreal est_infl = offset;
  myreal max_marg = 0.0;
  myint next_node = 0;
  myreal curr_tau = 0.0;

  //  vector< myint > tmp_sketch;
  vector< myreal > tmp_sketch;
  while (est_infl < T ) {

    //select the node with the max. marginal gain
    //for each u

    max_marg = 0.0;
    //    curr_tau = oracles.estimate_reachability_sketch( sketch );
    curr_tau = oracles.estimate_reachability_uniform_sketch( sketch );
    for (myint u = 0; u < n; ++u) {
      //tmp_sketch = merge( sketch, sketch_u )
      my_merge( sketch, oracles.uniform_global_sketches[ u ], tmp_sketch, oracles.k );

      myreal tmp_marg = oracles.estimate_reachability_uniform_sketch( tmp_sketch ) - curr_tau;
      if (tmp_marg > max_marg) {
        max_marg = tmp_marg;
        next_node = u;
      }
    }

    my_merge( sketch, oracles.uniform_global_sketches[ next_node ], tmp_sketch, oracles.k );
    sketch.swap( tmp_sketch );

    seeds.push_back( next_node );

    est_infl = offset + curr_tau + max_marg;// + seeds.size();

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
void my_merge( vector< myreal >& sk1, vector< myreal >& sk2,
               vector< myreal >& res_sk,
               myint k ) {
  //	       , vector< myint >& seeds ) {

  vector< myreal >::iterator it1 = sk1.begin();
  vector< myreal >::iterator it2 = sk2.begin();

  myint l = k;
  if (l > (sk1.size() + sk2.size())) {
    l = sk1.size() + sk2.size();
  }

  res_sk.assign( l, 0 );
  vector< myreal >::iterator res_it = res_sk.begin();

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

//takes place of 'construct_reachability_instance'
//does not store any of the reachability graphs
//works with the oracle online functions
void *compute_oracles_online( void* ptr ) {
  igraph_t& G = base_graph;
  vector< myreal >& IC = IC_weights;
  vector< myreal >& NP = node_probs;
  influence_oracles& O = my_oracles;

  myint N = *( (myint*) ptr );

  myint n = igraph_vcount( &G );
  myreal offset = 0.0;
  for (myint i = 0; i < N; ++i) {
    if (i % 1 == 0) {
      cout << "\r                                     \r"
           << ((myreal) i )/ N * 100
        //	   << i
           << "\% done";
      cout.flush();
    }

    igraph_t G_i;

    sample_independent_cascade( G, IC, G_i );
    //G_i now has an IC instance. Let's select the
    //externally activated nodes
    vector< myint > ext_act;
    //need to initialize ext_act
    sample_external_influence( G, NP, ext_act );
    //remove the reachable set in G_i, of the externally activated set.
    //actually, what we want is to remove all
    //edges incident to reachable
    //set.
    vector< myint > v_reach;

    igraph_vector_t v_edges_to_remove;
    igraph_vector_init( &v_edges_to_remove, 0 );
    igraph_vector_t v_tmp;
    igraph_vector_init( &v_tmp, 0 );

    forwardBFS( &G_i, ext_act, v_reach );

    offset = v_reach.size();

    for (myint iii = 0; iii < v_reach.size(); ++iii) {
      igraph_incident( &G_i, &v_tmp, v_reach[iii], IGRAPH_ALL );

      igraph_vector_append( &v_edges_to_remove, &v_tmp );

    }
    igraph_vector_destroy(&v_tmp);

    igraph_es_t es; //vertex selector for the reachable set to remove
    igraph_es_vector( &es, &v_edges_to_remove );
    igraph_delete_edges( &G_i, es );
    igraph_es_destroy( &es );
    igraph_vector_destroy( &v_edges_to_remove );

    //All edges incident to the reachable
    //set have been removed
    //That is, H_i has been created

    //use it for the oracle computation,
    //but do then discard it

    O.compute_uniform_oracles_online_step( &G_i, i, offset );
    igraph_destroy( &G_i );
  }

  cout << "\r                               \r100% done" << endl;



}


myreal construct_reachability_instance( igraph_t& G,
                                        vector< igraph_t* >& vgraphs,
                                        vector< myreal >& IC, //IC weights
                                        vector< myreal >& NP, //Node probabilities
                                        myint ell ) {
  //create a gen. reachability instance
  vgraphs.clear();

  myint n = igraph_vcount( &G );
  myreal offset = 0.0;
  for (myint i = 0; i < ell; ++i) {
    if (i % 1 == 0) {
      cerr << "\r                                     \r"
           << ((myreal) i )/ ell * 100
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
    //All edges incident to the reachable set
    //have been removed
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



void construct_independent_cascade( igraph_t& G, vector< myreal >& edge_weights, myreal int_maxprob) {
  edge_weights.clear();

  std::uniform_real_distribution<> dis(0, int_maxprob);

  myint m = igraph_ecount( &G );

  for (myint i = 0; i < m; ++i) {
    edge_weights.push_back( dis(gen) );
  }

}

void construct_external_influence( igraph_t& G, vector< myreal >& node_probs, myreal max_prob ) {
  node_probs.clear();
  myint n = igraph_vcount( &G );
  std::uniform_real_distribution<> dis(0, max_prob );

  for (myint i = 0; i < n; ++i) {
    node_probs.push_back( dis(gen) );
  }
}

void sample_independent_cascade( igraph_t& G, vector< myreal >& edge_weights, igraph_t& sample_graph ) {

  igraph_copy( &sample_graph, &G );

  std::uniform_real_distribution<> dis(0, 1);

  myint m = igraph_ecount( &G );

  myreal cvalue;
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
                                vector< myreal >& node_probs,
                                vector< myint >& nodes_sampled ) {


  std::uniform_real_distribution<> dis(0, 1);

  myint n = igraph_vcount( &G );

  myreal cvalue;

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

myint forwardBFS( igraph_t* G_i,
                  vector< myint >& vseeds,
                  vector< myint >& v_neighborhood ) {
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

myint forwardBFS( igraph_t* G_i,
                  vector< myint >& vseeds,
                  vector< myint >& vseeds2,
                  vector< myint >& v_neighborhood,
                  int max_distance ) {
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

  for (myint i = 0; i < vseeds2.size(); ++i) {
    if (dist[ vseeds2[i] ] != 0) {
      dist[ vseeds2[i] ] = 0;
      Q.push( vseeds2[i] );
      v_neighborhood.push_back( vseeds2[i] );
    }
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
      if (dist[i] <= max_distance) {
        //i is reachable from vertex in G_i
        ++count;
      }
    }
  }

  return count;

}

double kempe_greedy_max( igraph_t& G,
                         vector< myreal >& IC,
                         vector< myint >& seed_set,
                         unsigned L, //number of samples
                         unsigned kk //number of seed nodes
  ) {
  double eact = 0;
  seed_set.clear();
  seed_set.reserve( kk );
  double max_marge;
  vector< myint > v_reach;

  igraph_t* G_i;

  for (unsigned iter = 0; iter < kk; ++iter) {

    cerr << "Selecting seed node " << iter << endl;
    seed_set.push_back( 0 ); // we're going to be adding a new
    //seed node, identity to be determined
    myint next_node = 0;
    max_marge = 0.0;
    //try all nodes
    for (unsigned i = 0; i < n; ++i) {

      double tmp_eact = 0.0;
      for (unsigned j = 0; j < L; ++j) {
        v_reach.clear();
        G_i = new igraph_t;
        sample_independent_cascade(G, IC, *G_i );
        seed_set[ iter ] = i; // test the ith seed node
        forwardBFS( G_i, seed_set, v_reach );
        tmp_eact += v_reach.size();
        igraph_destroy( G_i );
      }

      tmp_eact /= L;
      double marge_i = tmp_eact - eact;
      if (marge_i > max_marge) {
        max_marge = marge_i;
        next_node = i;
      }
    }

    //we've found node with max marginal gain
    seed_set[iter] = next_node;
    eact += max_marge;
  }

  return eact;
}
