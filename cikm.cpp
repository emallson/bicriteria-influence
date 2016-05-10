#include <igraph.h>
#include "influence_oracles.cpp"
#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
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
myreal construct_reachability_instance( igraph_t& G,
                                        vector< igraph_t* >& vgraphs,
                                        vector< myreal >& IC, //IC weights
                                        vector< myreal >& NP, //Node probabilities
                                        myint ell );
vector<bool> sample_independent_cascade( igraph_t& G,
                                 vector< myreal >& edgne_weights);

myint forwardBFS( igraph_t* G_i, vector< myint >& vseeds, vector< myint >& v_neighborhood, vector<bool>& deleted);

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
  unsigned L ) {
  myreal activated = 0.0;

  for (unsigned i = 0; i < L; ++i) {
    // if (i % 1 == 0) {
    //   cout << "\r                                     \r"
    //        << ((myreal) i )/ L * 100
    //     //	   << i
    //        << "\% done";
    //   cout.flush();
    // }

    auto deleted = sample_independent_cascade( base_graph, IC);

    vector< myint > v_reach;
    forwardBFS( &base_graph, seed_set, v_reach , deleted);

    activated += v_reach.size();
    //    cout << ext_act.size() + seed_set.size() << ' ' << v_reach.size() << endl;

  }

  // cout << "\r                              \r100% done" << endl;

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

  std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
  start_time = std::chrono::system_clock::now();
  kempe_greedy_max(base_graph, IC_weights, seed_set, L, k);
  end_time = std::chrono::system_clock::now();

  for(int i = 0; i < k; i++) {
    cout << seed_set[i] << " ";
  }
  cout << endl;
  cout << "Size of seed set: " << seed_set.size() << endl;
  cout << "Finished in: " << (end_time - start_time).count() << " seconds" << endl;
  //compute "actual" influence of seed set
  cout << "Estimated Influence: ";
  myreal act_infl = actual_influence( seed_set, base_graph, IC_weights, L );
  cout << act_infl << endl;

  igraph_destroy(&base_graph);
  igraph_vector_destroy(&edgelist);

  return 0;
}

vector<bool> sample_independent_cascade( igraph_t& G, vector< myreal >& edge_weights) {

  std::uniform_real_distribution<> dis(0, 1);

  myint m = igraph_ecount( &G );
  vector<bool> deleted(m, false);

  myreal cvalue;
  for (myint i = 0; i < m; ++i ) {
    //toss a coin
    cvalue = dis(gen);
    if (cvalue < edge_weights[i]) {
      //keep edge i
    } else {
      //delete edge i
      deleted[i] = true;
    }
  }
  return deleted;
}



myint forwardBFS( igraph_t* G,
                  vector< myint >& vseeds,
                  vector< myint >& v_neighborhood,
                  vector< bool >& deleted ) {
  queue <myint> Q;
  vector < int > dist;
  myint n = igraph_vcount( G );
  dist.reserve( n );
  for (myint i = 0; i < n; ++i) {
    dist.push_back( -1 ); //infinity
  }

  for (myint i = 0; i < vseeds.size(); ++i) {
    dist[ vseeds[i] ] = 0;
    Q.push( vseeds[i] );
    v_neighborhood.push_back( vseeds[i] );
  }

  while (!(Q.empty())) {

    myint current = Q.front();
    Q.pop();
    //get forwards neighbors of current
    igraph_vector_t neis;
    igraph_vector_init( &neis, 0 );
    igraph_neighbors( G, &neis, current, IGRAPH_OUT );
    for (myint i = 0; i < igraph_vector_size( &neis ); ++i) {
      myint aneigh = VECTOR( neis )[i];
      igraph_integer_t eid;
      igraph_get_eid(G, &eid, current, aneigh, true, false);
      if (deleted[ eid ]) {
        continue;
      }
      if (dist[ aneigh ] == -1 ) { //if aneigh hasn't been discovered yet
        dist[ aneigh ] = dist[ current ] + 1;
        Q.push( aneigh );
        //	igraph_vector_push_back( &v_neighborhood, aneigh );
        v_neighborhood.push_back( aneigh );

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

      if(std::find(seed_set.begin(), seed_set.end(), i) != seed_set.end()) {
        // don't duplicate node
        continue;
      }

      double tmp_eact = 0.0;
      for (unsigned j = 0; j < L; ++j) {
        v_reach.clear();
        auto deleted = sample_independent_cascade(G, IC);
        seed_set[ iter ] = i; // test the ith seed node
        forwardBFS( &G, seed_set, v_reach , deleted);
        tmp_eact += v_reach.size();
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
