#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <stdint.h>
#include <iostream>
#include <random>
#include <thread>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <utility>
#include "progressBar.h"
#include "networkit/graph/Graph.hpp"
//#include "networkit/distance/Dijkstra.hpp"
#include <string>
#include "mytimer.h"
#include "networkit/centrality/DegreeCentrality.hpp"
#include "networkit/centrality/EstimateBetweenness.hpp"

using vertex = int;
using dist = uint32_t;
const vertex null_vertex = round(std::numeric_limits<vertex>::max()/2);
const dist null_distance = round(std::numeric_limits<dist>::max()/2);
//
// NOTE: Currently only unweighted and undirected graphs are supported.
//
class HighwayLabelling {
 public:
  HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type);
  ~HighwayLabelling();

    void ConstructUnweightedHighwayLabelling();
    void ConstructWeightedHighwayLabelling();
    void ConstructWeightedHighwayLabellingWithFlag();
    //void RecursivePathFlag(NetworKit::Dijkstra* sssp, vertex source, vertex current_vertex, uint8_t* pruning_flag);
    void ConstructUnweightedNaiveLabelling();
    void ConstructWeightedNaiveLabelling();

  int GetNumberOfNodes(){return V;};
  long LabellingSize();
  long NaiveLabellingSize();
  dist min(dist a, dist b);

  // Returns distance vetween vertices v and w if they are connected.
  dist QueryDistance(vertex s, vertex t);
    dist QueryDistanceQuadratic(vertex s, vertex t);
  dist NaiveQueryDistance(vertex s, vertex t);
  dist BFSQuery(vertex s, vertex t);
  dist DijkstraQuery(vertex s, vertex t);
  dist BoundedSearch(vertex s, vertex t);
  dist P2PBFS(vertex s, vertex landmark, dist ub, dist partial, dist* to_vertices);
  void GetBeerStores(std::vector<vertex>* beer_stores);

 private:
  vertex V;  // total number of vertices
  vertex E; // total number of edges
  vertex L; // total number of landmarks
  NetworKit::Graph &graph;
  vertex **vertices = nullptr;
    dist **distances;
    dist **query_distances;
  dist **highway;
  dist *C = nullptr;
  vertex *ordering;
  vertex *reverse_ordering;
  dist** naive_labeling = nullptr;

};


HighwayLabelling::~HighwayLabelling() {

  for(int i = 0; i < V; i++) {
    delete [] query_distances[i];
  }
  delete [] query_distances;

  if(vertices != nullptr) {
        for(int i = 0; i < V; i++) {
            delete [] vertices[i];
        }
        delete [] vertices;
    }
  delete [] C;
  for(int i = 0; i < L; i++)
    delete [] highway[i];
  delete [] highway;
  if(naive_labeling != nullptr) {
    for (int i = 0; i < V; i++)
        delete[] naive_labeling[i];
    delete[] naive_labeling;
  }

}

HighwayLabelling::HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type)
        : graph(g) {
  V = 0; E = 0; L = l;
  graph = g;
  V = this->graph.numberOfNodes();
  E = this->graph.numberOfEdges();
    auto *ordering_rank = new std::pair<double,vertex>[graph.numberOfNodes()];
    this->ordering = new vertex[graph.numberOfNodes()];
    this->reverse_ordering = new vertex[graph.numberOfNodes()];

    this->graph.parallelForNodes([&] (vertex i){
        assert(graph.hasNode(i));
        this->ordering[i]=null_vertex;
        this->reverse_ordering[i]=null_vertex;
        ordering_rank[i] = {0,i};
    });


    double centr_time = 0.0;

    if(ordering_type==0){

        INFO("BY DEGREE");
        NetworKit::DegreeCentrality* rank = new NetworKit::DegreeCentrality(graph);
        rank->run();
        this->graph.forNodes([&] (vertex i){
            assert(graph.hasNode(i));
            ordering_rank[i]=std::make_pair(rank->score(i),i);
        });
        delete rank;
    }

    if(ordering_type==1){

        INFO("RANDOM ORDERING");
        std::vector<vertex> t(graph.numberOfNodes());
        for(vertex i = 0; i < graph.numberOfNodes(); i++) t[i] = i;
        std::shuffle(t.begin(), t.end(), std::mt19937());
        for(vertex i = 0; i < graph.numberOfNodes(); i++) ordering_rank[i]=std::make_pair(t[i],i);
    }

    if(ordering_type==2){

        INFO("BY APX BETW");
        mytimer local_constr_timer;
        double max_time = 30.0;
        double cumulative_time = 0.0;
        double fract = 0.33;
        double n_samples =  round(std::pow((double)graph.numberOfNodes(),fract));

        while(cumulative_time<max_time && n_samples<(double)graph.numberOfNodes()){
            local_constr_timer.restart();

            std::cout<<"fract: "<<fract<<" "<<n_samples<<" SAMPLES\n";
            NetworKit::EstimateBetweenness* rank = new NetworKit::EstimateBetweenness(graph,n_samples,false,true);


            rank->run();

            this->graph.forNodes([&] (vertex i){
                assert(graph.hasNode(i));
                assert(i<graph.numberOfNodes());
                ordering_rank[i]=std::make_pair(rank->score(i),i);

            });
            delete rank;
            cumulative_time+=local_constr_timer.elapsed();
            n_samples*=2;
        }
    }
    if(ordering_type == 3){
        INFO("BY DISTANCE-" + std::to_string((vertex) log2(V)) +" BOUNDED DOMINATING SET");
        std::set<vertex> to_cover;
        for(vertex v: graph.nodeRange()) to_cover.insert(v);
        vertex landmark_counter = 0;
        while(!to_cover.empty() && landmark_counter < L){
            std::queue<std::pair<vertex,dist>> que;
            auto it = to_cover.begin();
            std::advance(it,random() % to_cover.size());
            vertex source = *it;
            landmark_counter ++;
            que.push(std::make_pair(source, 0));
            while(!que.empty()){
                auto p = que.front();
                que.pop();
                if(p.second > ((vertex) log2(V))) break;
                ordering_rank[p.first] = std::make_pair(to_cover.size(), p.first);
                to_cover.erase(p.first);
                for(vertex v: graph.neighborRange(p.first)){
                    if(to_cover.find(v) != to_cover.end())
                        que.push(std::make_pair(v, p.second + 1));
                }
            }

        }
    }
    std::sort(ordering_rank, ordering_rank+graph.numberOfNodes(), [](const std::pair<double,vertex>  &a, const std::pair<double,vertex>  &b) {
        if(a.first == b.first)
            return a.second > b.second;
        else{
            return a.first > b.first;
        }
    });


    for(size_t count = 0; count < graph.numberOfNodes();count++){
        this->reverse_ordering[count]=ordering_rank[count].second;
        this->ordering[ordering_rank[count].second]=count;
    }

    delete[] ordering_rank;
}

void HighwayLabelling::GetBeerStores(std::vector<vertex> *beer_stores) {
    for(auto v: graph.nodeRange()){
        if(ordering[v] < L){
            beer_stores->push_back(v);
        }
    }
}

long HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < C[i]; j++) {
      if(query_distances[i][j] != null_distance)
        size++;
    }
  }

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      if(highway[i][j] != null_distance)
       size++;
    }
  }

  return size;
}

long HighwayLabelling::NaiveLabellingSize() {
    long size = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < L; j++) {
            if(naive_labeling[i][j] != null_distance)
                size++;
        }
    }

    return size;
}

void HighwayLabelling::ConstructUnweightedHighwayLabelling() {
  // Initialization
    distances = new dist*[V];
  for(int i = 0; i < V; i++) {
      distances[i] = new dist[L];
    for(int j = 0; j < L; j++)
      distances[i][j] = null_distance;
  }

  highway = new dist*[L];
  for(int i = 0; i < L; i++)
    highway[i] = new dist[L];

  // Start computing Highway Labelling (HL)
  ProgressStream hl_bar(L);

  hl_bar.label() << "Unweighted highway labeling construction";
  for (int i = 0; i < L; i++) {
    dist *P = new dist[V];
    for(int j = 0; j < V; j++)
      P[j] = null_distance;

    std::queue<int> que[2];

    que[0].push(reverse_ordering[i]); que[0].push(-1);
    distances[reverse_ordering[i]][i] = 0; P[reverse_ordering[i]] = 0; int use = 0;
    while (!que[0].empty()) {
      int u = que[use].front();
      que[use].pop();

      if(u == -1) {
        use = 1 - use;
        que[use].push(-1);
        continue;
      }

      for (vertex w : graph.neighborRange(u)) {
        if (P[w] == null_distance) {
          P[w] = P[u] + 1;
          if(use == 1 || ordering[w] < L)
            que[1].push(w);
          else {
            que[0].push(w);
            distances[w][i] = P[w];
          }
        }
      }
    }

    for(int j = 0; j < L; j++) {
      if(P[reverse_ordering[j]] != null_distance) {
        highway[i][j] = P[reverse_ordering[j]];
        highway[j][i] = P[reverse_ordering[j]];
      }
    }

    delete [] P;
    ++hl_bar;
  }
    C = new dist[V];
    vertices = new vertex*[V];
    query_distances = new dist*[V];

    for(vertex i = 0; i < V; i++) {
        C[i] = 0;
        std::vector<vertex> pos;
        for (vertex j = 0; j < L; j++) {
            if (distances[i][j] != null_distance) {
                C[i]++;
                pos.push_back(j);
            }
        }
        vertices[i] = new vertex[C[i]];
        query_distances[i] = new dist[C[i]];
        for (vertex j = 0; j < C[i]; j++){
            vertices[i][j] = pos[j];
            query_distances[i][j] = distances[i][pos[j]];
        }
        vertex vc = 0;
        for(vertex l = 0; l < L; l++){
            if(distances[i][l] != null_distance){
                if(distances[i][l] != query_distances[i][vc]){
                    throw new std::runtime_error("query distances error");
                }
                vc++;
            }
        }
        delete [] distances[i];
    }
    delete [] distances;
}

struct PQComparator
{
bool operator() (const std::pair<dist, vertex>& p1, const std::pair<dist, vertex>& p2)
    {
        return p1.first > p2.first;
    }
};

void HighwayLabelling::ConstructWeightedHighwayLabellingWithFlag() {
    // Initialization
    distances = new dist*[V];
    for(int i = 0; i < V; i++) {
        distances[i] = new dist[L];
        for(int j = 0; j < L; j++)
            distances[i][j] = null_distance;
    }

    highway = new dist*[L];
    for(int i = 0; i < L; i++)
        highway[i] = new dist[L];

    // Start computing Highway Labelling (HL)
    ProgressStream hl_bar(L);

    hl_bar.label() << "Weighted highway labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* pruning_flag = new bool[V];
        dist* dij_distances = new dist[V];
        vertex* predecessors = new vertex[V];
        for(vertex i = 0; i < V; i++){
            pruning_flag[i] = false;
            dij_distances[i] = null_distance;
            predecessors[i] = null_vertex;
        }
        highway[b][b] = 0;
        dij_distances[reverse_ordering[b]] = 0;
        distances[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();
            if(ordering[v] < L && ordering[v] != b){
                pruning_flag[v] = true;
                highway[b][ordering[v]] = dij_distances[v];
                highway[ordering[v]][b] = dij_distances[v];
            }
            if(!pruning_flag[v]){
                distances[v][b] = dij_distances[v];
            }
            for(auto w: graph.neighborRange(v)){
                if(dij_distances[w] > dij_distances[v] + graph.weight(v,w)){
                    dij_distances[w] = dij_distances[v] + (dist)graph.weight(v,w);
                    predecessors[w] = v;
                    pruning_flag[w] = pruning_flag[v];
                    pq.push(std::make_pair(dij_distances[w], w));
                }

            }
        }
        delete [] pruning_flag;
        delete [] dij_distances;
        delete [] predecessors;
        ++hl_bar;
    }
    C = new dist[V];
    vertices = new vertex*[V];
    query_distances = new dist*[V];

    for(vertex i = 0; i < V; i++) {
        C[i] = 0;
        std::vector<vertex> pos;
        for (vertex j = 0; j < L; j++) {
            if (distances[i][j] != null_distance) {
                C[i]++;
                pos.push_back(j);
            }
        }
        vertices[i] = new vertex[C[i]];
        query_distances[i] = new dist[C[i]];
        for (vertex j = 0; j < C[i]; j++){
            vertices[i][j] = pos[j];
            query_distances[i][j] = distances[i][pos[j]];
        }
        vertex vc = 0;
        for(vertex l = 0; l < L; l++){
            if(distances[i][l] != null_distance){
                if(distances[i][l] != query_distances[i][vc]){
                    throw new std::runtime_error("query distances error");
                }
                vc++;
            }
        }
        delete [] distances[i];
    }
    delete [] distances;
}

void HighwayLabelling::ConstructUnweightedNaiveLabelling() {
    // Initialization
    naive_labeling = new dist*[V];
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i] = new dist[L];
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }


    // Start computing Naive Labelling
    ProgressStream na_bar(L);

    na_bar.label() << "Unweighted naive labeling construction";
    for(vertex s = 0; s < L; s++){
        dist *P = new dist[V];
        for(vertex j = 0; j < V; j++)
            P[j] = null_distance;
        std::queue<vertex> que;
        que.push(reverse_ordering[s]);
        naive_labeling[reverse_ordering[s]][s] = 0; P[reverse_ordering[s]] = 0;

        while(!que.empty()){
            vertex v = que.front();
            que.pop();
            naive_labeling[v][s] = P[v];

            for(vertex w: graph.neighborRange(v)){
                if(P[w]==null_distance){
                    P[w] = P[v] + 1;
                    que.push(w);
                }
            }
        }
    ++na_bar;
    }

}

void HighwayLabelling::ConstructWeightedNaiveLabelling() {
    // Initialization
    naive_labeling = new dist*[V];
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i] = new dist[L];
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }
    ProgressStream na_bar(L);

    na_bar.label() << "Weighted naive labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* pruning_flag = new bool[V];

        naive_labeling[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();

            for(auto w: graph.neighborRange(v)){
                if(naive_labeling[w][b] > naive_labeling[v][b] + graph.weight(v,w)){
                    naive_labeling[w][b] = naive_labeling[v][b] + (dist)graph.weight(v,w);

                    pruning_flag[w] = pruning_flag[v];
                    pq.push(std::make_pair(naive_labeling[w][b], w));
                }

            }
        }
        ++na_bar;
        delete [] pruning_flag;
    }
    
}

dist HighwayLabelling::min(dist a, dist b) {
  return (a < b) ? a : b;
}

dist HighwayLabelling::QueryDistance(vertex s, vertex t) {

  dist m = null_distance, uni1[C[s]] = {0}, uni2[C[t]] = {0}; vertex i = 0, j = 0, i1 = 0, j1 = 0;
  while (i < C[s] && j < C[t]) {
    if (vertices[s][i] < vertices[t][j]) {
      uni1[i1] = i; i++; i1++;
    } else if (vertices[t][j] < vertices[s][i]) {
      uni2[j1] = j; j++; j1++;
    } else {
      m = min(m, query_distances[s][i] + query_distances[t][j]);
      i++; j++;
    }
  }

  while (i < C[s]) {
    uni1[i1] = i; i++; i1++;
  }

  while (j < C[t]) {
    uni2[j1] = j; j++; j1++;
  }

  i = 0;
  while (i < i1) { j = 0;
    while (j < j1) {
      m = min(m, query_distances[s][uni1[i]] + highway[vertices[s][uni1[i]]][vertices[t][uni2[j]]] + query_distances[t][uni2[j]]);
      j++;
    }
    i++;
  }

  return m;
}

dist HighwayLabelling::QueryDistanceQuadratic(vertex s, vertex t) {
    dist m = null_distance;
    vertex i, j;
    for(i = 0; i < L; i++) {
        for (j = 0; j < L; j++)
            m = min(m, distances[s][i] + highway[vertices[s][i]][vertices[t][j]] + distances[t][j]);
    }

    return m;
}

dist HighwayLabelling::NaiveQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    for(vertex i = 0; i < L; i++){
            m = min(m, naive_labeling[s][i]+naive_labeling[t][i]);
    }
    return m;
}

dist HighwayLabelling::BFSQuery(vertex s, vertex t) {
    // Initialization
    dist *s_to_vertices = new dist[V];
    dist *t_to_vertices = new dist[V];
    for(vertex j = 0; j < V; j++) {
        s_to_vertices[j] = null_distance;
        t_to_vertices[j] = null_distance;
    }
    std::queue<vertex> que;
    que.push(s);
    s_to_vertices[s] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(s_to_vertices[w]==null_distance){
                s_to_vertices[w] = s_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    que.push(t);
    t_to_vertices[t] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(t_to_vertices[w]==null_distance){
                t_to_vertices[w] = t_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    dist m = null_distance;
    for(vertex l = 0; l < L; l++){
        m = min(s_to_vertices[reverse_ordering[l]] + t_to_vertices[reverse_ordering[l]], m);
    }

    delete [] s_to_vertices;
    delete [] t_to_vertices;
    return m;
}

dist HighwayLabelling::DijkstraQuery(vertex s, vertex t) {
    dist* s_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        s_distances[i] = null_distance;
    }
    s_distances[s] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0,s));
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        for(auto w: graph.neighborRange(v)){
            if(s_distances[w] > s_distances[v] + graph.weight(v,w)){
                s_distances[w] = s_distances[v] + (dist)graph.weight(v,w);

                pq.push(std::make_pair(s_distances[w], w));
            }

        }
    }

    dist* t_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        t_distances[i] = null_distance;
    }
    t_distances[t] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pqt;
    pqt.push(std::make_pair(0,t));
    while(!pqt.empty()){
        vertex v = pqt.top().second;
        pqt.pop();
        for(auto w: graph.neighborRange(v)){
            if(t_distances[w] > t_distances[v] + graph.weight(v,w)){
                t_distances[w] = t_distances[v] + (dist)graph.weight(v,w);

                pqt.push(std::make_pair(t_distances[w], w));
            }

        }
    }

//    auto s_to_vertices = NetworKit::Dijkstra(graph, s, false, false);
//    s_to_vertices.run();
//    auto t_to_vertices = NetworKit::Dijkstra(graph, t, false, false);
//    t_to_vertices.run();
//    dist m = null_distance;
//    for(vertex l = 0; l < L; l++){
//        m = min(s_to_vertices.distance(reverse_ordering[l]) + t_to_vertices.distance(reverse_ordering[l]), m);
//    }
    dist m = null_distance;
    for(vertex l = 0; l < L; l++){
        m = min(m, s_distances[reverse_ordering[l]] + t_distances[reverse_ordering[l]]);
    }
    delete [] s_distances;
    delete [] t_distances;
    return m;
}

dist HighwayLabelling::P2PBFS(vertex s, vertex landmark, dist ub, dist partial, dist* to_vertices){
    std::vector<vertex> visited;
    std::queue<vertex> que;
    que.push(s);
    to_vertices[s] = 0;
    visited.push_back(s);

    while(!que.empty()){
        vertex v = que.front();
        que.pop();
        if(to_vertices[v] >= ub || to_vertices[v] + partial >= ub){
            for(auto u: visited)
                to_vertices[u] = null_distance;
            return null_distance;
        }
        if(v == landmark){
            dist m = to_vertices[v];
            for(auto u: visited)
                to_vertices[u] = null_distance;
            return m + partial;
        }
        for(vertex w: graph.neighborRange(v)){
            if(to_vertices[w]==null_distance){
                to_vertices[w] = to_vertices[v] + 1;
                visited.push_back(w);

                que.push(w);
            }
        }
    }

    for(auto u: visited)
        to_vertices[u] = null_distance;
    return null_distance;
}

dist HighwayLabelling::BoundedSearch(vertex s, vertex t) {
    // Initialization
    dist *to_vertices = new dist[V];
    for(vertex j = 0; j < V; j++) {
        to_vertices[j] = null_distance;
    }
    dist m = null_distance;
    for(vertex v = 0; v < L; v++){
        dist s_to_v = P2PBFS(s, reverse_ordering[v], m, 0, to_vertices);
        if(s_to_v == null_distance)
            continue;
        dist s_to_t = P2PBFS(t, reverse_ordering[v], m, s_to_v, to_vertices);
        m = min(m, s_to_t);
    }

    delete [] to_vertices;
    return m;
}

#endif  // HGHWAY_LABELING_H_
