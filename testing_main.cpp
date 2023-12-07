#include <iostream>
#include "highway_cover_labelling.h"
#include <vector>
#include <boost/program_options.hpp>
#include "networkit/components/ConnectedComponents.hpp"
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/graph/GraphTools.hpp>
#include "mytimer.h"
#include <networkit/distance/Diameter.hpp>

using namespace std;
double median(std::vector<double>& arr) { //SORTS
    size_t n = arr.size() / 2;
    if (n % 2 == 0) {
        std::nth_element(arr.begin(),arr.begin() + n/2,arr.end());
        std::nth_element(arr.begin(),arr.begin() + (n - 1) / 2,arr.end());
        return (double) (arr[(n-1)/2]+ arr[n/2])/2.0;
    }

    else{
        std::nth_element(arr.begin(),arr.begin() + n / 2,arr.end());
        return (double) arr[n/2];
    }
    assert(false);
}



double average(std::vector<double> & arr) {

    auto const count = static_cast<double>(arr.size());
    double sum = 0;
    for(double value: arr) sum += value;
    return sum / count;
}

int main(int argc, char **argv) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");

    desc.add_options()
            ("graph_location,g", po::value<std::string>(), "Input Graph File Location")
            ("landmarks,l", po::value<int>(), "Number of Landmarks (beer shop)")
            ("num_queries,q", po::value<int>(), "Number of Queries to Be Performed")
            ("ordering,o",po::value<int>(), "Type of Node Ordering [DEGREE(0) RANDOM(1)]")
            ("index_file,i",po::value<string>(), "Name of the file that will store the index")
            ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if(vm.empty()){
        std::cout << desc << "\n";
        throw std::runtime_error("empty argument array");
    }
    std::string graph_location;

    if(vm.count("graph_location")){
        graph_location = vm["graph_location"].as<std::string>();
    }


    int L = -1;

    if (vm.count("landmarks")){
        L = vm["landmarks"].as<int>();
    }

    if (L < 1){
        std::cout << desc << "\n";
        throw std::runtime_error("L must be at least 1");
    }

    int ordering = -1;

    if (vm.count("ordering")){
        ordering = vm["ordering"].as<int>();
    }
    if(ordering != 0 && ordering != 1 && ordering != 2 && ordering != 3){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong ordering selection (0 or 1 or 2)");
    }

    int num_queries = -1;
    if (vm.count("num_queries")){
        num_queries = vm["num_queries"].as<int>();
    }
    if(num_queries < 2){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong num queries");
    }

   std::string index_file;
    if(vm.count("index_file")){
        index_file = vm["index_file"].as<std::string>();
    }
    std::cout << "Reading " << graph_location << " building with L = " << L << " num_queries = "<<num_queries<<" ordering "<<ordering<<" index_file " << index_file << "\n";
    NetworKit::EdgeListReader edl(' ', 1);
    NetworKit::Graph graph = edl.read(graph_location);
    const NetworKit::Graph &graph_handle = graph;
    NetworKit::ConnectedComponents *cc = new NetworKit::ConnectedComponents(graph_handle);
    graph = cc->extractLargestConnectedComponent(graph_handle, true);
    graph.shrinkToFit();
    graph.indexEdges();
    std::cout << "Graph after CC has " << graph.numberOfNodes() << " vertices and " << graph.numberOfEdges()
              << " edges\n";
//    double density = NetworKit::GraphTools::density(graph);
//
//    NetworKit::Diameter *dm = new NetworKit::Diameter(graph);
//    dm->run();
//
//
//   double diameter = dm->getDiameter().first;
//    delete dm;
//   std::cout << "Density: " << density << "\n";
//    std::cout << "Diameter: " << diameter << "\n";
    vertex diameter = 0;
    HighwayLabelling *hl = new HighwayLabelling(graph, L, ordering);
    mytimer timer;
    int topl[L];
    sort(topl, topl + L);
    // construct labelling
    double hl_constr_time = 0.0;
    timer.restart();
    if(graph.isWeighted())
        hl->ConstructWeightedHighwayLabellingWithFlag();
    else
        hl->ConstructUnweightedHighwayLabelling();
    hl_constr_time = timer.elapsed();
    long hl_size = hl->LabellingSize();
    cout << "HL constr time " << hl_constr_time << " | HL size " << hl_size << "\n";
    //hl->StoreIndex(index_file);

    //delete hl;
    //auto *hl2 = new HighwayLabelling(graph, L, ordering);
    //hl2->LoadIndex(index_file);
    double naive_cnstr_time = 0.0;
    timer.restart();
    if(graph.isWeighted())
        hl->ConstructWeightedNaiveLabelling();
    else
        hl->ConstructUnweightedNaiveLabelling();
    naive_cnstr_time += timer.elapsed();
    long hl_naive_size = hl->NaiveLabellingSize();
    cout << "Naive constr time " << naive_cnstr_time << " | Naive size " << hl_naive_size << "\n";
    vector<double> hl_query_time;
    vector<double> naive_query_time;
    vector<double> sssp_query_time;
    vector<double> bs_query_time;
    long q = num_queries;
    ProgressStream query_bar(q);

    query_bar.label() << "Query computation";
    while(num_queries--){
        int s = NetworKit::GraphTools::randomNode(graph);
        int t = NetworKit::GraphTools::randomNode(graph);
//        int s = 4;
//        int t = 0;
        double hl_q_time = 0.0;
        timer.restart();
        dist  hld= hl->QueryDistance(s,t);
        hl_q_time += timer.elapsed();
        hl_query_time.push_back(hl_q_time);
        double naive_q_time = 0.0;
        timer.restart();
        dist nad = hl->NaiveQueryDistance(s,t);
        naive_q_time += timer.elapsed();
        naive_query_time.push_back(naive_q_time);
        dist ssspd;
        double sssp_time = 0.0;
        timer.restart();
        if(graph.isWeighted())
             ssspd = hl->DijkstraQuery(s, t);
        else
            ssspd = hl->BFSQuery(s, t);
        sssp_time += timer.elapsed();
        sssp_query_time.push_back(sssp_time);
        if(ssspd != hld){
            cout << "Error between " << s << " and " << t << "\n";
            cout << "HL distance " << (int)hld << " | SSSP distance " << (int)ssspd << "\n";
            throw new std::runtime_error("experiment fails");
        }

//            double bs_time = 0.0;
//            timer.restart();
//            uint8_t bsd = hl2->BoundedSearch(s, t);
//            bs_time += timer.elapsed();
//            bs_query_time.push_back(bs_time);
//            if(bsd != hld){
//                cout << "Error between " << s << " and " << t << "\n";
//                cout << "HL distance " << (int)hld << " | BiDir distance " << (int)bsd << "\n";
//                throw new std::runtime_error("experiment fails");
//            }
        if(nad != hld){
            cout << "Error between " << s << " and " << t << "\n";
            cout << "HL distance " << (int)hld << " | Naive distance " << (int)nad << "\n";
            throw new std::runtime_error("experiment fails");
        }

        ++query_bar;
    }
    cout << "HL median time " << median(hl_query_time) << " | mean time " << average(hl_query_time) << "\n";
    cout << "Naive median time " << median(naive_query_time) << " | mean time " << average(naive_query_time) << "\n";
    cout << "SSSP median time " << median(sssp_query_time) << " | mean time " << average(sssp_query_time) << "\n";
    //cout << "BoundedSearch median time " << median(bs_query_time) << " | mean time " << average(bs_query_time) << "\n";

    std::string order_string;

    switch(ordering){
        case (0):
            order_string = "DEG";
            break;
        case (1):
            order_string = "RND";
            break;
        case (2):
            order_string = "BET";
            break;
        case (3):
            order_string = "DTS";
            break;
        default:
            throw new std::runtime_error("problem on order string");


    }

    std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
    std::string tmp_time = buffer;
    std::string shortenedName = graph_location.substr(0, 16);
    stringstream string_object_name;
    string_object_name<<L;
    std::string l_string;
    string_object_name>>l_string;
    std::string timestampstring = "custom_beer_"+shortenedName+"_"+l_string+"_"+"_"+order_string+"_"+"_"+tmp_time;

    std::string logFile = timestampstring +".csv";


    std::ofstream ofs(logFile);

    ofs << "G,V,E,L,ORD,Q,DIAM,HLCT,HLIS,HLQT,NACT,NAIS,NAQT,NPQT,BDQT\n";
    ofs << graph_location << "," << graph.numberOfNodes() << "," << graph.numberOfEdges() << "," << L << "," <<
        order_string << "," << q << "," << diameter << "," << hl_constr_time << "," << hl_size << "," <<
        average(hl_query_time) << "," << naive_cnstr_time << "," << hl_naive_size << "," << average(naive_query_time) <<
        "," << average(sssp_query_time) << "," << 0 << "\n";
    ofs.close();
    delete hl;
    exit(EXIT_FAILURE);
}
