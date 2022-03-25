#include <limits.h>
#include <stdbool.h>

#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<std::vector<int>> adj_lists;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

void make_adj_lists(Graph & g, const std::vector<bool> & active_vertices);
