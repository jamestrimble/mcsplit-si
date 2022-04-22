#include <limits.h>

#include <vector>

#define BYTES_PER_WORD sizeof(unsigned long long)
#define BITS_PER_WORD (CHAR_BIT * BYTES_PER_WORD)

struct SparseGraph {
    int n;
    std::vector<std::vector<int>> adj_lists;
    std::vector<std::vector<int>> in_edge_lists;
    std::vector<std::vector<int>> filtered_adj_lists;
    std::vector<unsigned int> label;
    SparseGraph(unsigned int n);
    bool has_edge(int v, int w) const;
    bool has_filtered_edge(int v, int w) const;
    int bitset_word_count() const { return (n+BITS_PER_WORD-1)/BITS_PER_WORD; }
};

SparseGraph induced_subgraph(struct SparseGraph& g, std::vector<int> vv);

SparseGraph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);
