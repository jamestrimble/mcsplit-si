#include "sparse_graph.h"

#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <list>

#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { heur_A, heur_B, heur_C, heur_RI };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Subgraph isomorphism\vHEURISTIC can be A or B";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"vertex-labelled-only", 'x', 0, 0, "Use vertex labels, but not edge labels"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    {"deg-heur", 'D', "heur", 0, "Degree heuristic (0, 1, or 2)"},
    { 0 }
};

static struct {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool gfd;
    bool vf;
    bool edge_labelled;
    bool vertex_labelled;
    Heuristic heuristic;
    char *filename1;
    char *filename2;
    int timeout;
    int deg_heur;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.edge_labelled = false;  // TODO: remove (unused)
    arguments.vertex_labelled = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.arg_num = 0;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'x':
            arguments.vertex_labelled = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case 'D':
            arguments.deg_heur = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "A")
                    arguments.heuristic = heur_A;
                else if (std::string(arg) == "B")
                    arguments.heuristic = heur_B;
                else if (std::string(arg) == "C")
                    arguments.heuristic = heur_C;
                else if (std::string(arg) == "RI")
                    arguments.heuristic = heur_RI;
                else
                    fail("Unknown heuristic (try A, B, C, or RI)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes{ 0 };

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    bool active;
    int mod_index;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent),
            active(true),
            mod_index(INT_MAX)
            { };
    Bidomain(): Bidomain(0, 0, 0, 0, false) { };
};

using It = vector<int>::iterator;

struct NewBidomain;

//using BDLL = std::list<NewBidomain>;

using BdIt = NewBidomain *;  //std::list<NewBidomain>::iterator;

struct NewBidomain {
    It l, r;
    It l_end, r_end;
    int l_adj_count, r_adj_count;   // used when partitioning
    union {
        NewBidomain *next_in_split_list;
        NewBidomain *next_in_free_list;
    };
    NewBidomain *next_in_deleted_list;
    NewBidomain *prev;
    NewBidomain *next;
    bool active;
    bool undergoing_split;

    void initialise(It l, It r, It l_end, It r_end)
    {
        this->l = l;
        this->r = r;
        this->l_end = l_end;
        this->r_end = r_end;
        this->active = true;
        this->undergoing_split = false;
    }

    void insert_before(BdIt p)
    {
        BdIt new_prev = p->prev;
        this->prev = new_prev;
        this->next = p;
        p->prev = this;
        new_prev->next = this;
    }
    void insert_after(BdIt p)
    {
        BdIt new_next = p->next;
        this->next = new_next;
        this->prev = p;
        p->next = this;
        new_next->prev = this;
    }
    void move_to_before(BdIt p)
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
        this->insert_before(p);
    }
    void move_to_after(BdIt p)
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
        this->insert_after(p);
    }
    void remove()
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
    }
    void reinsert()
    {
        this->prev->next = this;
        this->next->prev = this;
    }
    int l_size() const { return l_end - l; }
    int r_size() const { return r_end - r; }
};

struct BDLL {
    // A circular doubly-linked list
    // (It's just circular to avoid the need for separate head and tail)
    NewBidomain head;
    BDLL() {
        this->head.next = &this->head;
        this->head.prev = &this->head;
    }

    NewBidomain *begin() { return head.next; }
    NewBidomain *end() { return &head; }
    NewBidomain & back() { return *head.prev; }
    int size()
    {
        int result = 0;
        for (BdIt bd_it=begin(); bd_it!=end(); bd_it=bd_it->next) {
            ++result;
        }
        return result;
    }
    bool empty() { return head.next == &head; }
};


// Some temporary storage space
struct Workspace {
    vector<BdIt> split_bds;
    NewBidomain *get_from_free_list()
    {
#ifdef WITHOUT_OPTIMISATIONS
        return new NewBidomain;
#endif
        if (bd_free_list == nullptr) {
            bd_memory_pools.push_back(vector<NewBidomain>(100));
            for (NewBidomain & bd : bd_memory_pools.back()) {
                bd.next_in_free_list = bd_free_list;
                bd_free_list = &bd;
            }
        }
        NewBidomain *bd = bd_free_list;
        bd_free_list = bd->next_in_free_list;
        return bd;
    }
    void add_to_free_list(NewBidomain * bd)
    {
#ifdef WITHOUT_OPTIMISATIONS
        delete bd;
        return;
#endif
        bd->next_in_free_list = bd_free_list;
        bd_free_list = bd;
    }
    vector<int> left_swap_vv;
    vector<int> right_swap_vv;
    vector<int> vv0;
    vector<int> vv1;
    vector<int> vv0_inverse;
    vector<int> vv1_inverse;
    vector<int> order;
    Workspace(vector<int> & vv0, vector<int> & vv1, vector<int> order) :
            vv0(vv0), vv1(vv1), vv0_inverse(vv0.size()), vv1_inverse(vv1.size()), order(order)
    {
        for (unsigned i=0; i<vv0.size(); i++) {
            vv0_inverse[vv0[i]] = i;
        }
        for (unsigned i=0; i<vv1.size(); i++) {
            vv1_inverse[vv1[i]] = i;
        }
    }
private:
    NewBidomain *bd_free_list = nullptr;
    std::list<vector<NewBidomain>> bd_memory_pools;
};

struct Ptrs
{
    BdIt bd_it;
    It vtx_it;
};

void show(const vector<Ptrs> & left_ptrs, const vector<Ptrs> & right_ptrs,
        const SparseGraph & g0, const SparseGraph & g1, const vector<VtxPair>& current,
        BDLL & bdll)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    if (!current.empty() > 0) {
        std::cout << "Adjacent to " << current.back().v << " in g0: ";
        for (int x : g0.adj_lists[current.back().v]) std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Adjacent to " << current.back().w << " in g1: ";
        for (int x : g1.adj_lists[current.back().w]) std::cout << x << " ";
        std::cout << std::endl;
    }
    cout << "---------------------" << std::endl;
    for (BdIt bd_it=bdll.begin(); bd_it!=bdll.end(); bd_it=bd_it->next) {
        cout << "Left  ";
        for (It it=bd_it->l; it!=bd_it->l_end; it++)
            cout << *it << " ";
        cout << std::endl;
        cout << "Right  ";
        for (It it=bd_it->r; it!=bd_it->r_end; it++)
            cout << *it << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
    cout << "\n" << std::endl;
}

bool check_sol(const SparseGraph & g0, const SparseGraph & g1 , const vector<VtxPair> & solution) {
    return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if (g0.has_edge(p0.v, p1.v) != g1.has_edge(p0.w, p1.w))
                return false;
        }
    }
    return true;
}

// Return the position in vv0 of the pattern vertex on which we'd like to branch from the range [begin,end)
int find_best_v_score(It begin, It end, Workspace & workspace)
{
    int result_score = workspace.vv0_inverse[*begin];
    for (It it=std::next(begin); it!=end; it++) {
        int v_score = workspace.vv0_inverse[*it];
        if (v_score < result_score) {
            result_score = v_score;
        }
    }
    return result_score;
}

// Return the pattern vertex on which we'd like to branch from the range [begin,end)
int find_best_v(It begin, It end, Workspace & workspace)
{
    return workspace.vv0[find_best_v_score(begin, end, workspace)];
}

int select_branch_v_heur_A(BDLL & domains, Workspace & workspace)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int size = std::max(bd.l_size(), bd.r_size());
        if (size < min_size)
            min_size = size;
    }
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int size = std::max(bd.l_size(), bd.r_size());
        if (size != min_size)
            continue;
        int tie_breaker = find_best_v_score(bd.l, bd.l_end, workspace);
        if (tie_breaker < min_tie_breaker) {
            min_tie_breaker = tie_breaker;
            best = workspace.vv0[tie_breaker];
        }
    }
    return best;
}

int select_branch_v_heur_B(BDLL & domains, const SparseGraph & g0, Workspace & workspace)
{
    throw 1;
    double best_score = INT_MIN;
    int best = -1;
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int best_v = find_best_v(bd.l, bd.l_end, workspace);
        int right_len = bd.r_end - bd.r;
        if (right_len == 1) {
            // Special case where no branching is required
            return best_v;
        }
        int deg = g0.adj_lists[best_v].size();
        double score = double(deg) / right_len;
        if (score > best_score) {
            best_score = score;
            best = best_v;
        }
    }
    return best;
}

int select_branch_v_heur_C(BDLL & domains, const SparseGraph & g0, const vector<int> & g0_remaining_deg, Workspace & workspace)
{
    throw 1;
    double best_score = INT_MIN;
    int best = -1;
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int best_v = find_best_v(bd.l, bd.l_end, workspace);
        int right_len = bd.r_end - bd.r;
        if (right_len == 1) {
            // Special case where no branching is required
            return best_v;
        }
        int remaining_deg = g0_remaining_deg[best_v];
        double score = double(remaining_deg) / right_len;
        if (score > best_score) {
            best_score = score;
            best = best_v;
        }
    }
    return best;
}

int select_branch_v(BDLL & domains, const SparseGraph & g0, const vector<int> & g0_remaining_deg,
        Workspace & workspace, int current_len)
{
    if (arguments.heuristic == heur_A)
        return select_branch_v_heur_A(domains, workspace);
    else if (arguments.heuristic == heur_B)
        return select_branch_v_heur_B(domains, g0, workspace);
    else if (arguments.heuristic == heur_C)
        return select_branch_v_heur_C(domains, g0, g0_remaining_deg, workspace);
    else
        return workspace.order[current_len];
}

struct SplitAndDeletedLists {
    bool quit_early;
    unsigned new_bound;
    NewBidomain *split_bds_list;
    NewBidomain *deleted_bds_list;
    SplitAndDeletedLists(bool quit_early, unsigned new_bound, NewBidomain *split_bds_list, NewBidomain *deleted_bds_list)
        : quit_early(quit_early), new_bound(new_bound), split_bds_list(split_bds_list), deleted_bds_list(deleted_bds_list) {}
};

void partition_left(const vector<int> & left_vv, vector<Ptrs> & left_ptrs,
        vector<BdIt> & split_bds, vector<int> & left_swap_vv)
{
    for (int u : left_vv) {
        auto bd_it = left_ptrs[u].bd_it;
        if (bd_it == nullptr) continue;
        if (!bd_it->active) continue;
        if (!bd_it->undergoing_split) {
            bd_it->l_adj_count = 0;
            bd_it->r_adj_count = 0;
            bd_it->undergoing_split = true;
            split_bds.push_back(bd_it);
        }
        ++bd_it->l_adj_count;
        left_swap_vv.push_back(u);
    }
}

void partition_right(const vector<int> & right_vv, vector<Ptrs> & right_ptrs,
        vector<BdIt> & split_bds, vector<int> & right_swap_vv)
{
    for (int u : right_vv) {
        auto bd_it = right_ptrs[u].bd_it;
        if (bd_it == nullptr) continue;
        if (!bd_it->active) continue;
        if (!bd_it->undergoing_split) {
            bd_it->l_adj_count = 0;
            bd_it->r_adj_count = 0;
            bd_it->undergoing_split = true;
            split_bds.push_back(bd_it);
        }
        ++bd_it->r_adj_count;
        right_swap_vv.push_back(u);
    }
}

NewBidomain * do_splits(Workspace & workspace, vector<BdIt> & split_bds,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs)
{
    NewBidomain *split_bds_list = nullptr;
    for (auto bd_it : split_bds) {
        BdIt new_elem = workspace.get_from_free_list();
        new_elem->insert_after(bd_it);
        new_elem->initialise(bd_it->l_end - bd_it->l_adj_count,
                bd_it->r_end - bd_it->r_adj_count, bd_it->l_end, bd_it->r_end);

        // Insert the new BD at the head of the linked list of split BDs
        new_elem->next_in_split_list = split_bds_list;
        split_bds_list = new_elem;
    }

    for (int u : workspace.left_swap_vv) {
        auto & u_ptrs = left_ptrs[u];
        auto & bd = *u_ptrs.bd_it;
        It u_it = u_ptrs.vtx_it;
        It m = --bd.l_end;
        int t = *m;
        *u_it = t;
        *m = u;
        u_ptrs.vtx_it = m;
        u_ptrs.bd_it = bd.next;
        left_ptrs[t].vtx_it = u_it;
    }

    for (int u : workspace.right_swap_vv) {
        auto & u_ptrs = right_ptrs[u];
        auto & bd = *u_ptrs.bd_it;
        It u_it = u_ptrs.vtx_it;
        It m = --bd.r_end;
        int t = *m;
        *u_it = t;
        *m = u;
        u_ptrs.vtx_it = m;
        u_ptrs.bd_it = bd.next;
        right_ptrs[t].vtx_it = u_it;
    }

    return split_bds_list;
}

NewBidomain * do_deletions(vector<BdIt> & split_bds)
{
    NewBidomain *deleted_bds_list = nullptr;
    for (auto bd_it : split_bds) {
        for (int i=0; i<2; i++) {
            // Delete old and new BDs if necessary
            auto & bd = *bd_it;
            if (bd.l == bd.l_end || bd.r == bd.r_end) {
                bd.active = false;

                // add to deleted list
                bd.remove();
                bd.next_in_deleted_list = deleted_bds_list;
                deleted_bds_list = bd_it;
            }
            bd_it = bd_it->next;
        }
    }
    return deleted_bds_list;
}

SplitAndDeletedLists filter_domains(
        Workspace & workspace,
        BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        const vector<vector<int>> & g0_adj_lists, const vector<vector<int>> & g1_adj_lists,
        int v, int w, unsigned matching_size_goal, unsigned bound)
{
    vector<BdIt> & split_bds = workspace.split_bds;
    split_bds.clear();

    workspace.left_swap_vv.clear();
    workspace.right_swap_vv.clear();

    partition_left(g0_adj_lists[v], left_ptrs, split_bds, workspace.left_swap_vv);

    partition_right(g1_adj_lists[w], right_ptrs, split_bds, workspace.right_swap_vv);

    for (auto bd_it : split_bds) {
        bd_it->undergoing_split = false;
    }

    unsigned new_bound = bound;

    for (auto bd_it : split_bds) {
        int l0_size = bd_it->l_size() - bd_it->l_adj_count;
        int r0_size = bd_it->r_size() - bd_it->r_adj_count;
        int l1_size = bd_it->l_adj_count;
        int r1_size = bd_it->r_adj_count;
        unsigned old_bd_bound = std::min(l0_size + l1_size, r0_size + r1_size);
        unsigned new_bd_bound = std::min(l0_size, r0_size) + std::min(l1_size, r1_size);
        new_bound = new_bound - old_bd_bound + new_bd_bound;
        if (new_bound < matching_size_goal) {
            return { true, 0, nullptr, nullptr };
        }
    }

    NewBidomain *split_bds_list = do_splits(workspace, split_bds, left_ptrs, right_ptrs);

    NewBidomain *deleted_bds_list = do_deletions(split_bds);

    return {false, new_bound, split_bds_list, deleted_bds_list};
}

void unfilter_domains(
        Workspace & workspace,
        BDLL & bdll,
        vector<Ptrs> & left_ptrs,
        vector<Ptrs> & right_ptrs,
        NewBidomain *split_bds_list,
        NewBidomain *deleted_bds)
{
    while (deleted_bds != nullptr) {
        NewBidomain & bd = *deleted_bds;
        bd.active = true;
        bd.reinsert();
        deleted_bds = bd.next_in_deleted_list;
        //deleted_bds.back().move_to_before(bd.reinsertion_point);
        //bdll.splice(bd.reinsertion_point, deleted_bds, std::prev(deleted_bds.end()));
    }

    for (NewBidomain *p=split_bds_list; p!=nullptr; ) {
        // TODO: better variable names
        //BdIt nxt_it = p;
        auto & nxt = *p;
        BdIt bd_it = nxt.prev;
        for (It it=nxt.l; it<nxt.l_end; it++) {
            left_ptrs[*it].bd_it = bd_it;
        }
        for (It it=nxt.r; it<nxt.r_end; it++) {
            right_ptrs[*it].bd_it = bd_it;
        }
        bd_it->l_end = nxt.l_end;
        bd_it->r_end = nxt.r_end;

        p->remove();
        p = p->next_in_split_list;

        workspace.add_to_free_list(&nxt);
        //workspace.bd_free_list.splice(workspace.bd_free_list.end(), bdll, nxt_it);
        //bdll.erase(std::next(bd_it));
    }
}

void assign(int v, int w,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs)
{
    auto & bd = *left_ptrs[v].bd_it;
    It v_it = left_ptrs[v].vtx_it;
    It l_last = std::prev(bd.l_end);
    std::swap(*v_it, *l_last);
    left_ptrs[v].vtx_it = l_last;
    left_ptrs[v].bd_it = nullptr;
    left_ptrs[*v_it].vtx_it = v_it;
    --bd.l_end;

    It w_it = right_ptrs[w].vtx_it;
    It r_last = std::prev(bd.r_end);
    std::swap(*w_it, *r_last);
    right_ptrs[w].vtx_it = r_last;
    right_ptrs[w].bd_it = nullptr;
    right_ptrs[*w_it].vtx_it = w_it;
    --bd.r_end;
}

void unassign(int v, int w, BdIt bd_it,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs)
{
    left_ptrs[v].bd_it = bd_it;
    right_ptrs[w].bd_it = bd_it;
//    cout << (bd.l_end - bd.l) << "------------" << endl;
    ++bd_it->l_end;
    ++bd_it->r_end;
}

void delete_v(int v, vector<Ptrs> & left_ptrs)
{
    auto & bd = *left_ptrs[v].bd_it;
    It v_it = left_ptrs[v].vtx_it;
    It l_last = std::prev(bd.l_end);
    std::swap(*v_it, *l_last);
    left_ptrs[v].vtx_it = l_last;
    left_ptrs[v].bd_it = nullptr;
    left_ptrs[*v_it].vtx_it = v_it;
    --bd.l_end;
}

void undelete_v(int v, BdIt bd_it, vector<Ptrs> & left_ptrs)
{
    left_ptrs[v].bd_it = bd_it;
    ++bd_it->l_end;
}

void solve(Workspace & workspace, const SparseGraph & g0, const SparseGraph & g1, vector<VtxPair> & incumbent,
        vector<VtxPair> & current, BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        long long & solution_count, vector<int> & g0_remaining_deg, unsigned matching_size_goal, unsigned bound)
{
    if (abort_due_to_timeout)
        return;

//    std::cout << bound << std::endl;
    if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll);
    nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        //if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }

    if (current.size()==matching_size_goal) {
        solution_count++;
        return;
    }

    int v = select_branch_v(bdll, g0, g0_remaining_deg, workspace, current.size());
        
    if (v == -1)
        return;

    BdIt bd_it = left_ptrs[v].bd_it;

    std::vector<int> ww;
    ww.reserve(bd_it->r_end - bd_it->r);
    for (auto it = bd_it->r; it != bd_it->r_end; it++) {
        // Initially ww contains vertices' positions in the order.  After sorting,
        // we'll convert these to vertices
        ww.push_back(workspace.vv1_inverse[*it]);
    }
    std::sort(ww.begin(), ww.end());
    for (auto & w : ww) {
        w = workspace.vv1[w];
    }

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    for (int w : ww) {
        assign(v, w, left_ptrs, right_ptrs);

        bool removed_bd = bd_it->l_size() == 0 || bd_it->r_size() == 0;
        if (removed_bd) {
            bd_it->active = false;
            bd_it->remove();
        }

        auto filter_result = filter_domains(workspace,
                bdll, left_ptrs, right_ptrs, g0.adj_lists, g1.adj_lists, v, w, matching_size_goal, bound);

        if (!filter_result.quit_early) {
            current.push_back(VtxPair(v, w));

            solve(workspace, g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, solution_count,
                    g0_remaining_deg, matching_size_goal, filter_result.new_bound);

            current.pop_back();
            unfilter_domains(workspace, bdll, left_ptrs, right_ptrs,
                    filter_result.split_bds_list, filter_result.deleted_bds_list);
        }

        if (removed_bd) {
            bd_it->active = true;
            bd_it->reinsert();
        }
        unassign(v, w, bd_it, left_ptrs, right_ptrs);

        if (incumbent.size()==matching_size_goal)
            break;
    }

    {
        if (bd_it->l_size() <= bd_it->r_size())
            --bound;

        if (bound < matching_size_goal)
            return;

        delete_v(v, left_ptrs);

        bool removed_bd = bd_it->l_size() == 0;
        if (removed_bd) {
            bd_it->active = false;
            bd_it->remove();
        }

        solve(workspace, g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, solution_count,
                g0_remaining_deg, matching_size_goal, bound);

        if (removed_bd) {
            bd_it->active = true;
            bd_it->reinsert();
        }
        undelete_v(v, bd_it, left_ptrs);
    }

}

vector<int> get_RI_style_static_order(const SparseGraph & g0, double g1_density)
{
    bool g1_is_dense = g1_density > 0.5;

    vector<int> order;
    order.reserve(g0.n);

    // Set 0 is vertices that are in the order
    // Set 1 is vertices adjacent to those in set 0
    // Set 2 is other vertices
    vector<int> vtx_to_set(g0.n, 2);

    vector<vector<int>> scores(g0.n, {0,0,0});
    for (int i=0; i<g0.n; i++) {
        scores[i][2] = g0.adj_lists[i].size();
    }
    for (int i=0; i<g0.n; i++) {
        int v = -1;
        if (g1_is_dense) {
            vector<int> best_score = {INT_MAX, INT_MAX, INT_MAX};
            for (int j=0; j<g0.n; j++) {
                if (vtx_to_set[j] != 0 && scores[j] < best_score) {
                    best_score = scores[j];
                    v = j;
                }
            }
        } else {
            vector<int> best_score = {-1,-1,-1};
            for (int j=0; j<g0.n; j++) {
                if (vtx_to_set[j] != 0 && scores[j] > best_score) {
                    best_score = scores[j];
                    v = j;
                }
            }
        }
        order.push_back(v);
//        cout << order.back() << " ";
        int prev_set_of_v = vtx_to_set[v];
        vtx_to_set[v] = 0;
        for (int w : g0.adj_lists[v]) {
            --scores[w][prev_set_of_v];
            ++scores[w][0];
            if (vtx_to_set[w] == 2) {
                vtx_to_set[w] = 1;
                for (int u : g0.adj_lists[w]) {
                    --scores[u][2];
                    ++scores[u][1];
                }
            }
        }
    }
//    cout << endl;
//    for (int i=0; i<g0.n; i++) {
//        for (int x : scores[i]) cout << x << " ";
//        cout << endl;
//    }
    return order;
}

// Returns a common subgraph and the number of induced subgraph isomorphisms found
// vv0 and vv1 are vertex orders
vector<VtxPair> mcs(SparseGraph & g0, SparseGraph & g1, double g1_density,
        vector<int> & vv0, vector<int> & vv1)
{
    //for (int i=0; i<g0.n; i++) {
    //    cout << i << "   ";
    //    for (int v : g0.adj_lists[i]) {
    //        cout << " " << v;
    //    }
    //    cout << endl;
    //}
    //cout << endl;
    //for (int i=0; i<g0.n; i++) {
    //    cout << i << "   ";
    //    for (int v : g0.in_edge_lists[i]) {
    //        cout << " " << v;
    //    }
    //    cout << endl;
    //}
    //cout << endl;
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions
    left.reserve(g0.n);
    right.reserve(g1.n);

    vector<Ptrs> left_ptrs(g0.n);
    vector<Ptrs> right_ptrs(g1.n);

    BDLL bdll;

    vector<int> order;
    if (arguments.heuristic == heur_RI) {
        order = get_RI_style_static_order(g0, g1_density);
    }
    Workspace workspace {vv0, vv1, order};

    std::vector<unsigned int> left_labels(g0.label.begin(), g0.label.end());
    // sort and deduplicate
    std::sort( left_labels.begin(), left_labels.end() );
    left_labels.erase( std::unique( left_labels.begin(), left_labels.end() ), left_labels.end() );

    unsigned bound = 0;

    // Create a bidomain for each label that appears in the pattern graph
    for (unsigned int label : left_labels) {
        int left_len = 0;
        int right_len = 0;

        NewBidomain *new_elem = workspace.get_from_free_list();
        new_elem->insert_before(&bdll.head);

        for (int i=0; i<g0.n; i++) {
            if (g0.label[i]==label) {
                left.push_back(i);
                left_ptrs[i].vtx_it = std::prev(left.end());
                left_ptrs[i].bd_it = new_elem;
                ++left_len;
            }
        }
        for (int i=0; i<g1.n; i++) {
            if (g1.label[i]==label) {
                right.push_back(i);
                right_ptrs[i].vtx_it = std::prev(right.end());
                //cout << *right_ptrs[i].vtx_it << endl;
                right_ptrs[i].bd_it = new_elem;
                ++right_len;
            }
        }

        bound += std::min(left_len, right_len);

        new_elem->initialise(left.end() - left_len,
                             right.end() - right_len,
                             left.end(),
                             right.end());
    }

    vector<VtxPair> incumbent;
    for (int matching_size_goal = g1.n; matching_size_goal >= 0; matching_size_goal--) {
        if (int(incumbent.size()) == matching_size_goal)
            return incumbent;
        vector<VtxPair> current;
        long long solution_count = 0;
        vector<int> g0_remaining_deg {};  // used for heur_C
        for (int i=0; i<g0.n; i++)
            g0_remaining_deg.push_back(g0.adj_lists[i].size());
        solve(workspace, g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, solution_count,
                g0_remaining_deg, matching_size_goal, bound);

        if (int(incumbent.size()) == matching_size_goal)
            return incumbent;
    }
    return {};
}

vector<int> calculate_degrees(const SparseGraph & g) {
    vector<int> degree;
    degree.reserve(g.n);
    for (int v=0; v<g.n; v++) {
        degree.push_back(g.adj_lists[v].size());
    }
    return degree;
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

double calc_density(vector<int> & deg)
{
    int n = deg.size();
    long deg_total = 0;
    for (auto d : deg) {
        deg_total += d;
    }
    return (double) deg_total / (n * (n-1));
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : arguments.gfd ? 'G' :
        arguments.vf ? 'V' : 'B';
    struct SparseGraph g0 = readGraph(arguments.filename1, format, false,
            arguments.edge_labelled, arguments.vertex_labelled);
    struct SparseGraph g1 = readGraph(arguments.filename2, format, false,
            arguments.edge_labelled, arguments.vertex_labelled);

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }

    auto start = std::chrono::steady_clock::now();

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);
    double g0_density = calc_density(g0_deg);
    double g1_density = calc_density(g1_deg);

    bool reverse_deg_heur = false;
    if (arguments.deg_heur == 0 && g1_density > 0.5) reverse_deg_heur = true;
    if (arguments.deg_heur == 1 && g1_density > g0_density) reverse_deg_heur = true;

    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);

    if (reverse_deg_heur) {
        std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) { return g0_deg[a] < g0_deg[b]; });
        std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) { return g1_deg[a] < g1_deg[b]; });
    } else {
        std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) { return g0_deg[a] > g0_deg[b]; });
        std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) { return g1_deg[a] > g1_deg[b]; });
    }

    auto solution = mcs(g0, g1, g1_density, vv0, vv1);

    auto stop = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }


    cout << "Nodes:                      " << nodes << endl;
    cout << "Time (ms):                  " << time_elapsed << endl;
    if (aborted) {
        cout << "TIMEOUT" << endl;
    } else {
        if (!check_sol(g0, g1, solution))
            fail("*** Error: Invalid solution\n");

        cout << "Solution size " << solution.size() << std::endl;
        for (int i=0; i<g0.n; i++)
            for (unsigned int j=0; j<solution.size(); j++)
                if (solution[j].v == i)
                    cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
        cout << std::endl;
    }
}
