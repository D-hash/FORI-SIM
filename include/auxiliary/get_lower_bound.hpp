#pragma once
#include "graph.hpp"

template<typename T, typename U>
unsigned int compute_lower_bound(const graph<T,U>& g1, const graph<T,U>& g2,
                                 unsigned int upper_bound) {
    using node = uint32_t;

    unsigned int n1 = g1.number_of_nodes();
    unsigned int n2 = g2.number_of_nodes();
    unsigned int m1 = g1.number_of_edges();
    unsigned int m2 = g2.number_of_edges();

    // Step 1: Size-based lower bound
    unsigned int lb = std::abs((int)n1 - (int)n2) + std::abs((int)m1 - (int)m2) / 2;
    if (lb > upper_bound) return lb;

    // Step 2: Node label mismatch
    lb = std::max(n1, n2);
    std::unordered_map<T, int> vlabel_cnt;
    for (node i = 0; i < n1; ++i) ++vlabel_cnt[g1.get_node_label(i)];
    for (node i = 0; i < n2; ++i) {
        T vl = g2.get_node_label(i);
        if (vlabel_cnt[vl] > 0) {
            --vlabel_cnt[vl];
            --lb;
        }
    }
    if (lb > upper_bound) return lb;

    // Step 3: Degree mismatch
    std::vector<int> degree_q(n1), degree_g(n2);
    for (node i = 0; i < n1; ++i) degree_q[i] = g1.get_degree(i);
    for (node i = 0; i < n2; ++i) degree_g[i] = g2.get_degree(i);

    std::vector<int> degrees_cnt_q(n1 + 1, 0), degrees_cnt_g(n2 + 1, 0);
    int max_degree_q = 0, max_degree_g = 0;
    for (int d : degree_q) {
        ++degrees_cnt_q[d];
        if (d > max_degree_q) max_degree_q = d;
    }
    for (int d : degree_g) {
        ++degrees_cnt_g[d];
        if (d > max_degree_g) max_degree_g = d;
    }

    unsigned int de = 0, ie = 0;
    while (max_degree_q > 0 && max_degree_g > 0) {
        if (degrees_cnt_q[max_degree_q] == 0) { --max_degree_q; continue; }
        if (degrees_cnt_g[max_degree_g] == 0) { --max_degree_g; continue; }

        unsigned int td = std::min(degrees_cnt_q[max_degree_q], degrees_cnt_g[max_degree_g]);
        if (max_degree_q > max_degree_g) de += td * (max_degree_q - max_degree_g);
        else ie += td * (max_degree_g - max_degree_q);

        degrees_cnt_q[max_degree_q] -= td;
        degrees_cnt_g[max_degree_g] -= td;
    }
    while (max_degree_q > 0) {
        de += max_degree_q * degrees_cnt_q[max_degree_q];
        --max_degree_q;
    }
    while (max_degree_g > 0) {
        ie += max_degree_g * degrees_cnt_g[max_degree_g];
        --max_degree_g;
    }

    de = (de + 1) / 2;
    ie = (ie + 1) / 2;

    unsigned int edge_lb = de + ie;
    if (lb + edge_lb > upper_bound) return lb + edge_lb;

    // Step 4: Edge label matching
    std::unordered_map<U, int> elabel_cnt;
    for (auto e : g1.get_edgelist()) ++elabel_cnt[g1.get_edge_label(e)];

    unsigned int common_elabel_cnt = 0;
    for (auto e : g2.get_edgelist()) {
        U el = g2.get_edge_label(e);
        if (elabel_cnt[el] > 0) {
            --elabel_cnt[el];
            ++common_elabel_cnt;
        }
    }
    common_elabel_cnt /= 2;

    unsigned int e_cnt = std::max(m1, m2) / 2;
    unsigned int edge_label_lb = std::max({
        de + m2 / 2 - common_elabel_cnt,
        ie + m1 / 2 - common_elabel_cnt,
        e_cnt - common_elabel_cnt
    });

    edge_lb = std::max(edge_lb, edge_label_lb);

    return lb + edge_lb;
}


