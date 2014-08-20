#pragma once

#include <iostream>
#include "AllelesNode.h"


class AllelesTree {
public:
    AllelesTree();
    ~AllelesTree();

    AllelesTree(const string& newick_str);
    AllelesNode* get_root();
    ProbVec get_allele_node_probabilities(AllelesNode*);
    string convert_tree_to_str(AllelesNode*);
    void destroy_tree();

private:
    AllelesNode* convert_str_to_tree(const string& newick_str);
    AllelesNode* build_tree(const string& newick_str, unsigned from, unsigned to, AllelesNode* parent);
    ProbVec compute_prob_vec(AllelesNode*);
    void destroy_tree(AllelesNode*);

    void vec_sum(const ProbVec&, ProbVec&);
    void vec_on_number_multiply(const ProbVec&, double num, ProbVec&);
    ProbVec get_crossing_prob_vec(const ProbVec&, const ProbVec&); 

private:
    AllelesNode* root;
    string newick_str;
};

