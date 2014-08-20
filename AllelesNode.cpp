#include "AllelesNode.h"

AllelesNode::AllelesNode(AllelesNode *left, AllelesNode *right) {
    this->left = left; 
    this->right = right;
    is_computed_prob_vec = false;
    prob_vec.resize(ALLELES_NUM, 0);
}

void AllelesNode::init_probabilistic_vec(const string& allele_str) {
    int ord_num = allele_ord.at(allele_str);
    this->allele_str = allele_str;
    prob_vec[ord_num] = 1;
    is_computed_prob_vec = true;
}

const map<string, int> AllelesNode::allele_ord = init_allele_ord();
const vector< vector< ProbVec > > AllelesNode::prob_matrix = init_prob_matrix();

