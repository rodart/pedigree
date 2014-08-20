#include "AllelesTree.h"


AllelesTree::AllelesTree(const string& newick_str) {
    this->newick_str = newick_str.substr(0, newick_str.size() - 1);
}

AllelesTree::~AllelesTree() {
    destroy_tree();
}

AllelesNode* AllelesTree::convert_str_to_tree(const string& newick_str) {
    return build_tree(newick_str, 0, newick_str.size(), new AllelesNode(NULL, NULL));
}

string AllelesTree::convert_tree_to_str(AllelesNode *node) {
    if (node->left && node->right) {
        return "(" + convert_tree_to_str(node->left) + "," + convert_tree_to_str(node->right) + ")";
    } 
    return node->allele_str;
}

AllelesNode* AllelesTree::get_root() {
    root = convert_str_to_tree(newick_str);
    return root;
}

ProbVec AllelesTree::get_allele_node_probabilities(AllelesNode *node) {
    if (!node->is_computed_prob_vec) {
        compute_prob_vec(node);
    }
    return node->prob_vec;
}       

void AllelesTree::destroy_tree() {
    destroy_tree(root);
}

 
AllelesNode* AllelesTree::build_tree(const string& newick_str, unsigned from, unsigned to, AllelesNode* parent) {
    if (newick_str[from] == 'a' || newick_str[from] == 'A') {
        parent->init_probabilistic_vec(newick_str.substr(from, 2));
        return parent;
    }

    int open_bracket_num = 0;
    int cur_pos = from;

    for (vector<int>::size_type i = from; i < to; i++) {
        char cur_symbol = newick_str[i];

        if (cur_symbol == '(') {
            open_bracket_num++;
        }
        if (cur_symbol == ')') {
            open_bracket_num--;
        }

        if (open_bracket_num == 1 && cur_symbol == ',') {
            parent->left = build_tree(newick_str, cur_pos + 1, i, new AllelesNode(NULL, NULL));
            cur_pos = i;
        } 
        if (open_bracket_num == 0) {
            parent->right = build_tree(newick_str, cur_pos + 1, i, new AllelesNode(NULL, NULL));
            cur_pos = i;
        }
    }

    return parent;
}
    

ProbVec AllelesTree::compute_prob_vec(AllelesNode *node) {
    if (!node->left && !node->right) {
        return node->prob_vec;
    }
    if (node->left && !node->right) {
        return compute_prob_vec(node->left);
    }
    if (node->right && !node->left) {
        return compute_prob_vec(node->right);
    }

    node->is_computed_prob_vec = true;
    node->prob_vec = get_crossing_prob_vec(compute_prob_vec(node->left), compute_prob_vec(node->right));
    return node->prob_vec;
}

void AllelesTree::vec_sum(const ProbVec& a, ProbVec& b) {
    for (vector<int>::size_type i = 0; i < std::min(a.size(), b.size()); ++i) {
        b[i] += a[i];
    }
    return;
}

void AllelesTree::vec_on_number_multiply(const ProbVec& vec, double num, ProbVec& result) {
    for (vector<int>::size_type i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * num;
    }
    return;
}

ProbVec AllelesTree::get_crossing_prob_vec(const ProbVec& a, const ProbVec& b) {
    vector< ProbVec > crossing_prob_matrix(ALLELES_NUM, ProbVec(ALLELES_NUM, 0));
    for (vector<int>::size_type i = 0; i < a.size(); ++i) {
        for (vector<int>::size_type j = 0; j < b.size(); ++j) {
            crossing_prob_matrix[i][j] = a[i] * b[j];
        }
    }

    ProbVec res_prob_vec(ALLELES_NUM, 0);
    for (int i = 0; i < ALLELES_NUM; ++i) {
        for (int j = 0; j < ALLELES_NUM; ++j) {
            ProbVec tmp_res(ALLELES_NUM, 0);
            vec_on_number_multiply(AllelesNode::prob_matrix[i][j], crossing_prob_matrix[i][j], tmp_res);
            vec_sum(tmp_res, res_prob_vec);
        }
    }
    
    return res_prob_vec;
}

void AllelesTree::destroy_tree(AllelesNode *node) {
    if (node) {
        destroy_tree(node->left);
        destroy_tree(node->right);
        delete node;
    }
}
