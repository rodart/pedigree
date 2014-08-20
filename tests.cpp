#include <iostream>

#include "AllelesTree.h"
#include "AllelesNode.h"

const double EPS = 0.001;

bool elem_sum_is_one(const ProbVec& res_prob_vec) {
    double sum = 0;
    for (int i = 0; i < res_prob_vec.size(); ++i) {
        sum += res_prob_vec[i];
    }
    return (std::abs(sum - 1) < EPS);
}

bool is_right_probabilities(const ProbVec& res_prob_vec, const ProbVec& right_prob_vec) {
    int num = 0;
    for (int i = 0; i < res_prob_vec.size(); ++i) {
        if (std::abs(res_prob_vec[i] - right_prob_vec[i]) < EPS) {
            num++;
        }
    }
    return num == ALLELES_NUM;
}

void base_test(const string& newick_str, const vector<double>& right_prob_vec) {
    std::cout << "Current test is " << newick_str << std::endl;

    AllelesTree alleles_tree(newick_str);
    AllelesNode* root = alleles_tree.get_root();
    
    string tree_str = alleles_tree.convert_tree_to_str(root) + ";";
    if (newick_str.compare(tree_str) == 0) {
        std::cout << "OK: Correctly converting newick string to tree" << std::endl;
    } else {
        std::cout << "ERROR: Incorrectly converted newick string to tree" << std::endl;
    }

    ProbVec res_prob_vec = alleles_tree.get_allele_node_probabilities(root);
    
    if (elem_sum_is_one(res_prob_vec)) {
        std::cout << "OK: Sum of prob vector = 1" << std::endl;
    } else {
        std::cout << "ERROR: Sum of prob vector != 1" << std::endl;
    }

    if (right_prob_vec.size()) {
        if (is_right_probabilities(res_prob_vec, right_prob_vec)) {
            std::cout << "OK: All result probabilities are correct" << std::endl;
        } else {
            std::cout << "ERROR: Result probabilities are not correct" << std::endl;
        }
    }

    std::cout << std::endl;
}

string make_long_test_string(int depth) {
    if (depth > 0) {
        return "(" + make_long_test_string(depth - 1 - rand() % 3) + "," + make_long_test_string(depth - 1 - rand() % 3) + ")";
    }

    int rnd = rand() % 3;
    if (rnd == 0) return "AA";
    if (rnd == 1) return "Aa";
    return "aa";
}

int main() {
    double arr[] = {0.156, 0.5, 0.344};
    vector<double> prob_vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    base_test("((((Aa,aa),(Aa,Aa)),((aa,aa),(aa,AA))),Aa);", prob_vec);

    double arr2[] = {0, 0.875, 0.125};
    vector<double> prob_vec2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
    base_test("(aa,((Aa,AA),AA));", prob_vec2);

    string test_str = make_long_test_string(15) + ";";
    base_test(test_str, vector<double>());

    return 0;
}