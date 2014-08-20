#pragma once

#include <vector>
#include <map>
#include <string>

using std::string;
using std::vector;
using std::map;

typedef vector<double> ProbVec;

const int ALLELES_NUM = 3;

struct AllelesNode
{
    AllelesNode(AllelesNode *left, AllelesNode *right);   
    void init_probabilistic_vec(const string& allele_str);

private:
    static map<string, int> init_allele_ord() {
        map<string, int> m;
        m["AA"] = 0;
        m["Aa"] = 1;
        m["aA"] = 1;
        m["aa"] = 2;
        return m;
    }

    static vector< vector< ProbVec > > init_prob_matrix() {
        double const matrix[ALLELES_NUM][ALLELES_NUM][ALLELES_NUM] = { 
            { {1, 0, 0},     {0.5, 0.5, 0},     {0, 1, 0}     }, 
            { {0.5, 0.5, 0}, {0.25, 0.5, 0.25}, {0, 0.5, 0.5} }, 
            { {0, 1, 0},     {0, 0.5, 0.5},     {0, 0, 1}     }
        };

        vector< vector< ProbVec > > tmp_matrix(ALLELES_NUM, vector< ProbVec >(ALLELES_NUM, ProbVec(ALLELES_NUM, 0)));
        for (int i = 0; i < ALLELES_NUM; ++i) {
            for (int j = 0; j < ALLELES_NUM; ++j) {
                for (int k = 0; k < ALLELES_NUM; ++k) {
                    tmp_matrix[i][j][k] = matrix[i][j][k];
                }
            }
        }

        return tmp_matrix;
    }

 public:
    // alleles probabilistic vector (x, y, z)
    // AA = (1, 0, 0)
    // Aa = (0, 1, 0)
    // aa = (0, 0, 1)
    // Aa x Aa = (0.25, 0.5, 0.25), because Aa x Aa = {AA, Aa, Aa, aa}
    ProbVec prob_vec;
    bool is_computed_prob_vec;
    string allele_str;

    AllelesNode *left;
    AllelesNode *right;

    static const map<string, int> allele_ord;
    static const vector< vector< ProbVec > > prob_matrix;
};

