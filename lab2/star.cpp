#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <limits>
#include <utility>

using namespace std;

typedef pair<string, string> alin;

const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

int score(alin alignment){

    int score = 0;
    for(int i = 0; i < alignment.first.size(); i++ ){
        if(alignment.first[i] == '-' || alignment.second[i] == '-')
            score -= 2;
        else if(alignment.first[i] == alignment.second[i])
            score++;
        else 
            score--;
    }
    return score;
}

int contRup(const string& a, const string& b) {
    if (a.size() != b.size()) return -1; 

    int rupturas = 0;
    char estado_anterior;

    estado_anterior = (a[0] == '-' || b[0] == '-') ? 'G' : 'A';

    for (size_t i = 1; i < a.size(); ++i) {
        char estado_actual = (a[i] == '-' || b[i] == '-') ? 'G' : 'A';

        if (estado_actual != estado_anterior) {
            rupturas++;
        }
        estado_anterior = estado_actual;
    }
    return rupturas;
}

void traceback(const vector<vector<int>>& score, const string& seq1, const string& seq2, 
    int i, int j, string align1, string align2, set<alin>& results) {
    if (i == 0 && j == 0) {
        results.insert({align1, align2});
        return;
    }

    if (i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH)) {
        traceback(score, seq1, seq2, i - 1, j - 1, seq1[i - 1] + align1, seq2[j - 1] + align2, results);
    }
    if (i > 0 && score[i][j] == score[i - 1][j] + GAP) {
        traceback(score, seq1, seq2, i - 1, j, seq1[i - 1] + align1, "-" + align2, results);
    }
    if (j > 0 && score[i][j] == score[i][j - 1] + GAP) {
        traceback(score, seq1, seq2, i, j - 1, "-" + align1, seq2[j - 1] + align2, results);
    }
}

set<alin> needlemanWunsch(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));

    for (int i = 0; i <= m; i++) {
        score[i][0] = i * GAP;
    }
    for (int j = 0; j <= n; j++) {
        score[0][j] = j * GAP;
    }

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match = score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH);
            int deleteGap = score[i - 1][j] + GAP;
            int insertGap = score[i][j - 1] + GAP;
            score[i][j] = max({match, deleteGap, insertGap});
        }
    }

    set<alin> results;
    traceback(score, seq1, seq2, m, n, "", "", results);
    return results;
}

alin bestAlignment(const string& seq1, const string& seq2,vector<vector<int>>& scores, int i, int j) {

    set<alin>alignments = needlemanWunsch(seq1, seq2);

    alin bestAlign;
    int minRupturas = numeric_limits<int>::max(); 

    for (const auto& alignment : alignments) {
        int rupturas = contRup(get<0>(alignment), get<1>(alignment));
        
        if (rupturas < minRupturas) {
            bestAlign = alignment;
            minRupturas = rupturas;
        }
    }

    scores[i][j] = i == j ? 0: score(bestAlign);

    return bestAlign;
}

int getBest(vector<vector<int>> scores, int &bestScore) {

    int max = numeric_limits<int>::min(); 
    int best;
    for(int i = 0; i < scores.size(); i++){
        int score = accumulate(scores[i].begin(), scores[i].end(), 0);

        if (score > max) {
            max = score;
            best = i;
        } 
    }
    bestScore = max;
    return best;
}

void printTable(vector<vector<int>> table) {
    for (int i = 0; i < table.size(); i++) {
        for (int j = 0; j < table[i].size(); j++) {
            cout << table[i][j] << "\t";
        }
        cout << endl;
    }
}

void printMSA(vector<vector<alin>> alignments, int best){

    for(int i = 0; i < alignments[best].size(); i++){
        cout << "S" <<i << ": " << alignments[best][i].second << endl;
    }
}

void MSA(vector<string> cadenas){
    vector<vector<int>> scores(cadenas.size(), vector<int>(cadenas.size()));
    vector<vector<alin>> alignments(cadenas.size(), vector<alin>(cadenas.size()));

    for(int i = 0; i < cadenas.size(); i++){ 
        for(int j = 0; j < cadenas.size(); j++){
            alignments[i][j] = bestAlignment(cadenas[i], cadenas[j], scores, i, j);
        }    
    }

    printTable(scores);
    cout << endl;
    int bestScore;
    int best = getBest(scores, bestScore);
    cout << "Cadena centro: " << best << endl;
    cout << "Score: " << bestScore << endl << endl;

    printMSA(alignments, best);

    

}

int main() {
    vector<string> cadenas = 
    { "ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT", "ACTGACC"};

    vector<string> F = 
    {"TGCCGGCAGGGATGTGCTTG", "GTTTAGGTTTTTGCTTATGCAGCATCCA", "GGAAAAGCACAGAACTGGCCAACA",
        "GCCAGTTGGTTGATTTCCACCTCCA", "ACCCCCGACATGCAGAAGCTG", "TGACGTGTCTGCTCCACTTCCA"};

        vector<string> R = 
    {"TGCTTGCAGTTTGCTTTCACTGATGGA", "TCAGGTACCCTGACCTTCTCTGAAC", "GTGGGTTGTAAAGGTCCCAAATGGT",
        "TGCCTTGGGTCCCTCTGACTGG", "GTGGTGCATTGATGGAAGGAAGCA", "AGTGAGAGGAGCTCCCAGGGC"};

    MSA(F);


    return 0;
}

