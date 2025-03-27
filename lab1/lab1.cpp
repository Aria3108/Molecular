#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

void saveAlignmentToFile(const pair<string, string>& alignment, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: No se pudo abrir el archivo para escribir." << endl;
        return;
    }
    file << alignment.first << "\n";
    file << alignment.second << "\n";
    file.close();
    cout << "Alineacion guardada en " << filename << endl;
}

pair<string, string> backtrackUnique(const vector<vector<int>>& score, const string& seq1, const string& seq2) {
    string align1 = "", align2 = "";
    int i = seq1.length(), j = seq2.length();
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH)) {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--;
            j--;
        } else if (i > 0 && score[i][j] == score[i - 1][j] + GAP) {
            align1 = seq1[i - 1] + align1;
            align2 = "-" + align2;
            i--;
        } else {
            align1 = "-" + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        }
    }
    return {align1, align2};
}

void needlemanWunsch(const string& seq1, const string& seq2) {
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

    auto uniqueAlignment = backtrackUnique(score, seq1, seq2);
    cout << uniqueAlignment.first << endl;
    cout << uniqueAlignment.second << endl;

    saveAlignmentToFile(uniqueAlignment, "alineacion.txt");
}

int main() {
    string inputA = { "GCGGCCCCGGGGGGGGCGG" };
    string inputB = { "CGCCGGGGGCCCGGGC" };

    needlemanWunsch(inputA, inputB);
    return 0;
}
