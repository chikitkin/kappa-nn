/* 
 * File:   cvs2YamlParticles.cpp
 * Author: aspera
 *
 * Created on 4 марта 2016 г., 15:26
 */

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <memory>


using namespace std;


#include "yaml-cpp/yaml.h"

int getdir(string dir, vector<string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(dir.c_str())) == nullptr) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != nullptr) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

string getFileExt(const string& s) {
    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }
    return ("");
}

string getFileName(const string& s) {
    size_t i = s.rfind('.', s.length());
    if (i != 0) {
        return (s.substr(0, i));
    }
    return ("");
}

void splitIntoTokens(std::string &line, std::vector<std::string> &result, char delimiter = ',') {
    std::stringstream lineStream(line);
    std::string cell;

    while (std::getline(lineStream, cell, delimiter)) {
        result.push_back(cell);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr << " Use: " << argv[0] << " filename.csv" << endl;
        return -1;
    }
    string csvFile(argv[1]);
    // ---------------------
    cout << "Convert file:" << csvFile << " to yaml" << endl;
    ifstream in(csvFile.c_str());
    if (!in.is_open()) {
        cerr << "Error: file:" << csvFile << " can not open!!!" << endl;
        return -1;
    }

    // ------------
    string line;
    vector<string> headers;
    getline(in, line);
    splitIntoTokens(line, headers, '\t');
    if (headers[0].compare("Particle")) {
        cerr << "First row name must be Particle -- partice name  :" << headers[0] << endl;
        return -1;
    }
    // --------
    int numRow = headers.size();
    vector<vector<double>> data(numRow - 1);
    vector<string> name;

    while (getline(in, line)) {
        vector< string > vec;
        splitIntoTokens(line, vec, '\t');
        if (vec.size() > numRow) {
            cerr << "Num row in line:" << line << " must be:" << numRow << " or <"<< endl;
            return -1;
        }
        if (vec[0] != "")
            name.push_back(vec[0]);
        for (int i = 1; i < vec.size(); i++) {
            if (vec[i] != "") {
                data[i - 1].push_back(stod(vec[i]));
            }
        }
    }
    // --------
    //    cout << "NumName:" << name.size() << endl;
    //    for (int i = 0; i < name.size(); i++) {
    //        cout << "name[" << i + 1 << "]:" << name[i] << endl;
    //    }
    //        for (int i = 0; i < data.size(); i++) {
    //            cout << "RowInd:" << i + 1 << " Name:" << headers[i + 1] << " numLine:" << data[i].size() << endl;
    //        }
    // ------------------   
    YAML::Node particle;
    for (int i = 0; i < data.size(); i++) {
        string label = headers[i + 1];
        //        cout << "RowInd:" << i + 1 << " Name:" << label << " numLine:" << data[i].size() << endl;
        if (data[i].size() == 1) {
            particle[label.c_str()] = data[i][0];
        } else if (data[i].size() > 1) {
            particle[label.c_str()] = data[i];
        }
    }
    //---
    YAML::Node node;
    node[name[0].c_str()] = particle;

    // ---
    YAML::Emitter out;
    out.SetSeqFormat(YAML::Flow);
    //    out.SetMapFormat(YAML::Flow);
    out << node;
    // ---
    string josnFile = getFileName(csvFile) + ".yaml";
    cout << "Output to file:" << josnFile << endl;
    std::ofstream fout(josnFile.c_str());
    fout << out.c_str() << endl;
    return 0;
}



