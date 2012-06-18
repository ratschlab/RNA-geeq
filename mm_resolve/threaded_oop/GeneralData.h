#ifndef __GENERALDATA_H__
#define __GENERALDATA_H__

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

using namespace std;

struct GeneralData {

    // general coverage structures
    map <string, int> chr_num;
    vector <unsigned int> chr_size;
    map <int, vector<unsigned short> > coverage_map;

    // best hit maps (will stay empty for batch setting)
    unordered_map <string, size_t, hash<string> > best_left;
    unordered_map <string, size_t, hash<string> > best_right;

};
#endif
