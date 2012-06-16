
#include <assert.h>

#include "OnlineData.h"
#include "Utils.h"
#include "config.h"

extern Config* conf;

extern pthread_mutex_t mutex_coverage;
extern pthread_mutex_t mutex_fifo;
extern pthread_mutex_t mutex_done;
extern pthread_mutex_t mutex_best_left;
extern pthread_mutex_t mutex_best_right;

void OnlineData::process_data_online(GeneralData* genData) {

    vector<vector<Alignment>::iterator> active_left_reads;
    vector<vector<Alignment>::iterator> active_right_reads;

    vector<vector<Alignment>::iterator>::iterator lv_idx;
    vector<vector<Alignment>::iterator>::iterator rv_idx;
    vector<vector<Alignment>::iterator>::iterator curr_best;

    vector<Alignment>::iterator best_left_idx;
    vector<Alignment>::iterator best_right_idx;

    set<vector<Alignment>::iterator> ignore_idx_left;
    set<vector<Alignment>::iterator> ignore_idx_right;

    unsigned int num_changed = 0;
    bool found_pairs = false;

    if (conf->pre_filter) {
        ignore_idx_left = filter_alignments(this->left_reads);
        pthread_mutex_lock(&mutex_best_left);
        if (ignore_idx_left.find(this->left_reads.begin() + genData->best_left[this->last_id]) != ignore_idx_left.end()) {
            (this->left_reads.begin() + genData->best_left[this->last_id])->is_best = false;
            for (vector<Alignment>::iterator v_idx = this->left_reads.begin(); v_idx != this->left_reads.end(); v_idx++) {
                if (ignore_idx_left.find(v_idx) == ignore_idx_left.end()) {
                    genData->best_left[this->last_id] = v_idx - this->left_reads.begin();
                    v_idx->is_best = true;
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_left);

        ignore_idx_right = filter_alignments(this->right_reads);
        pthread_mutex_lock(&mutex_best_right);
        if (ignore_idx_right.find(this->right_reads.begin() + genData->best_right[this->last_id]) != ignore_idx_right.end()) {
            (this->right_reads.begin() + genData->best_right[this->last_id])->is_best = false;
            for (vector<Alignment>::iterator v_idx = this->right_reads.begin(); v_idx != this->right_reads.end(); v_idx++) {
                if (ignore_idx_right.find(v_idx) == ignore_idx_right.end()) {
                    genData->best_right[this->last_id] = v_idx - this->right_reads.begin();
                    v_idx->is_best = true;
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_right);
    }

    get_active_reads(this->last_id, ignore_idx_left, ignore_idx_right, active_left_reads, active_right_reads, genData, found_pairs);

    if (found_pairs) {
        bool best_found = false;
        for(lv_idx = active_left_reads.begin(), rv_idx = active_right_reads.begin(); lv_idx < active_left_reads.end() && rv_idx  < active_right_reads.end(); lv_idx++, rv_idx++) {
            if ((*rv_idx)->is_best && (*lv_idx)->is_best) {
                best_left_idx = (*lv_idx);
                best_right_idx = (*rv_idx);
                best_found = true;
                break;
            }
        }
        if (! best_found) 
            assert(best_found);
        bool changed = false;
        for(lv_idx = active_left_reads.begin(), rv_idx = active_right_reads.begin(); lv_idx < active_left_reads.end() && rv_idx  < active_right_reads.end(); lv_idx++, rv_idx++) {
            if ((*lv_idx)  == best_left_idx && (*rv_idx) == best_right_idx)
                continue;
            if (compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx)) {
                changed = true;
                best_left_idx->is_best = false;
                best_right_idx->is_best = false;
                best_left_idx->update_coverage_map(0);
                best_right_idx->update_coverage_map(0);
                best_left_idx = *lv_idx;
                best_right_idx = *rv_idx;
                best_left_idx->is_best = true;
                best_right_idx->is_best = true;
                best_left_idx->update_coverage_map(1);
                best_right_idx->update_coverage_map(1);
                
                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[this->last_id] = (best_left_idx - this->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);
                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[this->last_id] = (best_right_idx - this->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
        }
        if (changed) 
            num_changed++;
    } else {
        bool best_found = false;
        for (lv_idx = active_left_reads.begin(); lv_idx != active_left_reads.end(); lv_idx++) {
            if ((*lv_idx)->is_best) {
                curr_best = lv_idx;
                best_found = true;
                break;
            }
        }
        if (active_left_reads.size() > 0)
            assert(best_found);
        bool changed = false;
        for (lv_idx = active_left_reads.begin(); lv_idx != active_left_reads.end(); lv_idx++) {
            if (lv_idx == curr_best)
                continue;

            // check if lv_idx < curr_best
            if (compare_single(*lv_idx, *curr_best)) {
                changed = true;
                (*curr_best)->is_best = false;
                (*curr_best)->update_coverage_map(0);
                curr_best = lv_idx;          
                (*curr_best)->is_best = true;
                (*curr_best)->update_coverage_map(1);
                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[this->last_id] = (*curr_best - this->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);
            }
        }
        if (changed) 
            num_changed++;

        best_found = false;
        for (rv_idx = active_right_reads.begin(); rv_idx != active_right_reads.end(); rv_idx++) {
            if ((*rv_idx)->is_best) {
                curr_best = rv_idx;
                best_found = true;
                break;
            }
        }
        if (active_right_reads.size() > 0)
            assert(best_found);
        changed = false;
        for (rv_idx = active_right_reads.begin(); rv_idx != active_right_reads.end(); rv_idx++) {
            if (rv_idx == curr_best)
                continue;

            // check if rv_idx < curr_best
            if (compare_single(*rv_idx, *curr_best)) {
                changed = true;
                (*curr_best)->is_best = false;
                (*curr_best)->update_coverage_map(0);
                curr_best = rv_idx;          
                (*curr_best)->is_best = true;
                (*curr_best)->update_coverage_map(1);
                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[this->last_id] = (*curr_best - this->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
        }
        if (changed) 
            num_changed++;
    }
    delete this;
}

void OnlineData::get_active_reads(string read_id, set<vector<Alignment>::iterator> &ignore_reads_left, set<vector<Alignment>::iterator> &ignore_reads_right, vector<vector<Alignment>::iterator> &active_left_reads, vector<vector<Alignment>::iterator> &active_right_reads, GeneralData* genData, bool &found_pairs) {

    double insert1 = 0.0;
    double insert2 = 0.0;

    // if pair processing, identify all compatible pairs
    // pairs are compatible, if they show same chr, opposite strands and
    // have an inner distance within the insert size range

    active_left_reads.clear();
    active_right_reads.clear();

    if (conf->use_pair_info && ! (this->right_reads.size() == 0 || this->left_reads.size() == 0)) {
        bool best_pair = false;
        for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                if ((rv_idx->chr == lv_idx->chr) && (rv_idx->reversed != lv_idx->reversed)) {
                    insert1 = abs((double) lv_idx->get_end() - (double) rv_idx->start);
                    insert2 = abs((double) rv_idx->get_end() - (double) lv_idx->start);
                    if (insert1 <= (conf->insert_size * (1.0 + conf->insert_dev)) || insert2 <= (conf->insert_size * (1.0 + conf->insert_dev))) {
                        if (! best_pair) 
                            best_pair = (lv_idx->is_best && rv_idx->is_best);
                        active_left_reads.push_back(lv_idx);
                        active_right_reads.push_back(rv_idx);
                        found_pairs = true;
                    }
                }
            }
        }

        // no active pair is best alignment
        if (! best_pair && found_pairs) {
            for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
                if (lv_idx->is_best) {
                    lv_idx->is_best = false;
                    lv_idx->update_coverage_map(0);
                    break;
                }
            }
            active_left_reads.front()->is_best = true;
            active_left_reads.front()->update_coverage_map(1);

            pthread_mutex_lock(&mutex_best_left);
            genData->best_left[read_id] = (active_left_reads.front() - this->left_reads.begin());
            pthread_mutex_unlock(&mutex_best_left);

            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (rv_idx->is_best) {
                    rv_idx->is_best = false;
                    rv_idx->update_coverage_map(0);
                    break;
                }
            }
            active_right_reads.front()->is_best = true;
            active_right_reads.front()->update_coverage_map(1);

            pthread_mutex_lock(&mutex_best_right);
            assert(active_right_reads.front() - this->right_reads.begin() >= 0);
            genData->best_right[read_id] = (active_right_reads.front() - this->right_reads.begin());
            pthread_mutex_unlock(&mutex_best_right);
        }
        // did not find valid pairs
        if (! found_pairs) {
            active_left_reads.clear();
            for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
                if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                    continue;
                active_left_reads.push_back(lv_idx);
            }
            active_right_reads.clear();
            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                active_right_reads.push_back(rv_idx);
            }
        }
    } else {
        for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            active_left_reads.push_back(lv_idx);
        }
        for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
            if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                continue;
            active_right_reads.push_back(rv_idx);
        }
    }
}

char* OnlineData::parse_file(FILE* infile, char* last_line, GeneralData* genData, unsigned int &counter) {

    char line[1000] ;
    char cp_line[1000];
    char* ret = last_line;

    unsigned char pair_info = 0;

    Alignment curr_alignment;
    string id;
    this->last_id.clear();

    unordered_map <string, size_t, hash<string> >::iterator b_idx;

    this->left_reads.clear();
    this->right_reads.clear();

    while (true) {
        if (strlen(last_line) > 0 && strcmp(last_line, "samtools subprocess for reading terminated successfully\n")) { 
            strcpy(line, last_line);
            strcpy(last_line, "");
        }
        else {
            ret = fgets(line, sizeof(line), infile);
            counter++;

            if (conf->verbose && counter % 100000 == 0) 
                fprintf(stdout, "\n\t%i...", counter);
        }

        strcpy(cp_line, line);
        if (!ret) {
            break;
        }

        char* sl = strtok(line, "\t");

        id = curr_alignment.fill(sl, pair_info);

        if (id.size() == 0) {
            if (strcmp(cp_line, "samtools subprocess for reading terminated successfully\n")) {
                fprintf(stderr, "\nWARNING: SAM line incomplete! Ignoring line:\n%s\n", cp_line);
            } else {
                if (conf->verbose)
                    fprintf(stdout, "\n%s", cp_line);
            }
            continue ;
        }

        if (pair_info == 0) {
            // id == last_id or last_id is empty
            if ((! id.compare(this->last_id)) || this->last_id.size() == 0) {
                b_idx = genData->best_left.find(id);
                if (b_idx == genData->best_left.end()) {
                    genData->best_left.insert(pair<string, size_t>(id, this->left_reads.size()));
                    curr_alignment.is_best = true;
                }
                else if (b_idx->second == this->left_reads.size())
                    curr_alignment.is_best = true;
                this->left_reads.push_back(curr_alignment);
                this->last_id = id;
            } else {
                break ;
            }
        } else {
            // id == last_id or last_id is empty
            if ((! id.compare(this->last_id)) || this->last_id.size() == 0) {
                b_idx = genData->best_right.find(id);
                if (b_idx == genData->best_right.end()) {
                    genData->best_right.insert(pair<string, size_t>(id, this->right_reads.size()));
                    curr_alignment.is_best = true;
                }
                else if (b_idx->second == this->right_reads.size())
                    curr_alignment.is_best = true;
                this->right_reads.push_back(curr_alignment);
                this->last_id = id;
            } else {
                break ;
            }
        }
    }
    strcpy(last_line, cp_line);
    return ret;
}

