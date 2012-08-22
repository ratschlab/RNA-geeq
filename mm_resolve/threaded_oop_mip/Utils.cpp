#include <string>
#include <cstdio>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <assert.h>

#include "Alignment.h"
#include "GeneralData.h"
#include "config.h"

using namespace std;

extern GeneralData* genData;
extern Config* conf;

FILE* open_bam_pipe_in(std::string & in_fname)
{
	string command = string("/fml/ag-raetsch/share/software/samtools-0.1.7a/samtools view -h 2> /dev/null ") + in_fname  + " && echo samtools subprocess for reading terminated successfully";
	FILE* IN_FP=NULL ;
	fflush(stdout) ;
	IN_FP = popen(command.c_str(), "r") ;
	return IN_FP ;
}

FILE* open_bam_pipe_out(std::string & out_fname)
{
	string command = string("/fml/ag-raetsch/share/software/samtools-0.1.7a/samtools view -bS -o ") + out_fname +  " - 2> /dev/null" + " && echo samtools for writing subprocess terminated successfully";
	FILE* OUT_FP=NULL ;
	fflush(stdin) ;
	OUT_FP = popen(command.c_str(), "w") ;
	return OUT_FP ;
}

void parse_cigar(string cigar, vector<char> &operations, vector<int> &sizes) {
    
    size_t op_pos = cigar.find_first_of("HSMIDNXP");
    while (op_pos != string::npos) {
        sizes.push_back(atoi(cigar.substr(0, op_pos).c_str()));
        operations.push_back(cigar.at(op_pos));
        cigar = cigar.substr(op_pos + 1, cigar.size());
        op_pos = cigar.find_first_of("HSMIDNXP");
    }
}

double intron_penalty(vector<unsigned int> &intron_coverage) {
    double sum = 0.0;

    if (intron_coverage.size() > 0) {
        for (vector<unsigned int>::iterator it = intron_coverage.begin(); it != intron_coverage.end(); it++) {
            sum += (double) *it;
        }
        sum /= (double) intron_coverage.size();
    }

    return sum;
}

double get_variance(vector<unsigned int> &exon_coverage, vector<unsigned int> &intron_coverage) {

    if (exon_coverage.size() < 2) {
        return 0.0;
    } else {
        double sum = 0.0;
        double int_pnlty = 0.0;
/*        for ( vector<unsigned short>::iterator it = exon_coverage.begin(); it != exon_coverage.end(); it++) {
            sum += (double) *it;
            fprintf(stdout, "%i ", *it);
        }
        fprintf(stdout, "\n");*/
        double mean = sum / (double) exon_coverage.size();
        sum = 0.0;
        for ( vector<unsigned int>::iterator it = exon_coverage.begin(); it != exon_coverage.end(); it++) {
            sum += pow((double) *it - mean, 2.0);
        }

        if (intron_coverage.size() > 0) {
            int_pnlty = intron_penalty(intron_coverage);
        }
        
        //fprintf(stdout, "var: %f, intron: %f\n", (sum / ((double) exon_coverage.size() - 1.0)), int_pnlty);
        // TODO
        return (sum / ((double) exon_coverage.size() - 1.0)) + int_pnlty;  
    }
}

vector<unsigned int> alter_coverage(vector<unsigned int> &source, unsigned int window_left, unsigned int window_right, bool is_positive) {

    vector<unsigned int> target;

    for (size_t i = 0; i < source.size(); i++) {
        if (i < window_left || i >= source.size() - window_right) {
            target.push_back(source.at(i));
        } else {
            target.push_back((source.at(i) > 0 || is_positive)?source.at(i) + (2*is_positive) - 1:0);
        }
    }
    return target;
}

string update_line_flag(char* line, bool is_best) {

    char cp_line[1000];
    strcpy(cp_line, line);

    char* sl = strtok(line, "\t");

    int idx = 0;
    string return_line;
    while (sl != NULL) {
        if (idx == 1) {
            if ( (!is_best && (atoi(sl) & 256) == 256) || (is_best && (atoi(sl) & 256) == 0)) { 
                return_line = string(cp_line);
                return return_line.substr(0, return_line.size() - 1);
            }
            else {
                char new_flag[10]; 
                if (is_best) {
                    sprintf(new_flag, "%i", atoi(sl) - 256);        
                }
                else {
                    sprintf(new_flag, "%i", atoi(sl) + 256);        
                }
                return_line += (string("\t") + new_flag);
            }
        }
        else {
            return_line += (string("\t") + sl);
        }
        sl = strtok(NULL, "\t");
        idx++;
    }
    return return_line.substr(1, return_line.size());
}


void parse_header(char* sl) {

    int idx = 0;
    string chr_name;
    while(sl != NULL) {
        string tmp_sl = string(sl);
        
        if (idx == 1) {
            tmp_sl = tmp_sl.substr(3, tmp_sl.size());
            chr_name = tmp_sl ;
            if (genData->chr_num.find(tmp_sl) == genData->chr_num.end()) {
                genData->chr_num.insert( pair<string, int>(tmp_sl,  genData->chr_num.size() + 1) );
            }
            else {
                fprintf(stderr, "WARNING: Doubled contig names in header!\n Ignoring %s\n\n", tmp_sl.c_str());
            }
        }
        else if (idx == 2) {
            tmp_sl = tmp_sl.substr(3, tmp_sl.size());
            genData->chr_size.push_back(atoi(tmp_sl.c_str()));    
            vector<unsigned int> tmp_cov(genData->chr_size.back(), 0);
            if (conf->verbose) 
                fprintf(stdout, "\t...reserving memory for contig %s (+) of size %i\n", chr_name.c_str(), genData->chr_size.back());
            genData->coverage_map.insert( pair<pair<unsigned char, unsigned char>, vector<unsigned int> >(pair<unsigned char, unsigned char>(genData->chr_num[chr_name], '+'), tmp_cov) ); 
            if (conf->strand_specific) {
                if (conf->verbose) 
                    fprintf(stdout, "\t...reserving memory for contig %s (-) of size %i\n", chr_name.c_str(), genData->chr_size.back());
                genData->coverage_map.insert( pair<pair<unsigned char, unsigned char>, vector<unsigned int> >(pair<unsigned char, unsigned char>(genData->chr_num[chr_name], '-'), tmp_cov) ); 
            }
        }   
        idx ++;
        sl = strtok(NULL, "\t");    
    }

    if (genData->chr_size.size() != genData->chr_num.size() || (!conf->strand_specific && genData->chr_num.size() != genData->coverage_map.size()) || (conf->strand_specific && 2*genData->chr_num.size() != genData->coverage_map.size())) {
        fprintf(stderr, "\nERROR: Header information incomplete!");
        exit(-1);
    }
}

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right, double &loss, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    // compute overlap between candidate alignments
    set<unsigned long> already_covered_pos;
    set<unsigned long> not_covered_pos;
    set<unsigned long> empty_set;
    set<unsigned long> overlap;
    set<unsigned long> genome_pos_cl = candidate_left->get_genome_pos();
    set<unsigned long> genome_pos_cr = candidate_right->get_genome_pos();
    set<unsigned long> genome_pos_bl = best_left->get_genome_pos();
    set<unsigned long> genome_pos_br = best_right->get_genome_pos();

    pair<double, double> loss_cl, loss_cr, loss_bl, loss_br;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {
       
        // loss candidate_left
        // loss best_left
        if (candidate_left == best_left) {
            loss_cl = pair<double, double>(0.0, 0.0);
            loss_bl = pair<double, double>(0.0, 0.0);
        } else {
            if (candidate_left->chr == best_left->chr)
                set_intersection(genome_pos_bl.begin(), genome_pos_bl.end(), genome_pos_cl.begin(), genome_pos_cl.end(), inserter(overlap, overlap.begin()));

            loss_cl = genData->segments.get_exon_segment_loss(candidate_left, overlap, debug);
            loss_bl = genData->segments.get_exon_segment_loss(best_left, overlap, debug);
        }
        overlap.clear();

        // loss candidate_rigth
        // loss best_right
        if (candidate_right == best_right) {
            loss_cr = pair<double, double>(0.0, 0.0);
            loss_br = pair<double, double>(0.0, 0.0);
        } else {
            if (candidate_right->chr == best_right->chr)
                set_intersection(genome_pos_br.begin(), genome_pos_br.end(), genome_pos_cr.begin(), genome_pos_cr.end(), inserter(overlap, overlap.begin()));
            loss_cr = genData->segments.get_exon_segment_loss(candidate_right, overlap, debug);
            loss_br = genData->segments.get_exon_segment_loss(best_right, overlap, debug);
        }
        overlap.clear();

        // check if any loss is valid
        // first is loss_with and second is loss_without
        if (loss_cl.first >= 0.0 || loss_cr.first >= 0.0 || loss_br.first >= 0.0 || loss_bl.first >= 0.0) {
           used_mip = true;
        }
        if (debug) {
            fprintf(stdout, "candidate left:\n");
            candidate_left->print();
            fprintf(stdout, "candidate right:\n");
            candidate_right->print();
            fprintf(stdout, "best left:\n");
            best_left->print();
            fprintf(stdout, "best right:\n");
            best_right->print();
            fprintf(stdout, "cand loss left (w): %f  cand loss right (w): %f best loss left (w): %f best loss right (w):%f\n", loss_cl.first, loss_cr.first, loss_bl.first, loss_br.first);
            fprintf(stdout, "cand loss left (wo): %f  cand loss right (wo): %f best loss left (wo): %f best loss right (wo):%f\n", loss_cl.second, loss_cr.second, loss_bl.second, loss_br.second);
        }
    }

    // if we never were supposed to use the mip objective or did not use it
    // for other reasons, use the variance objective instead, but do not add it
    // to the total loss for mip-objective
    if (! conf->use_mip_objective || ! used_mip) {

        // loss candidate_left
        set_intersection(genome_pos_bl.begin(), genome_pos_bl.end(), genome_pos_br.begin(), genome_pos_br.end(), inserter(already_covered_pos, already_covered_pos.begin()));
        loss_cl = candidate_left->get_variance_loss(already_covered_pos, empty_set);

        // loss candidate_rigth
        loss_cr = candidate_right->get_variance_loss(already_covered_pos, empty_set);
        already_covered_pos.clear();
        not_covered_pos.clear();

        // loss best_left and loss best_right
        assert(!candidate_left->is_best || !candidate_right->is_best);
        if (!candidate_left->is_best && !candidate_right->is_best) {
            set_intersection(genome_pos_cl.begin(), genome_pos_cl.end(), genome_pos_cr.begin(), genome_pos_cr.end(), inserter(not_covered_pos, not_covered_pos.begin()));
            loss_bl = best_left->get_variance_loss(empty_set, not_covered_pos);
            loss_br = best_right->get_variance_loss(empty_set, not_covered_pos);
        } else if (candidate_left->is_best) {
            loss_bl = best_left->get_variance_loss(genome_pos_cl, genome_pos_cr);
            loss_br = best_right->get_variance_loss(genome_pos_cl, genome_pos_cr);
        } else if (candidate_right->is_best) {
            loss_bl = best_left->get_variance_loss(genome_pos_cr, genome_pos_cl);
            loss_br = best_right->get_variance_loss(genome_pos_cr, genome_pos_cl);
        }
    }

    // determine total loss
    candidate_loss += (loss_cl.first >= 0.0) ? loss_cl.first : 0.0;
    candidate_loss += (loss_cr.first >= 0.0) ? loss_cr.first : 0.0;
    candidate_loss += (loss_bl.second >= 0.0) ? loss_bl.second : 0.0;
    candidate_loss += (loss_br.second >= 0.0) ? loss_br.second : 0.0;
    best_loss += (loss_cl.second >= 0.0) ? loss_cl.second : 0.0;
    best_loss += (loss_cr.second >= 0.0) ? loss_cr.second : 0.0;
    best_loss += (loss_bl.first >= 0.0) ? loss_bl.first : 0.0;
    best_loss += (loss_br.first >= 0.0) ? loss_br.first : 0.0;
    if (debug) {
        fprintf(stdout, "best loss: %f\n", best_loss);
        fprintf(stdout, "candidate loss: %f\n", candidate_loss);
        fprintf(stdout, "delta (best - cand): %f\n", best_loss - candidate_loss);
    }
 
    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
    } else {
        loss = min(candidate_loss, best_loss);
    }
    return (candidate_loss < best_loss);
}

bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best, double &loss, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    // compute overlap between candidate and best alignments
    set<unsigned long> genome_pos_candidate = candidate->get_genome_pos();
    set<unsigned long> genome_pos_best = best->get_genome_pos();
    set<unsigned long> empty_set;

    pair<double, double> loss_candidate, loss_best;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {

        loss_candidate = genData->segments.get_exon_segment_loss(candidate, genome_pos_best, debug);
        loss_best = genData->segments.get_exon_segment_loss(best, genome_pos_candidate, debug);

        // check if any loss is valid
        // first is loss_with and second is loss_without
        if (loss_candidate.first >= 0.0 || loss_best.first >= 0.0) {
            used_mip = true;
        }
        if (debug) {
            fprintf(stdout, "candidate (single):\n");
            candidate->print();
            fprintf(stdout, "best (single):\n");
            best->print();
            fprintf(stdout, "cand loss (w): %f best loss (w): %f\n", loss_candidate.first, loss_best.first);
            fprintf(stdout, "cand loss (wo): %f best loss (wo): %f\n", loss_candidate.second, loss_best.second);
        }

    }

    // if we never were supposed to use the mip objective or did not use it
    // for other reasons, use the variance objective instead, but do not add it
    // to the total loss for mip-objective
    if (! conf->use_mip_objective || ! used_mip) {

        loss_candidate = candidate->get_variance_loss(genome_pos_best, empty_set);
        loss_best = best->get_variance_loss(empty_set, genome_pos_candidate);
    }

    candidate_loss +=  (loss_candidate.first >= 0.0) ? loss_candidate.first : 0;
    candidate_loss +=  (loss_best.second >= 0.0) ? loss_best.second : 0;
    best_loss +=  (loss_candidate.second >= 0.0) ? loss_candidate.second : 0;
    best_loss +=  (loss_best.first >= 0.0) ? loss_best.first : 0;

    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
    } else {
        loss = min(candidate_loss, best_loss);
    }
//    fprintf(stderr, "single cand loss: %f  single best loss: %f\n", candidate_loss, best_loss);
    return (candidate_loss < best_loss);
}

set<vector<Alignment>::iterator> filter_alignments(vector<Alignment> &aligns) {

    set<vector<Alignment>::iterator> to_erase;
    vector<Alignment>::iterator v_idx;
    //sort(aligns.begin(), aligns.end(), alignment_comparator);
    if (aligns.size() == 0)
        return to_erase;
    unsigned int min_ops = aligns.begin()->edit_ops;
    for (v_idx = aligns.begin() + 1; v_idx < aligns.end(); v_idx++)
        min_ops = v_idx->edit_ops < min_ops?v_idx->edit_ops:min_ops;

    for (v_idx = aligns.begin(); v_idx < aligns.end(); v_idx++) {
        if (v_idx->edit_ops - min_ops > conf->filter_distance) {
            to_erase.insert(v_idx);
        }
    }
    return to_erase;
}

void get_plifs_from_file() {
    
    FILE* infile = fopen(conf->lossfile.c_str(), "r");    
	if (!infile){
		fprintf(stderr, "Could not open file %s for reading!\n", conf->lossfile.c_str());
		exit(2);
	}
    char* ret;
    char line[1000];
    int i_count;
    float c, sl1, sl2, sr1, sr2;

    map< double, vector<double> > plifs;

    while (true) {
        ret = fgets(line, sizeof(line), infile);

        if (!ret)
            break;
        
        if (line[0] == '#')
            continue;

        i_count = sscanf(line, "%f %f %f %f %f", &c, &sl1, &sl2, &sr1, &sr2);

        if (i_count != 5) {
            fprintf(stderr, "Number of items per line in PLIF-Input file is different from 5! - Bail out!\n");
            exit(1);
        }

        vector<double> tmp_plif;
        tmp_plif.push_back((double) sl1);
        tmp_plif.push_back((double) sl2);
        tmp_plif.push_back((double) sr1);
        tmp_plif.push_back((double) sr2);
        plifs.insert(plifs.begin(), pair<double, vector<double> >((double) c, tmp_plif));
    }

    genData->plifs = plifs;
    
}

double compute_mip_loss(double observed_cov, double predicted_cov, unsigned long segment_len) {


    // in the case of exonic segments, the observed coverage is a count value and the
    // predicted coverage is a mean that needs to be transfered into counts
    // in the case of intronic segments, the coverages are already counts, segment_len is 0 then 
    if (segment_len > 0 && conf->read_len > 0) {
        predicted_cov *= segment_len / conf->read_len; 
        observed_cov /= conf->read_len; 
    }

    if (observed_cov > 30000)
        observed_cov = 30000;
    if (observed_cov < 0)
        observed_cov = 0;
    
    // get plif iterators to interpolate entries
    map< double, vector<double> >::iterator upper = genData->plifs.lower_bound(observed_cov);
    map< double, vector<double> >::iterator lower = upper;
    if (observed_cov > 0)
        lower--;

    double loss = 0.0;

    double diff = abs(predicted_cov - observed_cov);

    // determine if left or right part of function needs to be queried
    if (predicted_cov > observed_cov) {
        //loss = ((lower->second.at(0) + upper->second.at(0)) / 2)*predicted_cov*predicted_cov + ((lower->second.at(1) + upper->second.at(1)) / 2)*predicted_cov; 
        loss = ((lower->second.at(0) + upper->second.at(0)) / 2)*diff*diff + ((lower->second.at(1) + upper->second.at(1)) / 2)*diff; 
    } else {
        //loss = ((lower->second.at(2) + upper->second.at(2)) / 2)*predicted_cov*predicted_cov + ((lower->second.at(3) + upper->second.at(3)) / 2)*predicted_cov; 
        loss = ((lower->second.at(2) + upper->second.at(2)) / 2)*diff*diff + ((lower->second.at(3) + upper->second.at(3)) / 2)*diff; 
    }

    return loss;
}

void prepare_mip_objective() {

        if (conf->segmentfile.size() == 0) {
            fprintf(stderr, "For using the mip objective, you need to specify a segments file via -s/--segmentfile. Bailing out!\n");
            exit(-1);
        } else {
            if (conf->verbose) fprintf(stdout, "Parsing segments from %s ...\n", conf->segmentfile.c_str());
            genData->segments.get_from_file(); 
            if (conf->verbose) fprintf(stdout, "... done.\n\n");
        }

       // fprintf(stdout, "chr_size:%i\tchr_num:%i\tcov_map:%i\n", genData->chr_num.size(), genData->chr_num.size(), genData->coverage_map.size());

        if (conf->lossfile.size() == 0) {
            fprintf(stderr, "For using the mip objective, you need to specify a loss function parameter file via -l/--lossfile. Bailing out!\n");
            exit(-1);
        } else {
            if (conf->verbose) fprintf(stdout, "Parsing loss function parameters from %s ...\n", conf->lossfile.c_str());
            get_plifs_from_file();
            if (conf->verbose) fprintf(stdout, "... done.\n\n");
        }
}
