#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>

#include <timer.h>
#include <src/superread_index.hpp>

#define BASE_WEIGHT 1.0

const std::string SR_DELIM = "_\n"; 
const std::string READ_DELIM = " \t\n";

Read::Read(std::string read_str) {
    char *cstr = new char[read_str.length()+1];
    strcpy(cstr, read_str.c_str());
    
    if (read_str.length() == 0) {
        is_valid = false;
        return;
    }


    int delim_ct = 0;
    for (unsigned int i = 0; i < read_str.length(); i++) {
        for (unsigned int j = 0; j < READ_DELIM.size(); j++) {
            if (cstr[i] == READ_DELIM[j]) {
                delim_ct++;
                break;
            }
        }
    }

    ku_ids.reserve((delim_ct-1)/3);
    ku_strands.reserve(ku_ids.capacity());
    
    char *read_name = strtok(cstr, READ_DELIM.c_str());

    if (strlen(read_name) < PREFIX_LEN+1) {
        is_valid = false;
        return;
    }

    for (unsigned int i = 0; i < PREFIX_LEN; i++) {
        prefix[i] = read_name[i];
    }
    read_id = atoi(&read_name[PREFIX_LEN]);

    length = atoi(strtok(NULL, READ_DELIM.c_str()));
    
    char *ku_id;

    while ( (ku_id = strtok(NULL, READ_DELIM.c_str())) != NULL ) {
        ku_ids.push_back(atoi(ku_id));
        strtok(NULL, READ_DELIM.c_str());
        ku_strands.push_back(*strtok(NULL, READ_DELIM.c_str()));
    }

    delete[] cstr;

    is_valid = true;
}

Read::Read(const Read &r) {
    read_id = r.read_id;
    ku_ids = r.ku_ids;
    ku_strands = r.ku_strands;
    for (unsigned int i = 0; i < PREFIX_LEN; i++)
        prefix[i] = r.prefix[i];
}

bool Read::is_mate(Read r) {
    if (read_id % 2 == 0)
        return read_id+1 == r.read_id;
    return read_id-1 == r.read_id;
}

bool Read::same_prefix(Read r) {
    for (unsigned int i = 0; i < PREFIX_LEN; i++)
        if (prefix[i] != r.prefix[i])
            return false;
    return true;
}

IndexQuery::IndexQuery(int buffer_size, bool paired) {
    reads.reserve(buffer_size);
    results = new std::vector< std::vector<int>* >();    
    matched_reads = new std::vector<Read*>();
    finished = true;
}

IndexQuery::~IndexQuery() {
    delete results;
    delete matched_reads;
}

bool IndexQuery::add_read(Read *read1, Read *read2) {
    if (is_full())
        return false;

    if (!is_empty() && !reads.back()->same_prefix(*read1))
        return false;

    reads.push_back(read1);

    if (read2 != NULL)
        reads.push_back(read2);

    return true;
}

bool IndexQuery::is_full() {
    return reads.size() >= reads.capacity()-1;
}

bool IndexQuery::is_empty() {
    return reads.size() == 0;
}

void IndexQuery::clear() {
    for (unsigned int r = 0; r < results->size(); r++)
        delete results->at(r);                                    
    results->clear();
    for (unsigned int i = 0; i < matched_reads->size(); i++)
        delete matched_reads->at(i);
    matched_reads->clear();
    reads.clear();
}


ReadIO::ReadIO(std::string read_ku_fname, std::vector<std::string> read_files) {
    read_ku_file.open(read_ku_fname.c_str());
    read1 = NULL;
    read2 = NULL;

    prefixes.resize((PREFIX_LEN+1)*read_files.size(), ' ');
    is_paired.resize(read_files.size());
    fasta_fnames.resize(read_files.size());
    matched_reads.resize(read_files.size());

    for (unsigned int i = 0; i < read_files.size(); i++) {
        std::string s, file1, file2;
       
        //Generate std::string "p1 p2 p3 " where each pn is a prefix
        //Used to index prefixes
        for (unsigned int j = 0; j < PREFIX_LEN; j++) 
            prefixes[(PREFIX_LEN+1)*i+j] = read_files[i][j];

        if (read_files[i].find(" ", PREFIX_LEN+1) != std::string::npos) {
            is_paired[i] = true;
        } else {
            is_paired[i] = false;
        }

        
        fasta_fnames[i] = read_files[i].substr(PREFIX_LEN+1);
    }
}

int ReadIO::get_prefix_id(Read r) {
    for (unsigned int i = 0; i < prefixes.size(); i += PREFIX_LEN+1) {
        for (unsigned int j = 0; j < PREFIX_LEN; j++) {
            if (prefixes[i+j] == r.prefix[j]) {
                if (j == PREFIX_LEN-1) {
                    return i / (PREFIX_LEN+1);
                }
            } else {
                break;
            }
        }
    }

    return -1;
}

ReadIO::~ReadIO() {
    if (read1 != NULL)
        delete read1;
    if (read2 != NULL)
        delete read2;
}

bool ReadIO::fill_buffer(IndexQuery *q) {
    std::string read_ku_line;

    while (read1 == NULL) {
        if (!getline(read_ku_file, read_ku_line))
            return false;
        read1 = new Read(read_ku_line);
        if (!read1->is_valid)
            read1 = NULL;
    }

    q->prefix_id = get_prefix_id(*read1);
    q->is_paired = is_paired[q->prefix_id];

    if (is_paired[q->prefix_id]) {
        if (read2 != NULL) {
            q->add_read(read1, read2);
            read1 = read2 = NULL;
        }
    } else {
        q->add_read(read1);
        if (read2 != NULL)
            q->add_read(read2);
        read1 = read2 = NULL;
    }
        

    //Fill read buffer  
    while (!(q->is_full()) && getline(read_ku_file, read_ku_line)) {

        if (is_paired[q->prefix_id]) {  

            //Set mate 1 if unset, then move on to next read                    
            if (read1 == NULL) {                                                
                read1 = new Read(read_ku_line);
                if (!read1->is_valid)
                    read1 = NULL;
                continue;                                                       
            }                                                                   

            //Read mate 2                                                       
            read2 = new Read(read_ku_line);
            if (!read2->is_valid) {
                read2 = NULL;
                continue;
            }

            //Add reads to buffer if they are mates                             
            if (read1->is_mate(*read2)) {  

                //Only add reads with same prefixes
                if (!q->add_read(read1, read2))
                    break;

                read1 = read2 = NULL;                                           

            //read1 has no mate                                             
            } else {
                delete read1;                                               
                read1 = read2;                                                  
                read2 = NULL;                                                   
            }   
 
        //Single-end reads
        } else {
            if (read1 == NULL)
                read1 = new Read(read_ku_line);

            if (!read1->is_valid) {
                read1 = NULL;
                continue;
            }

            if (!q->add_read(read1))
                break;

            read1 = NULL;
        }
    }
   
    if (q->is_empty())
        return false;

    q->finished = false;

    return true;
}

void ReadIO::update_matched(IndexQuery *q) {
    if (q->prefix_id < 0) {
        std::cerr << "Error: empty buffer\n";
        return;
    }
    
    int step = is_paired[q->prefix_id] ? 2 : 1;

    for (unsigned int i = 0; i < q->matched_reads->size(); i += step) 
        matched_reads[q->prefix_id].push_back(q->matched_reads->at(i)->read_id);
    
}

void ReadIO::write_unmatched(std::string pe1_fname, 
                             std::string pe2_fname, 
                             std::string se_fname) {

    std::ofstream pe1_out, pe2_out, se_out;
    std::ifstream pe1_in, pe2_in, se_in;

    for (unsigned int i = 0; i < fasta_fnames.size(); i++) {
        
        std::sort(matched_reads[i].begin(), matched_reads[i].end());

        unsigned int j = 0;
        int read_id = 0;

        if (is_paired[i]) {

            if (!pe1_out.is_open()) {
                pe1_out.open(pe1_fname.c_str());
                pe2_out.open(pe2_fname.c_str());
            }

            unsigned int split = fasta_fnames[i].find(" ");
            pe1_in.open(fasta_fnames[i].substr(0, split).c_str());
            pe2_in.open(fasta_fnames[i].substr(split+1).c_str());
        
            while (j < matched_reads[i].size()) {

                if (matched_reads[i][j] == read_id) {
                    skip_fastq(pe1_in);
                    skip_fastq(pe2_in);
                    j++;
                } else {
                    copy_fastq(pe1_in, pe1_out);
                    copy_fastq(pe2_in, pe2_out);
                }

                read_id += 2;
            }

            pe1_in.close();
            pe2_in.close();

        } else {

            if (!se_out.is_open()) 
                se_out.open(se_fname.c_str());
            
            se_in.open(fasta_fnames[i].c_str());

            while (j < matched_reads[i].size()) {

                if (matched_reads[i][j] == read_id) {
                    skip_fastq(se_in);
                    j++;
                } else {
                    copy_fastq(se_in, se_out);
                }

                read_id += 1;
            }

            se_in.close();
        }

    }

    if (pe1_out.is_open()) {
        pe1_out.close();
        pe2_out.close();
    }

    if (se_out.is_open())
        se_out.close();
}

bool ReadIO::copy_fastq(std::ifstream &in, std::ofstream &out) {
    std::string s;

    if (!getline(in, s))
        return false;
    
    out << s << "\n";

    for (unsigned int i = 0; i < 3; i++) {
        getline(in, s);
        out << s << "\n";
    }
    
    return true;
}

bool ReadIO::skip_fastq(std::ifstream &in) {
    bool not_eof;
    std::string s;

    for (unsigned int i = 0; i < 4; i++) 
        not_eof = getline(in, s);
    
    return not_eof;
}

std::vector< std::vector<SuperreadIndex::KUnitig> > SuperreadIndex::ku_list;
std::vector<int> SuperreadIndex::sr_lens;
std::vector<int> SuperreadIndex::sr_ids;
std::map< std::vector<int>, int > SuperreadIndex::sr_frag_counts;
std::vector<double> SuperreadIndex::weights;
std::string SuperreadIndex::sr_fasta_fname;

bool SuperreadIndex::load_files(std::string sr_name_path, 
                                std::string sr_fasta_path, 
                                int max_ku) {

    //Init superread index
    if (max_ku > 0) 
        ku_list.resize(max_ku);
    
    //Save fasta filename for later output
    sr_fasta_fname = sr_fasta_path;

    //Only want to index SRs present in sr_seq_file
    //SR KUs read from sr_name_file
    std::ifstream sr_name_file(sr_name_path.c_str()), 
             sr_seq_file(sr_fasta_path.c_str());
    std::string sr_name_line, sr_seq_line;

    if (!sr_name_file.is_open()) {
        std::cerr << "Error: unable to open '" << sr_name_path << "'\nExiting\n";
        return false;
    } else if (!sr_seq_file.is_open()) {
        std::cerr << "Error: unable to open '" << sr_fasta_path << "'\nExiting\n";
        return false;
    }

    //Iterate thru sr_seq file to find SRs to store
    int sr_line = 0, sr_id;
    while (getline(sr_seq_file, sr_seq_line)) {
        sr_id = atoi(sr_seq_line.substr(1).c_str());

        //Find next superread K Unitigs
        while (sr_line <= sr_id) {
            getline(sr_name_file, sr_name_line);
            sr_line++;
        }

        ///Add superread to ku index
        add_superread_kus(sr_ids.size(), sr_name_line);
       
        //Store SR length
        getline(sr_seq_file, sr_seq_line);
        sr_lens.push_back(sr_seq_line.size());        

        sr_ids.push_back(sr_id);
    }

    weights = std::vector<double>(sr_lens.size(), 0.0);


    sr_name_file.close();
    sr_seq_file.close();

    return true;
}
 
void SuperreadIndex::destroy() {
    for (unsigned int i = 0; i < ku_list.size(); i++) {
        for (unsigned int j = 0; j < ku_list[i].size(); j++) {
            delete ku_list[i][j];
        }
        ku_list[i].clear();
    }
    ku_list.clear();
}

int SuperreadIndex::match_reads(IndexQuery *q) {
    if (q->is_paired)
        match_paired_reads(q);
    else
        match_unpaired_reads(q);

    return 0;
}

void SuperreadIndex::update_frag_counts(IndexQuery *q) {
    for (unsigned int i = 0; i < q->results->size(); i++) {
        std::vector<int> &sr_list = *(q->results->at(i));

        int read_len;
        if (q->is_paired) {
            read_len = q->matched_reads->at(2*i)->length 
                     + q->matched_reads->at(2*i+1)->length;
        } else {
            read_len = q->matched_reads->at(i)->length;
        }

        //std::cout << read_len << "\n";

        if (sr_frag_counts.count(sr_list) == 0) {
            sr_frag_counts[sr_list] = read_len;
        } else {
            sr_frag_counts[sr_list] += read_len;
        }

    
        if (sr_list.size() == 1) { //??
            weights[sr_list[0]] += (double) read_len / (double) sr_lens[sr_list[0]];
        }
    }
}

double SuperreadIndex::iter_weights() {
    int count;
    double weight_sum;
    std::vector<double> new_weights(weights.size(), 0.0);

    std::map<std::vector<int>, int>::iterator frag_iter = sr_frag_counts.begin();
    for (; frag_iter != sr_frag_counts.end(); ++frag_iter) {
        const std::vector<int> &sr_list = frag_iter->first;
        count = frag_iter->second;
        
        weight_sum = 0;
        for (unsigned int i = 0; i < sr_list.size(); i++)
            weight_sum += sr_weight_share(sr_list[i]);

        for (unsigned int i = 0; i < sr_list.size(); i++) {
            new_weights[sr_list[i]] += count*sr_weight_share(sr_list[i])/weight_sum;
        }
    }
    
    double diff, max_diff = 0;
    for (unsigned int i = 0; i < weights.size(); i++) {
        
        //Normalize
        new_weights[i] /= double(sr_lens[i]);

        diff = weights[i]>new_weights[i] ? 
                  weights[i]-new_weights[i] 
                : new_weights[i]-weights[i];

        if (diff > max_diff)
            max_diff = diff;
    }

    weights.swap(new_weights);

    return max_diff;
}

void SuperreadIndex::write_fastq(std::string fname, std::string delim) {
    std::ifstream sr_seq_in(sr_fasta_fname.c_str());
    std::ofstream sr_seq_out(fname.c_str());

    std::string sr_header, sr_seq;

    int i = 0;
    while (getline(sr_seq_in, sr_header)) {
        getline(sr_seq_in, sr_seq);
        
        if (weights[i] > 0) {
            sr_header[0] = '@';

            std::string sr_quals(sr_seq.size(), 'J');
            sr_seq_out << sr_header << delim << weights[i] << "\n"
                       << sr_seq << "\n+\n" << sr_quals << "\n";
        }
        i++;
    }
}

//Private functions

void SuperreadIndex::add_superread_kus(int sr_index, std::string sr_str) {
    char *cstr = new char[sr_str.length()+1];
    strcpy(cstr, sr_str.c_str());

    KUnitig ku;

    char *ku_str;   //Stores each std::string representing a KU (ex 12345F)
    unsigned int ku_id,      //Stores int of KU (ex 12345)
        ku_len,     //Stores length of KU std::string
        sr_pos = 0; //Stores position order of KU within superread

    ku_str = strtok(cstr, SR_DELIM.c_str());

    do {
        ku_len = strlen(ku_str);

        //Create new K-Unitig record
        ku = new struct kunitig;
        ku->sr = sr_index;
        ku->strand = ku_str[ku_len-1]; //Strand is last char
        ku->sr_pos = sr_pos;
        
        ku_str[ku_len] = '\0'; //Remove strand char
        ku_id = atoi(ku_str);

        //Make room for KU (if needed) and store it
        if (ku_id > ku_list.size())
            ku_list.resize(ku_id+1);
        ku_list[ku_id].push_back(ku);
        
        sr_pos++;

    } while ( (ku_str = strtok(NULL, SR_DELIM.c_str())) != NULL );

    delete[] cstr;
}

void SuperreadIndex::match_paired_reads(IndexQuery *q) { 
    
    #ifdef DEBUG_TIMER
    Timer timer;
    #endif

    std::vector<int> read1_srs, read2_srs;

    int reads_checked = 0;

    Read *read1, *read2;
    std::vector<int> *joined;
    for (unsigned int r = 0; r < q->reads.size(); r += 2) {
        read1 = q->reads[r];
        read2 = q->reads[r+1];
        
        read1_srs = match_read(read1);
        read2_srs = match_read(read2);
        
        reads_checked += 2;
        
        joined = new std::vector<int>();

        unsigned int i = 0, j = 0;
        while (i < read1_srs.size() && j < read2_srs.size()) {
            if (read1_srs[i] < read2_srs[j]) {
                i++;
            } else if (read1_srs[i] > read2_srs[j]) {
                j++;
            } else {
                joined->push_back(read1_srs[i]);
                i++;
                j++;
            }
        }

        if (joined->size() > 0) {
            q->results->push_back( joined );
            q->matched_reads->push_back( read1 );
            q->matched_reads->push_back( read2 );
        } else {
            delete joined;
            delete read1;
            delete read2;
        }

        read1_srs.clear();
        read2_srs.clear();
        
    }

    q->reads.clear();
    
    #ifdef DEBUG_TIMER
    q->thread_time = timer.get();
    #endif
     
    q->finished = true;
}

void SuperreadIndex::match_unpaired_reads(IndexQuery *q) { 
    
    #ifdef DEBUG_TIMER
    Timer timer;
    #endif

    std::vector<int> *read_srs;

    Read *read;
    for (unsigned int r = 0; r < q->reads.size(); r += 1) {
        read = q->reads[r];
        read_srs = new std::vector<int>(match_read(read));
        
        if (read_srs->size() > 0) {
            q->results->push_back( read_srs );
            q->matched_reads->push_back( read );
        } else {
            delete read_srs;
            delete read;
        }

    }

    q->reads.clear();
    
    #ifdef DEBUG_TIMER
    q->thread_time = timer.get();
    #endif
     
    q->finished = true;
}

std::vector<int> SuperreadIndex::match_read(Read *read) {
    std::vector<KUnitig> prev_kus, next_kus;
    prev_kus = ku_list[read->ku_ids[0]]; //
    next_kus.reserve(prev_kus.size());
    
    //Loops thru every read KU
    for (unsigned int ku = 1; ku < read->ku_ids.size(); ku++) {
        std::vector<KUnitig> &sr_kus = ku_list[read->ku_ids[ku]];
        unsigned int i, j, k;
        i = j = 0;
        while (i < prev_kus.size() && j < sr_kus.size()) {
            if (prev_kus[i]->sr < sr_kus[j]->sr) {
                i++;
            } else if (prev_kus[i]->sr > sr_kus[j]->sr) {
                j++;
            } else {
                k = j;
                while (k < sr_kus.size() && sr_kus[k]->sr == prev_kus[i]->sr) {
                    int pos_shift = 1;
                    if (sr_kus[k]->strand == read->ku_strands[ku]) 
                        pos_shift = -1;

                    if (prev_kus[i]->sr_pos == sr_kus[k]->sr_pos + pos_shift) {
                        next_kus.push_back(sr_kus[k]);
                        break;
                    }
                    k++;
                }
                if (k == j)
                    j++;
                i++;
            }
        }

        if (next_kus.size() == 0) {
            prev_kus.clear();
            break;
        }

        prev_kus.swap(next_kus);
        next_kus.clear();
    }
    
    std::vector<int> result(prev_kus.size());
    for (unsigned int i = 0; i < prev_kus.size(); i++)
        result[i] = prev_kus[i]->sr;
    return result;

}

inline double SuperreadIndex::sr_weight_share(int sr) {
    return BASE_WEIGHT + weights[sr];
}
