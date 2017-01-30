#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>

#include "superread_index.hpp"

#define QUEUE_SIZE 10000

using namespace std;

const char SR_DELIM[] = "_\n"; 
const string READ_DELIM = " \t\n";


int main(int argc, char **argv) {
    SuperreadIndex::init(8130595);
    
    string sr_name_line, sr_seq_line, read_ku_line;
    ifstream sr_name_file(argv[1]), sr_seq_file(argv[2]), read_ku_file(argv[3]);
    ofstream unmatched_out("unmatched_reads.txt");
    
    
    int sr_line = 0, sr_id;
    while (getline(sr_seq_file, sr_seq_line)) {
        sr_id = atoi(sr_seq_line.substr(1).c_str());
        while (sr_line <= sr_id) {
            getline(sr_name_file, sr_name_line);
            sr_line++;
        }
        SuperreadIndex::add_superread(sr_id, sr_name_line);
        getline(sr_seq_file, sr_seq_line);
    }
    
    //vector<Read*> read_list(QUEUE_SIZE, NULL);

    IndexQuery *query;
    query->reads.resize(QUEUE_SIZE);
    query->results = new ReadResults();
    
    Read *read1 = NULL, *read2 = NULL;

    vector< vector<int> > result;
    
    int i = 0;
    while(getline(read_ku_file, read_ku_line)) {
        if (read1 == NULL) {
            read1 = new Read(read_ku_line);
            continue;
        }

        read2 = new Read(read_ku_line);

        if (read1->is_mate(*read2)) {    
            query->reads[i] = read1;
            query->reads[i+1] = read2;
            read1 = read2 = NULL;
            i += 2;
        } else {
            if (read1 != NULL)
                unmatched_out << read1->read_id << endl;
            read1 = read2;
            read2 = NULL;
        }

        if (i >= query->reads.size()-1) {
            SuperreadIndex::match_pe_reads(query);

            for (int r = 0; r < query->results->size(); r++) {
                if (query->results->at(r).size() == 2) {
                    unmatched_out << query->results->at(r)[0] << endl << query->results->at(r)[1] << endl;
                } else {
                
                    printf("%d ", query->results->at(r)[0]);
                    for (int s = 2; s < query->results->at(r).size(); s++) {
                        printf("%d ", query->results->at(r)[s]);
                    }
                    printf("\n");
                    
                }
            }
            
            query->results->clear();
            query->reads.assign(QUEUE_SIZE, NULL);

            i = 0;
        }
    }

    unmatched_out.close();

    return 0;
}




