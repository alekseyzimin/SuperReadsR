#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <unistd.h>

#include <thread_pool.hpp>
#include <src/superread_index.hpp>
#include <timer.h>


//Command line arguments
const char OPTIONS[] = "w:r:o:p:b:k:";
const char WORK_DIR_OPT       = 'w',
           READ_KU_OPT        = 'r',
           OUT_PREFIX_OPT     = 'o',
           NUM_THREADS_OPT    = 'p',
           BUFFER_SIZE_OPT    = 'b',
           KU_NUM_OPT         = 'k'; //8130595

//Argument defaults
const std::string DEF_WORK_DIR = "./work1/",
             DEF_READ_KU_FNAME = "newTestOutput.nucmerLinesOnly",
             SR_FASTA_FNAME    = "superReadSequences.fasta",
             SR_NAME_FNAME     = "superReadNames.txt"; 

int BUFFER_SIZE = 1000;
int NUM_THREADS = 8;


bool build_index(std::string sr_name_fname, std::string sr_seq_fname, int max_ku);

int main(int argc, char **argv) {
    
    //Keeps track of total time for various tasks
    //DEBUG_TIMER shoud be set in "superread_index.hpp"
    #ifdef DEBUG_TIMER
    double index_time = 0,
           input_time = 0,
           submit_time = 0,
           output_time = 0,
           total_thread_time = 0;
    Timer timer;
    int thread_count = 0;
    #endif

    //Stores option arguments
    std::string work_dir = DEF_WORK_DIR, 
           read_ku_path,
           sr_fasta_path, 
           sr_name_path,
           //pe1_fname,
           //pe2_fname,
           //unpaired_fname,
           out_prefix;
    
    int max_ku = 0; 
    
    //Parse arguments
    char o;
    while ( (o = getopt(argc, argv, OPTIONS)) != -1 ) {
        switch (o) {
            case WORK_DIR_OPT:
                work_dir = std::string(optarg);
                if (work_dir[work_dir.size()-1] != '/')
                    work_dir += "/";
                break;
            case READ_KU_OPT:
                read_ku_path = std::string(optarg);
                break;
            case OUT_PREFIX_OPT:
                out_prefix = std::string(optarg);
                break;
            case NUM_THREADS_OPT:
                NUM_THREADS = atoi(optarg);
                break;
            case BUFFER_SIZE_OPT:
                BUFFER_SIZE = atoi(optarg);
                break;
            case KU_NUM_OPT:
                max_ku = atoi(optarg);
                break;
        }
    }

    //Set default filenames
    sr_fasta_path = work_dir + SR_FASTA_FNAME;
    sr_name_path = work_dir + SR_NAME_FNAME;


    if (read_ku_path.empty())
        read_ku_path = work_dir + DEF_READ_KU_FNAME;

    std::vector<std::string> short_read_fnames;
    for (int i = optind; i < argc; i++) {
        short_read_fnames.push_back(argv[i]);
    }
    
    ReadIO short_read_io(read_ku_path, short_read_fnames);

    //Build index    
    std::cerr << "Building super-read index...\n";
    if (!SuperreadIndex::load_files(sr_name_path, sr_fasta_path, max_ku))
        return 1;

    #ifdef DEBUG_TIMER
    index_time = timer.lap();
    #endif

    //Init thread pool
    int thread_ret; //Unused thread return value
    thread_pool<IndexQuery*, int> threads(NUM_THREADS, 
                                  SuperreadIndex::match_reads);

    //Init thread input queries
    std::vector<IndexQuery*> queries(NUM_THREADS);
    for (int q = 0; q < NUM_THREADS; q++) {
        queries[q] = new IndexQuery(BUFFER_SIZE, true);
    }
    
    //Query buffer
    IndexQuery *query_buf = new IndexQuery(BUFFER_SIZE, true);

    //true if buffer is ready to be filled
    //false if buffer points to thread parameter (not yet swapped)
    bool buffer_swapped = true;

    //Stores each mate pair
    //Read *read1 = NULL, *read2 = NULL;
    
    //Stores IDs of matched reads (only mate 1, mate 2 assumed)
    //std::vector<int> matched_reads;

    //True when all reads have been read and all threads are finished
    bool all_reads_assigned = false;

    std::cerr << "Index built, assigning reads...\n";

    //ifstream read_ku_file(read_ku_path.c_str());
    //std::string read_ku_line;

    //if (!read_ku_file.is_open()) {
    //    std::cerr << "Error: unable to open '" << read_ku_path << "'\nExiting\n";
    //    return 1;
    //}

    //Read K Unitig file line-by-line
    while(!all_reads_assigned) {
    
        #ifdef DEBUG_TIMER
	    timer.reset();
	    #endif

        all_reads_assigned = !short_read_io.fill_buffer(query_buf);
        //all_reads_assigned = !query_buf->fill_buffer(read_ku_file, read1, read2);
        
        #ifdef DEBUG_TIMER
        input_time += timer.lap();
        #endif

    	//Submit job
        if (!all_reads_assigned) {
            threads.submit_job(&query_buf, &thread_ret);
            buffer_swapped = false; 
        }
        
        #ifdef DEBUG_TIMER
        submit_time += timer.lap();
        #endif

        //Check status of running threads
        for (int q = 0; q < NUM_THREADS; q++) {
            if (queries[q]->finished) {
                #ifdef DEBUG_TIMER
        		total_thread_time += queries[q]->thread_time;
                #endif

                //Output results if present
                if (queries[q]->results->size() > 0) {
                    #ifdef DEBUG_TIMER
                    thread_count++;
                    #endif

                    SuperreadIndex::update_frag_counts(queries[q]);

                    short_read_io.update_matched(queries[q]);

                    queries[q]->clear();
                }

                //Store buffer results if needed
                if (!buffer_swapped) {
                    std::swap(query_buf, queries[q]);
                    query_buf->finished = false;
                    buffer_swapped = true;
                }
            
            //Threads still running
            } else {
                all_reads_assigned = false;
            }
        }
        
        #ifdef DEBUG_TIMER
        output_time += timer.get();
        #endif
    }
    
    std::cerr << "All reads assigned, calculating weights..." << "\n";
    //std::vector< std::vector<double> > all_weights;
    //all_weights.push_back(SuperreadIndex::weights0);


    
    for (int i = 0; i < 100; i++) {
        if (SuperreadIndex::iter_weights() < 0.01)
            break;
    }
    
    //for (int i = 0; i < SuperreadIndex::weights.size(); i++) {
    //    cout<< SuperreadIndex::sr_ids[i] << " " << SuperreadIndex::weights[i] << " ";
    //}

    std::cerr << "Weights calculated, writing superread fastq\n";
    SuperreadIndex::write_fastq(out_prefix + "superreads.fastq");


    std::cerr << "Writing unmatched read fastq\n";
    short_read_io.write_unmatched(out_prefix + "pe1.fastq",
                                  out_prefix + "pe2.fastq",
                                  out_prefix + "unpaired.fastq");

    std::cerr << "Cleaning up" << "\n";
    
    threads.release_workers();
    SuperreadIndex::destroy();

    //delete short_read_io;
    delete query_buf;

    for (int q = 0; q < NUM_THREADS; q++) {
        delete queries[q];
    }
    
    #ifdef DEBUG_TIMER
    std::cerr << "\nMain thread\n";
    std::cerr << "Index time:  " << index_time  << "s\n"; 
    std::cerr << "Parse time:  " << input_time  << "s\n";
    std::cerr << "Wait time:   " << submit_time << "s\n";
    std::cerr << "Output time: " << output_time << "s\n";
    
    std::cerr << "\nIndex threads\n";
    std::cerr << "Total thread time: " << total_thread_time <<  "\n";
    std::cerr << "Threads run:       " << thread_count << "\n";
    std::cerr << "Time per thread:   " << (total_thread_time / thread_count) << "\n";
    #endif

    return 0;
}


