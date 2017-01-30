#ifndef SUPERREAD_INDEX_H
#define SUPERREAD_INDEX_H

#include <string.h>
#include <vector>
#include <map>

//#define DEBUG_TIMER 1

#define PREFIX_LEN 2

class Read {
    public:
    Read() {}
    Read(const Read &r); 

    Read(std::string read_str);
    bool is_mate(Read r);
    bool same_prefix(Read r);
    
    bool is_valid;
    int read_id, length;
    char prefix[PREFIX_LEN];
    std::vector<int> ku_ids;
    std::vector<char> ku_strands;
};

class IndexQuery {
    public:
    IndexQuery(int buffer_size, bool paired);
    ~IndexQuery();

    //bool fill_buffer(std::ifstream &infile, Read *read1=NULL, Read *read2=NULL);
    bool add_read(Read *read1, Read *read2=NULL);
    bool is_full();
    bool is_empty();
    void clear();
	
    std::vector<Read*> reads;
    std::vector< std::vector<int>* > *results;
    std::vector<Read*> *matched_reads;
    bool finished, is_paired;
    int prefix_id;

    #ifdef DEBUG_TIMER
    float thread_time;
    #endif
};

class ReadIO {
    public:
    // [MEAN:[STDEV:]]FRAG1[:FRAG2]
    ReadIO(std::string read_ku_fname, std::vector<std::string> read_files);
    ~ReadIO();

    bool fill_buffer(IndexQuery *q);

    void update_matched(IndexQuery *q);

    void write_unmatched(std::string pe1_fname, 
                         std::string pe2_fname, 
                         std::string unpaired_fname);
    
    private:
    int get_prefix_id(Read r);
    bool copy_fastq(std::ifstream &in, std::ofstream &out);
    bool skip_fastq(std::ifstream &in);

    std::ifstream read_ku_file;
    std::vector<std::string> fasta_fnames;
    std::vector< std::vector<int> > matched_reads;
    std::vector<bool> is_paired;
    Read *read1, *read2;
    std::string prefixes;
};

class SuperreadIndex {
    public:

    //Creates a new index given a list of superread names (K-unitig sequences) 
    //and a fasta of reduced SRs (only SRs present in fasta will be included).
    //Will reserve space for 'max_ku' K-unitigs if included
    //Returns true if successful, false otherwise
    static bool load_files(std::string sr_name_path, 
                           std::string sr_fasta_path, 
                           int max_ku=0);

    //Frees all memory used by index
    static void destroy();

    //Matches reads to superreads and stores results in 'q'
    //Designed to be called as a thread
    //Always returns 0 (required for 'thread_pool')
    static int match_reads(IndexQuery *q);
    
    //Increments fragment (reads or mates) counts for matched superread groups
    static void update_frag_counts(IndexQuery *q);
    
    //Perform iteration of estimation/maximization to update superread weights
    //Weights initially set to number of reads/fragments uniquely assigned
    //Returns maximum absolute difference of weight estimates for each superread
    static double iter_weights();

    //Writes superreads to fastq with superread weights
    //Should be called after iter_weights()
    static void write_fastq(std::string fname, std::string delim=":");

    private:

    //Stores a record of a KU in a superread
    typedef struct kunitig {
        int sr,     //Which superread it is in
            sr_pos; //Location in superread, starting at 0 from the leftmost KU
        char strand;//Forward (F) or reverse (R) strand
    } *KUnitig;
   
    static std::string sr_fasta_fname;
    
    //Stores locations of all K-unitigs
    //Index in list is equivilent to KU number
    static std::vector< std::vector<KUnitig> > ku_list;
    
    //Stores superread IDs corrasponding to other SR lists
    static std::vector<int> sr_ids;

    //Stores lengths of superreads
    static std::vector<int> sr_lens;

    //Stores counts of reads/fragments assigned to sets of superreads
    //Used to calculate weights
    static std::map< std::vector<int>, int > sr_frag_counts;

    //Stores current suprread weight estimates
    //Before iter_weights() is called, stores number of reads/fragments uniquely
    //assigned to each superread (initial estimate)
    static std::vector<double> weights;

    //Parses superread std::string and indexes locations K-unitigs
    //'sr_index' should corraspond to element in 'sr_ids' and 'sr_lens'
    static void add_superread_kus(int sr_index, std::string sr_str);

    //Helper functions for 'match_reads'
    static void match_paired_reads(IndexQuery *q);
    static void match_unpaired_reads(IndexQuery *q);

    //Returns all superreads that contain given read
    static std::vector<int> match_read(Read *read);

    //Returns denominator for updated superread weight
    //based on previous weight and superread length
    static inline double sr_weight_share(int sr);
};


#endif


