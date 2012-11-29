#include<ftw.h>
#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/wait.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <pthread.h>
#include <err.hpp>
#include <misc.hpp>
#include <src2/collectReadSequencesForLocalGapClosing_cmdline.hpp>

typedef std::string stdString;
typedef std::vector<stdString> vectorOfStrings;
typedef std::set<stdString> setOfStrings;
typedef std::unordered_map<stdString, stdString> readNameToReadSequence;
typedef std::unordered_map<stdString, int> stringToIntMap;
typedef std::unordered_map<stdString, stdString> stringToStringMap;

struct arguments {
     int dirNum;
     charb fauxReadFileDataStr;
//     basic_charb<remaper<char> > readFileDataStr;
     std::vector<vectorOfStrings> *readsInGroup;
     std::vector<vectorOfStrings> *mateReadsInGroup;
     readNameToReadSequence *readSeq; // Read name to read sequence     
};

int remove_function(const char *fpath, const struct stat *sb,
                          int typeflag, struct FTW *ftwbuf){
return(remove(fpath));
}

FILE *fopen_wait(int timeout_, char *filename) {
     int timeout=0;
     FILE *infile;
     do {
	  infile = fopen (filename, "r");
	  if(infile != NULL) break;
	  sleep(1);
	  timeout++;
     }
     while(timeout < timeout_);
     if(infile == NULL)
          printf("Timed out on file %s\n", filename);
     return(infile);
}

int reportNumGaps (const char *contigEndSeqFile);
int analyzeGap(struct arguments threadArg); // Returns int for now

cmdline_parse args;

std::string exeDir;

// Names of files and directories

struct numAndOriStruct {
     int groupNum;
     int ori;
};

stdString getReadMateName (stdString readName);
void loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq);
void checkArgs (void);


int main(int argc, char **argv)
{
     FILE *infile; // , *outfile;
     FILE *contigEndSeqFile;
     struct arguments threadArgs;
     charb line(100);
     vectorOfStrings fauxReadGroups;
     stringToIntMap group;
     std::unordered_map<stdString, char> end;
     setOfStrings readIsNeeded;
     readNameToReadSequence readSeq; // Read name to read sequence
     struct stat     statbuf;
     charb tempBuffer(100);
     args.parse (argc, argv);
     checkArgs ();
     //thread_pool<struct arguments, int> pool (args.num_threads_arg, analyzeGap);
     char *tempPtr = strrchr (argv[0], '/');
     if (tempPtr == NULL)
	  exeDir = std::string(".");
     else {
	  unsigned int diff = tempPtr - argv[0];
	  strcpy (tempBuffer, argv[0]);
	  tempBuffer[diff] = 0;
	  exeDir = std::string (tempBuffer); }

     contigEndSeqFile = fopen (args.contig_end_sequence_file_arg, "r");
     if (contigEndSeqFile == NULL) {
	  fprintf (stderr, "File '%s' doesn't exist! Bye!\n", args.contig_end_sequence_file_arg);
	  exit (1); }
     fclose (contigEndSeqFile);
     int numGaps = reportNumGaps (args.contig_end_sequence_file_arg);

     infile = fopen(args.faux_reads_file_arg, "r");
     while (fgets (line, 100, infile)) {
	  fauxReadGroups.push_back(stdString(line));
	  char *cptr = line+1;
	  char *fauxReadName;
	  char *saveptr;
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  stdString fauxReadNameStr = stdString(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'F';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
	  cptr = line+1;  // Don't know if it was reset (what happens in scanf)
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  fauxReadNameStr = stdString(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'R';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
     }
     fclose (infile);

     int numFauxReadGroups = (int) fauxReadGroups.size();
     
     infile = fopen (args.faux_read_matches_to_kunis_file_arg, "r");
     stdString fauxRead;
     ExpBuffer<char *> flds;
     char matchToKeepMate;

     typedef std::vector<numAndOriStruct> numAndOriList;
     numAndOriList emptyNumAndOriList;
     std::vector<numAndOriList> kUniMatchStructs;
     std::vector<vectorOfStrings> readsInGroup, mateReadsInGroup;
     vectorOfStrings emptyStringVector;
     for (int i=0; i<numFauxReadGroups; i++) {
	  readsInGroup.push_back (emptyStringVector);
	  mateReadsInGroup.push_back (emptyStringVector);
     }

     int maxKUniSeen = -1;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  // We will keep the info on what kind of a match we need to bring
	  // in a mate as well as a read
	  fauxRead = stdString (flds[0]);
	  int fauxReadGroup = group[fauxRead];
	  char fauxReadGapEnd = end[fauxRead];
	  for (unsigned int i=2; i<flds.size(); i+=3) {
	       int kUni = atoi (flds[i]);
	       char relOri = *flds[i+2];
	       while (kUni > maxKUniSeen) {
		    ++maxKUniSeen;
		    kUniMatchStructs.push_back(emptyNumAndOriList); }
	       if (fauxReadGapEnd == relOri)
		    matchToKeepMate = 'F';
	       else
		    matchToKeepMate = 'R';
	       numAndOriStruct numAndOri;
	       numAndOri.groupNum = fauxReadGroup;
	       numAndOri.ori = matchToKeepMate;
	       kUniMatchStructs[kUni].push_back(numAndOri);
	  }
     }
     fclose (infile);

     infile = fopen (args.read_matches_to_kunis_file_arg, "r");
     // Ex. line: pe2836 101 610 10 F
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  stdString readName = stdString (flds[0]);
	  int flds0len = strlen(flds[0]);
	  if (flds[0][flds0len-1] % 2 == 0)
	       ++flds[0][flds0len-1];
	  else
	       --flds[0][flds0len-1];
	  stdString readMate = stdString (flds[0]);
	  int kUni = atoi (flds[2]);
	  char relOri = *(flds[4]);
	  for (unsigned int i=0; i<kUniMatchStructs[kUni].size(); i++) {
	       readsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readName);
	       readIsNeeded.insert (readName);
	       if (kUniMatchStructs[kUni][i].ori == relOri) {
		    readIsNeeded.insert (readMate);
		    mateReadsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readMate); }
	  }
     }

     if (stat (args.dir_for_gaps_arg, &statbuf) != 0)
          mkdir((char *)args.dir_for_gaps_arg, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
     contigEndSeqFile = fopen (args.contig_end_sequence_file_arg, "r");

     // Here we make the output directory
     if (stat (args.output_dir_arg, &statbuf) != 0)
          mkdir((char *)args.output_dir_arg, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//     int outputGroupNum = 0;

     charb tempLine(1000);
     kUniMatchStructs.clear();    
     pid_t pid;
     int status;
     loadNeededReads (readIsNeeded, readSeq);

     for (unsigned int dirNum=0; dirNum<readsInGroup.size(); dirNum++) {
     if(dirNum>=args.num_threads_arg)
	     wait(&status);
     //as soon as wait releases, submit another job
      
          threadArgs.fauxReadFileDataStr.clear();
          for (int j=0; j<4; j++) {
               fgets (tempLine, 1000, contigEndSeqFile);
               strcat (threadArgs.fauxReadFileDataStr, tempLine); }
          threadArgs.readsInGroup = &readsInGroup;
          threadArgs.mateReadsInGroup = &mateReadsInGroup;
          threadArgs.readSeq = &readSeq;
          threadArgs.dirNum = dirNum;
          fflush(stdout);
          pid=fork();
          if(pid<0)
                perror("Fork failed");
          else if(pid==0){
                	analyzeGap(threadArgs);
                	exit(0);
                }
	}

     for(unsigned int ttt=0;ttt<readsInGroup.size();ttt++){
	wait(&status);
	}   

        if (! args.keep_directories_flag)
            nftw((char *) args.output_dir_arg,remove_function,2048,FTW_DEPTH|FTW_PHYS);
	return 0;
}

void loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq)
{
     FILE *infile;
     static unsigned int currentFileNum = 0;
//     static off_t currentFileOffset = (off_t) 0;
     const char *localReadsFile;
     uint64_t numReadSeqsLoaded;
     
     readSeq.clear();
     numReadSeqsLoaded = 0;
     while (currentFileNum < args.reads_file_arg.size()) {
	  localReadsFile = args.reads_file_arg[currentFileNum];
	  infile = fopen (localReadsFile, "r");
//	  fseeko (infile, currentFileOffset, SEEK_SET);
	  int offsetStart = strlen (localReadsFile)-5;
	  charb line(100), readNameStr(100);
	  char *cptr = line+1;
	  stdString readName;
	  bool isNeeded = false;
	  if (offsetStart < 0) offsetStart = 0;
	  const char *cptr2 = localReadsFile + offsetStart;
	  if (strcmp (cptr2, "fastq") == 0)
	       goto fastQ;
	  // If we get here it's a fasta file
	  while (fgets (line, 100, infile)) {
	       cptr = line+1;
	       if (line[0] == '>') {
		    char *saveptr;
		    readNameStr = strtok_r(cptr, " \t\n", &saveptr);
		    readName = stdString (readNameStr);
		    if (readIsNeeded.find(readName) != readIsNeeded.end())
			 isNeeded = true;
		    else
			 isNeeded = false;
		    continue; }
	       if (isNeeded) {
		    line[strlen(line)-1] = 0; // chop (line);
		    if (readSeq.find(readName) == readSeq.end()) {
			 readSeq[readName] = stdString ("");
			 ++numReadSeqsLoaded; }
		    readSeq[readName] += stdString (line);
	       }
	  }
	  goto endProcessingReads;
     fastQ:
	  while (fgets (line, 100, infile)) {
	       cptr = line+1;
	       char *saveptr;
	       isNeeded = false;
	       readNameStr = strtok_r(cptr, " \t\n", &saveptr);
	       readName = stdString (readNameStr);
	       if (readIsNeeded.find(readName) != readIsNeeded.end())
		    isNeeded = true;
	       else
		    isNeeded = false;

	       fgets (line, 100, infile);
	       line[strlen(line)-1] = 0; // chop (line);
	       if (isNeeded) {
		    readSeq[readName] = stdString (line);
		    ++numReadSeqsLoaded;
	       }
	       fgets (line, 100, infile); // Skip quality bases
	       fgets (line, 100, infile);
	  }
     endProcessingReads:
	  fclose (infile);
	  ++currentFileNum;
//	  currentFileOffset = (off_t) 0;
     }
     return;
}

stdString getReadMateName (stdString readName)
{
     stdString readMate;
     int readNameLen = readName.length();
     charb readNameStr(100);
     if (readNameLen == 0)
	  return (stdString (""));
     strcpy (readNameStr, readName.c_str());
     if (readNameStr[readNameLen-1] % 2 == 0)
	       ++readNameStr[readNameLen-1];
	  else
	       --readNameStr[readNameLen-1];
     readMate = stdString (readNameStr);

     return (readMate);
}

void checkArgs (void)
{
     struct stat     statbuf;
     int fail = 0;
     if (stat (args.faux_reads_file_arg, &statbuf) != 0) {
	  std::cerr << "Faux reads file '" << args.faux_reads_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (args.faux_read_matches_to_kunis_file_arg, &statbuf) != 0) {
	  std::cerr << "File of faux read matches to k-unitigs '" << args.faux_read_matches_to_kunis_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (args.read_matches_to_kunis_file_arg, &statbuf) != 0) {
	  std::cerr << "File of read matches to k-unitigs '" << args.read_matches_to_kunis_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     for (unsigned int i=0; i<args.reads_file_arg.size(); i++) {
	  const char *readFile = args.reads_file_arg[i];
	  if (stat (readFile, &statbuf) != 0) {
	       std::cerr << "Read file '" << readFile << "' doesn't exist.\n";
	       fail = 1; }
     }
     
     if (fail)
	  exit (1);
}

int analyzeGap(struct arguments threadArg)
{
     struct stat statbuf;
     FILE *infile, * outfile;
     std::vector<vectorOfStrings> *readsInGroup = threadArg.readsInGroup;
     std::vector<vectorOfStrings> *mateReadsInGroup = threadArg.mateReadsInGroup;
     readNameToReadSequence *readSeq = threadArg.readSeq; // Read name to read sequence     

     // Doing the actual work of the worker thread
     charb outDirName(100), tempFileName(10), cmd(100), line(100);
     unsigned int dirNum = threadArg.dirNum;

     sprintf (outDirName, "%s/gap%09ddir", args.output_dir_arg, dirNum);
    
     if (stat (outDirName, &statbuf) != 0) 
	  mkdir((char *)outDirName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

     sprintf (tempFileName, "%s/fauxReads.fasta", (char *) outDirName);
     outfile = fopen (tempFileName, "w");
     fputs (threadArg.fauxReadFileDataStr, outfile);
     fclose (outfile);
     sprintf (tempFileName, "%s/reads.fasta", (char *) outDirName);
     outfile = fopen (tempFileName, "w");

     stringToIntMap willBeOutput;
     willBeOutput.clear();
     vectorOfStrings readsToOutput;
     readsToOutput.clear();
     for (unsigned int readNum=0; readNum<(*readsInGroup)[dirNum].size(); readNum++) {
	  stdString readName = (*readsInGroup)[dirNum][readNum];
	  stringToIntMap::iterator it = willBeOutput.find (readName);
	  if (it != willBeOutput.end())
	       continue;
	  readsToOutput.push_back (readName);
	  willBeOutput[readName] = 2; }
     for (unsigned int readNum=0; readNum<(*mateReadsInGroup)[dirNum].size(); readNum++) {
	  stdString readName = (*mateReadsInGroup)[dirNum][readNum];
	  stringToIntMap::iterator it = willBeOutput.find (readName);
	  if (it != willBeOutput.end())
	       continue;
	  willBeOutput[readName] = 1; }

     // Variables used only if we keep the directories; put here for the compiler
     int outCount = 0;
     setOfStrings alreadyOutput;
     vectorOfStrings bothMatesHaveMatches, mateBroughtInViaMatePair, readOnlys;

     if (! args.keep_directories_flag) {
	  fputs (">Roeroeroeyourboat\n", outfile);
	  stringToIntMap::iterator it;
	  for (it = willBeOutput.begin(); it != willBeOutput.end(); it++) {
	       stdString readName = (*it).first;
	       stringToStringMap::iterator readSeqIt =
		    readSeq->find (readName);
	       if (readSeqIt == readSeq->end())
		    continue;
	       fputs ((*readSeq)[readName].c_str(), outfile);
	       fputc ('N', outfile); }
	  fputc ('\n', outfile);
          willBeOutput.clear();
          readsToOutput.clear();
	  goto closeTheReadsFileAndGoOn; }

     // If we get here we are keeping the directories and doing a better reads file
     // We used the goto above to keep this large section from being overly indented
     alreadyOutput.clear();
     bothMatesHaveMatches.clear();
     mateBroughtInViaMatePair.clear();
     readOnlys.clear();
     for (unsigned int readNum=0; readNum<readsToOutput.size(); readNum++) {
	  stdString readName = readsToOutput[readNum];
	  if (alreadyOutput.find(readName) != alreadyOutput.end())
	       continue;
	  stdString readMate = getReadMateName (readName);
	  alreadyOutput.insert (readName);
	  stringToIntMap::iterator it2 = willBeOutput.find (readMate);
	  if (it2 != willBeOutput.end()) {
	       alreadyOutput.insert(readMate);
	       if (it2->second == 2) {
		    bothMatesHaveMatches.push_back (readName);
		    bothMatesHaveMatches.push_back (readMate); }
	       else {
		    mateBroughtInViaMatePair.push_back (readName);
		    mateBroughtInViaMatePair.push_back (readMate); }
	  }
	  else
	       readOnlys.push_back (readName);
     }
     for (unsigned int readNum=0; readNum<bothMatesHaveMatches.size(); readNum++) {
	  stringToStringMap::iterator readSeqIt =
	       readSeq->find (bothMatesHaveMatches[readNum]);
	  if (readSeqIt == readSeq->end())
	       continue;
	  fprintf (outfile, ">%s B\n%s\n", bothMatesHaveMatches[readNum].c_str(), (readSeqIt->second).c_str());
     }
     
     for (unsigned int readNum=0; readNum<mateBroughtInViaMatePair.size(); readNum++) {
	  ++outCount;
	  stdString extraStr;
	  stdString readName = mateBroughtInViaMatePair[readNum];
	  stringToStringMap::iterator readSeqIt =
	       readSeq->find (readName);
	  if (readSeqIt == readSeq->end())
	       continue;
	  if (outCount % 2 == 1)
	       extraStr = stdString ("match");
	  else
	       extraStr = stdString ("mate");
	  fprintf (outfile, ">%s %s M\n%s\n", readName.c_str(), extraStr.c_str(), ((*readSeq)[readName]).c_str());
     }
     
     for (unsigned int readNum=0; readNum<readOnlys.size(); readNum++) {
	  stdString readName = readOnlys[readNum];
	  stringToStringMap::iterator readSeqIt =
	       readSeq->find (readName);
	  if (readSeqIt == readSeq->end())
	       continue;
	  stdString mateRead = getReadMateName (readName);
	  fprintf (outfile, ">%s O\n%s\n>%s N\nN\n", readName.c_str(), ((*readSeq)[readName]).c_str(), mateRead.c_str());
     }
     
closeTheReadsFileAndGoOn:
     fclose (outfile);
     sprintf (cmd, "%s/closeGaps.oneDirectory.perl --dir-to-change-to %s --Celera-terminator-directory %s --reads-file reads.fasta --output-directory outputDir --max-kmer-len %d --min-kmer-len %d --maxnodes %d --mean-for-faux-inserts %d --stdev-for-faux-inserts %d --use-all-kunitigs --noclean 1>%s/out.err 2>&1", exeDir.c_str(), (char *) outDirName, args.Celera_terminator_directory_arg, args.max_kmer_len_arg, args.min_kmer_len_arg, args.max_nodes_arg, args.mean_for_faux_inserts_arg, args.stdev_for_faux_inserts_arg, (char *) outDirName, (char *) outDirName);
     sprintf (tempFileName, "%s/passingKMer.txt", (char *) outDirName);
     int passingKMerValue = 0;

     system (cmd);
     /* Now, if "passingKMer.txt" exists in outDirName, copy the files
	superReadSequences.fasta and
	readPlacementsInSuperReads.final.read.superRead.offset.ori.txt (after appropriate
	modifications), copy to the desired output directory from (e.g.)
	'outDirName'/work_localReadsFile_41_2 */

     while(passingKMerValue <11){
     	infile = fopen_wait (60,tempFileName);
	if (infile == NULL) {
	     passingKMerValue = 11;
	     break; }
	fgets(line,100,infile);
	passingKMerValue=atoi(line);
     	fclose (infile);
     }

     if(passingKMerValue == 11){
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
        }
     charb superReadFastaString(1000);     
     // If we get here we have found a join and passingKMer.txt exists
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/superReadSequences.fasta", (char *) outDirName, passingKMerValue);
     infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
	  fprintf(stderr,"failed to open superReadSequences.fasta %s passingKMerValue %d\n", (char *) outDirName,passingKMerValue);
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
        } 
            
     fgets (line, 100, infile);
     sprintf (superReadFastaString, "%d ", threadArg.dirNum);
     while (fgets (line, 100, infile)){
	  line.chomp();
	  strcat (superReadFastaString, line);
	}
//     while (fgets_append (*(threadArg.superReadFastaString), infile))
//          ;
     fclose (infile);
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt", (char *) outDirName, passingKMerValue);
        infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
          fprintf(stderr,"failed to open readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s passingKMerValue %d\n", (char *) outDirName,passingKMerValue);
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
	}
     ExpBuffer<char *>flds;
     charb readPlacementLine,readPlacementLines;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
          sprintf (readPlacementLine, " %s %d %s %s", flds[0], threadArg.dirNum, flds[2], flds[3]);
          strcat(readPlacementLines,readPlacementLine);
     }
     fclose (infile);
     fprintf(stderr,"%s%s\n",(char*)superReadFastaString,(char*)readPlacementLines);
#if 0
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
#endif

//printf("Gap %s passing kmer %d super read %s\n",(char *) outDirName,passingKMerValue,(char*)*(threadArg.superReadFastaString));
     if (! args.keep_directories_flag)
            nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
     
     return (0);
}

int reportNumGaps (const char *fn)
{
     charb cmd(100), line(100);
     FILE *infile;

     sprintf (cmd, "tail -2 %s | head -1", fn);
     infile = popen (cmd, "r");
     fgets (line, 100, infile);
     char *cptr = line;
     while (! isdigit(*cptr))
          cptr++;
     int lastFauxContig = atoi (cptr);
     int numGaps = (lastFauxContig+1)/2;
     return (numGaps);
}

