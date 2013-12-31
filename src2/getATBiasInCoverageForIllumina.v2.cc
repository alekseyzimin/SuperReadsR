/* Take in a delta-type file (as generated by sam2_to_delta) and
   report coverage of bases up to 255 (1 unsigned char per base)
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <charb.hpp>
#include <misc.hpp>
#include <exp_vector.hpp>

typedef std::pair<int, int> intPair;

FILE *Fopen (const char *fn, const char *mode);

// Change this to half the length of the interval of interest
const int LENGTH_ON_EACH_SIDE_OF_REGION = 100;
const int MIN_SEQUENCE_TO_CONTINUE = 40000;
const int FASTQ = 1;
const int FASTA = 2;
const int LENGTHS_FILE = 3;

int main (int argc, char **argv)
{
     charb consensusSequence (100000000), outputLine (10000), line(1000);
     // Process the args ...
     char seqFilename[200], placementFile[200], readLengthsFile[200];
     int minOffsetFromEndToAccept = 0;
     int defaultReadLength = 101;
     int minimumReadLength = 10;
     std::vector<char *> flds;
     strcpy (readLengthsFile, "");
     for (int i=1; i<argc; ++i) { // Processing the args here
	  if (strcmp (argv[i], "--sequence-file") == 0) {
	       ++i;
	       strcpy (seqFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "--placement-file") == 0) {
	       ++i;
	       strcpy (placementFile, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "--min-from-end") == 0) {
	       ++i;
	       minOffsetFromEndToAccept = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "--read-lengths-file") == 0) {
	       ++i;
	       strcpy (readLengthsFile, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "--default-read-length") == 0) {
	       ++i;
	       defaultReadLength = atoi (argv[i]);
	       continue; }
     }
     // Done processing the args
     if (minOffsetFromEndToAccept < LENGTH_ON_EACH_SIDE_OF_REGION)
	  minOffsetFromEndToAccept = LENGTH_ON_EACH_SIDE_OF_REGION;
     // Processing the consensus sequence ...
     FILE *infile = Fopen (seqFilename, "r");
     char *currentSequenceStart = consensusSequence;
     char *cptr = consensusSequence;
     std::string seqNameHold;
     std::vector<std::string> seqNames;
     std::map<std::string, int> seqLen;
     std::map<std::string, char *> sequenceStart;
     charb localConsensusSequence(1000);
     unsigned long long totalLengthOfUsableSequence = 0;
     int minSequenceLength = 2 * minOffsetFromEndToAccept;
     bool isFirstLine = true;
     while (fgets (line, 100, infile)) {
	  if (line[0] == '>') {
	       if (! isFirstLine) {
		    int tempLen = localConsensusSequence.len();
		    if (tempLen > minSequenceLength) {
			 seqLen[seqNameHold] = localConsensusSequence.len();
			 sequenceStart[seqNameHold] = currentSequenceStart;
			 seqNames.push_back (seqNameHold);
			 totalLengthOfUsableSequence += (seqLen[seqNameHold] - minSequenceLength);
			 strcat (consensusSequence, localConsensusSequence);
			 cptr = consensusSequence + consensusSequence.len();
		    }
	       }
	       else
		    isFirstLine = false;
	       localConsensusSequence.clear();
	       getFldsFromLine (line, flds);
	       seqNameHold = std::string(flds[0]+1);
	       currentSequenceStart = cptr;
	       continue;
	  }
	  line.chomp();
	  strcat (localConsensusSequence, line);
     }
     fclose (infile);
     if (currentSequenceStart != consensusSequence) {
	  int tempLen = localConsensusSequence.len();
	  if (tempLen > minSequenceLength) {
	       seqLen[seqNameHold] = localConsensusSequence.len();
	       sequenceStart[seqNameHold] = currentSequenceStart;
	       seqNames.push_back (seqNameHold);
	       totalLengthOfUsableSequence += (seqLen[seqNameHold] - minSequenceLength);
	       strcat (consensusSequence, localConsensusSequence);
	       cptr = consensusSequence + consensusSequence.len();
	  }
     }
     unsigned long long totalSequenceLength = cptr - consensusSequence;
     // Done processing the consensus sequence
     // Check we have enough usable sequence
     if (totalLengthOfUsableSequence < MIN_SEQUENCE_TO_CONTINUE) {
	  fprintf (stderr, "Need to have at least %d bases of sequence at least %d from the ends. Only have %llu. Bye!\n", MIN_SEQUENCE_TO_CONTINUE, minOffsetFromEndToAccept, totalLengthOfUsableSequence);
	  exit (-1); }
     
     // Get the read lengths from the file if specified
     std::map<std::string, int> readLengths;
     if (strlen (readLengthsFile) != 0) { // We're reading in the read lengths file
	  infile = Fopen (readLengthsFile, "r");
	  fgets (line, 100, infile);
	  int fileType = 0;
	  if (line[0] == '@')
	       fileType = FASTQ;
	  else if (line[0] == '>')
	       fileType = FASTA;
	  else {
	       int numFlds = getFldsFromLine (line, flds);
	       if (numFlds == 2) {
		    fileType = LENGTHS_FILE;
		    for (char *cptr2=flds[1]; *cptr2; ++cptr2)
			 if (! isdigit (*cptr2))
			      fileType = 0;
	       }
	  }
	  if (! fileType) {
	       fprintf (stderr, "Unknown file type for the lengths file '%s'. Bye!\n", readLengthsFile);
	       exit (-1); }
	  rewind (infile);
	  if (fileType == FASTQ) {
	       while (fgets (line, 100, infile)) {
		    getFldsFromLine (line, flds);
		    std::string name = std::string (flds[0]+1);
		    fgets (line, 100, infile);
		    int tlen = strlen(line)-1;
		    if (tlen >= minimumReadLength) // Only do reads longer than a minimum
			 readLengths[name] = strlen (line)-1;
		    fgets (line, 100, infile);
		    fgets (line, 100, infile);
	       }
	  }
	  else if (fileType == FASTA) {
	       std::string name;
	       while (fgets (line, 100, infile)) {
		    if (line[0] == '>') {
			 getFldsFromLine (line, flds);
			 name = std::string (flds[0]+1);
			 readLengths[name] = 0;
			 continue; }
		    readLengths[name] += (strlen (line) - 1);
	       }
	  }
	  else { // It's a read lengths file
	       while (fgets (line, 100, infile)) {
		    getFldsFromLine (line, flds);
		    int tlen = atoi (flds[1]);
		    if (tlen > minimumReadLength)
			 readLengths[std::string(flds[0])] = atoi (flds[1]); }
	  }
	  fclose (infile);
     }
	  
//     char *seqFilename = "/genome2/raid/tri/rhodobacter/finished_sequence/Rsphaeroides.1con.fa";
//     char *filename = "/genome2/raid/tri/rhodobacter/pe.renamed_genome.finishedSeqMatches_aln.k5.fromSamFile.bestAlignsLe6errs.delta";
     std::vector<int> beginAndEndNetCounts (totalSequenceLength,0);
     std::vector<unsigned long long> Acounts (totalSequenceLength,0), Ccounts (totalSequenceLength,0), Gcounts (totalSequenceLength,0), Tcounts (totalSequenceLength,0), invalidCounts (totalSequenceLength,0);
     std::map<intPair, int> countsForPctCoveragePairs;
     std::vector<int> &coverageCounts = beginAndEndNetCounts;
     Acounts[0] = Ccounts[0] = Gcounts[0] = Tcounts[0] = invalidCounts[0] = 0;
     printf ("totalSequenceLength = %llu\n", totalSequenceLength);
     printf ("consensusSequence length = %d\n", (int) strlen (consensusSequence));
     for (unsigned long long i=0; i<totalSequenceLength; ++i) {
	  Acounts[i+1] = Acounts[i];
	  Ccounts[i+1] = Ccounts[i];
	  Gcounts[i+1] = Gcounts[i];
	  Tcounts[i+1] = Tcounts[i];
	  invalidCounts[i+1] = invalidCounts[i];
	  switch (consensusSequence[i]) {
	  case 'A': case 'a': ++Acounts[i+1]; break;
	  case 'C': case 'c': ++Ccounts[i+1]; break;
	  case 'G': case 'g': ++Gcounts[i+1]; break;
	  case 'T': case 't': ++Tcounts[i+1]; break;
	  default:
	       printf ("%llu charNum = %d\n", i, consensusSequence[i]); fflush (stdout);
	       ++invalidCounts[i+1]; break; }
     }
     
     fprintf (stderr, "totalSequenceLength = %llu\n", totalSequenceLength);
     fprintf (stderr, "counts at end of consensus sequence: 'A': %llu, 'C': %llu, 'G': %llu, 'T': %llu, invalid: %llu\n", Acounts[totalSequenceLength], Ccounts[totalSequenceLength], Gcounts[totalSequenceLength], Tcounts[totalSequenceLength], invalidCounts[totalSequenceLength]);
     // Now get the read placement info
     infile = Fopen (placementFile, "r");
     while (fgets (line, 100, infile)) {
//	  printf ("%s", (char *) line);
	  getFldsFromLine (line, flds);
//	  printf ("rd = %s\n", flds[0]); fflush (stdout);
	  std::string kUnitig = std::string (flds[1]);
	  int offset = atoi (flds[2]);
	  char ori = flds[3][0];
	  int localLen = defaultReadLength;
	  if (strcmp (readLengthsFile, "") != 0)
	       if (readLengths.find (std::string(flds[0])) != readLengths.end())
		    localLen = readLengths[std::string(flds[0])];
//	  printf ("localLen = %d\n", localLen); fflush (stdout);
	  int beginOffset, endOffset;
	  if (ori == 'F') {
	       beginOffset = offset;
	       endOffset = offset + localLen; }
	  else {
	       endOffset = offset;
	       beginOffset = offset - localLen; }
	  if (beginOffset < 0) // printf ("At 1\n"),
	       beginOffset = 0;
	  if (endOffset > seqLen[kUnitig]) // printf ("At 2, endOffset = %d, seqLen = %d, kUnitig = %s\n", endOffset, seqLen[kUnitig], kUnitig.c_str()),
	       endOffset = seqLen[kUnitig];
	  if (sequenceStart.find(kUnitig) == sequenceStart.end())
	       continue;
	  unsigned long long globalBeginOffset = (sequenceStart[kUnitig]-consensusSequence)  + beginOffset;
	  unsigned long long globalEndOffset = (sequenceStart[kUnitig]-consensusSequence) + endOffset;
//	  printf ("globalEndOffset = %llu\n", globalEndOffset);
	  beginAndEndNetCounts[globalBeginOffset]++;
	  beginAndEndNetCounts[globalEndOffset]--; }
     fclose (infile);
     for (int i=1; i<beginAndEndNetCounts.size(); ++i)
	  coverageCounts[i] = coverageCounts[i-1] + beginAndEndNetCounts[i];
     for (int i=0; i<coverageCounts.size(); ++i)
	  printf ("%d %d T\n", i, (int) coverageCounts[i]);
     
     for (int seqNum=0; seqNum<seqNames.size(); ++seqNum) {
	  std::string seqName = seqNames[seqNum];
	  if (seqLen[seqName] < minSequenceLength)
	       continue;
	  unsigned long long startIndex = sequenceStart[seqName] - consensusSequence;
	  unsigned long long endIndex = startIndex + seqLen[seqName];
	  printf ("startIndex = %llu\n", startIndex);
	  for (unsigned long long i=startIndex; i<=endIndex; ++i) {
	       int count = coverageCounts[i];
	       if ((i>=startIndex+LENGTH_ON_EACH_SIDE_OF_REGION) &&
		   (i+LENGTH_ON_EACH_SIDE_OF_REGION <= endIndex)) {
		    unsigned long long beginOfInterval = i-LENGTH_ON_EACH_SIDE_OF_REGION;
		    unsigned long long endOfInterval = i+LENGTH_ON_EACH_SIDE_OF_REGION;
		    if (invalidCounts[beginOfInterval] != invalidCounts[endOfInterval])
			 continue;
//		    printf ("i = %llu count = %d, lastCount = %d\n", i, count, lastCount);
//	  printf ("i = %d, count = %d, lastCount = %d, beginAndEndNetCounts = %d\n", i, count, lastCount, beginAndEndNetCounts[i]);
		    int GCcount;
		    if (count < 0)
			 printf ("i = %d, count = %d\n", (int) i, (int) count);
		    if (invalidCounts[beginOfInterval] == invalidCounts[endOfInterval]) {
			 int Ccount = Ccounts[endOfInterval] - Ccounts[beginOfInterval];
			 int Gcount = Gcounts[endOfInterval] - Gcounts[beginOfInterval];
			 GCcount = Gcount + Ccount;
			 countsForPctCoveragePairs[intPair (GCcount, count)] += 1;
		    }
	       }
	  }
     }

//     printf ("At 0\n"); fflush (stdout);
     for (std::map<intPair, int>::iterator it=countsForPctCoveragePairs.begin(); it != countsForPctCoveragePairs.end(); ++it)
	  std::cout << (it->first).first << " " << (it->first).second << " " << it->second << '\n';
//     printf ("At 1\n"); fflush (stdout);
     
     return(0);
}

FILE *Fopen (const char *fn, const char *mode)
{
     FILE *result;
     result = fopen (fn, mode);
     if (result == NULL)
     {
          fprintf (stderr, "Couldn't open file '%s' for ", fn);
          switch (mode[0])
          {
          case 'r':
               fprintf (stderr, "reading");
               break;
          case 'w':
               fprintf (stderr, "writing");
               break;
          case 'a':
               fprintf (stderr, "appending");
               break;
          default:
               fprintf (stderr, "unknown operation code '%c'", mode[0]);
               break;
          }
          fprintf (stderr, ". Bye!\n");
          exit (-1);
     }

     return (result);
}

