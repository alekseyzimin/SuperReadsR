// Flags for this program:
// -mean-and-stdev-by-prefix-file fn (the file with mean and stdev for each prefix)
// -unitig-lengths-file fn (the file with the unitig lengths)
// -num-kunitigs-file fn (the file with the largest k-unitig number (+1))
// -overlaps-file fn (the file with the overlaps (of length kmer-1))
// For this exec we are using the unitig numbers starting from 0
// 1 (optional) arg:
// 1) The file of overlaps between k-unitigs and reads
#define NEW_STUFF // Put in to get node-to-node connections
// #define KILLED111115
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "priorityHeap.h"
#include "redBlackTreesInsertOnly.h"

#define NUM_UNITIGS_FILE "numKUnitigs.txt"
#define UNITIG_LENGTHS_FILE "kUnitigLengths.txt"
#define OVERLAPS_FILE "overlap.overlaps"
#define MEAN_AND_STDEV_BY_PREFIX_FILE "meanAndStdevByPrefix.txt"
#define READ_VS_KUNITIG_FILE "testOutput.nucmerLinesOnly.reorderedAndRenamed.txt"

#define DEFAULT_MAX_OFFSET_CONSIDERED_SAME 5
#define MAX_OFFSET_TO_TEST 13000

#define MAX_NODES_ALLOWED 20000

#define NEWLY_ADDED_MAXIMAL_UNITIG 2

#define FRONT_END 1
#define BACK_END 2

struct overlapDataStruct
{
    int unitig1;
    int unitig2;
    int ahg;
    int bhg;
    char ori;
};

struct unitigLocStruct
{
    int unitig2;
    int frontEdgeOffset;
    char ori;
};

// In the following, lastOverlappingOffset gives us the offset for the node
//  where we have the largest overlap into the node
// It is set artificially to the length of the first unitig to start
struct abbrevUnitigLocStruct
{
    int frontEdgeOffset;
    int lastOffsetAtOverlap;
    unsigned short pathNum; 
    char ori;
};

struct unitigPathPrintStruct
{
     int unitig1;
     int frontEdgeOffset;
     int numOverlapsIn;
     int numOverlapsOut;
     char ori;
};

struct kuniToReadMatchStruct
{
     int matchRgnBegin;
     int matchRgnEnd;
     int ahg;
     int bhg;
     int kUnitigMatchBegin;
     int kUnitigMatchEnd;
     int orientedReadMatchBegin;
     int orientedReadMatchEnd;
     int readMatchBegin;
     int readMatchEnd;
     int kUnitigLength;
     int readLength;
     int kUnitigNumber;
     char ori;
} evenReadMatchStructs[10000], oddReadMatchStructs[10000], *kUTRMSptr;
int numEvenReadMatches, numOddReadMatches;
unsigned char matchStructIsUsed[10000];
struct unitigConnectionsForPathStruct
{
     int unitig1;
     int unitig2;
     int frontEdgeOffset1;
     int frontEdgeOffset2;
     char ori1;
     char ori2;
} unitigConnectionsForPathData[1000000];
int numUnitigConnectionsForPathData;

struct augmentedUnitigPathPrintStruct
{
     int unitig1;
     int frontEdgeOffset;
     int numOverlapsIn;
     int numOverlapsOut;
     int beginOffset;
     int endOffset;
     char ori;
} augmentedUnitigPathPrintData[100000];
int numUnitigPathPrintRecsOnPath;


int *startOverlapByUnitig, *unitigLengths;
unsigned char *isUnitigMaximal;
struct overlapDataStruct *overlapData;
struct unitigLocStruct *unitigLocData1, *unitigLocData2;
int maxOffsetToConsiderTheSame;
int *treeReinitList, numTreesUsed;
int *startOverlapIndexByUnitig2, *unitig2OverlapIndex;
int mateUnitig1, mateUnitig2;
unsigned char mateUnitig1ori, mateUnitig2ori;
struct heapStruct priorityQueue;
int curPathNum;
int treeSize;
int minOverlapLength;
int maxDiffInsertSizesForPrinting;
int maxTotAllowableMissingOnEnds;
FILE *outfile;
int reportPaths;
char *flds[1000];
double mean[256][256], stdev[256][256];
char fileName[500];
char rdPrefix[3], rdPrefixHold[3];
long long readNum, readNumHold;
int approxNumPaths;
double insertLengthMeanBetweenKUnisForInsertGlobal, insertLengthStdevGlobal;
// The following keeps track of the distance the 2 read mates are from the
// ends of the k-unitigs at the end
int lengthAdjustment1, lengthAdjustment2;
char superReadName[100000];

int joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev);
FILE *Fopen (const char *fn, const char *mode);
FILE *Popen (const char *fn, const char *mode);
int getFldsFromLine (char *cptrHold);
int getInt (char *fname);
int unitigLocStructCompare (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2);
int unitigLocStructCompareReversed (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2);
// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd);
int findOtherOverlapIndex (int ovlIndex1);
void printIfGood (struct abbrevUnitigLocStruct *ptr);
void completePathPrint (struct abbrevUnitigLocStruct *ptr);
void printPathNode (struct unitigPathPrintStruct *ptr);
int setSuperReadNameFromAugmentedPath (void);
int getSuperReadLength(void);
void funcToGetTreeSize (void *ptr); // Adds 1 to treeSize (a global) each time
void findSingleReadSuperReads(char *readName);
void getSuperReadsForInsert (void);

// RB tree data stuff
struct RBTreeStruct *treeArr, *treeArr2;
struct dataArrayStruct dataArr, dataArr2;
int abbrevLocStructCompForSearch (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2);
int abbrevLocStructCompForSort (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2);
int unitigPathPrintStructComp (struct unitigPathPrintStruct *ptr1,
      struct unitigPathPrintStruct *ptr2);

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

#ifndef updateMatchRecords
#define updateMatchRecords if (readNum % 2 == 0) kUTRMSptr = evenReadMatchStructs+numEvenReadMatches; else kUTRMSptr = oddReadMatchStructs+numOddReadMatches; \
cptr = flds[1]+1; kUTRMSptr->matchRgnBegin = atoi (cptr); \
kUTRMSptr->matchRgnEnd = atoi (flds[2]); \
cptr = flds[3]+1; kUTRMSptr->ahg = atoi (cptr); \
while(*cptr != ',') ++cptr;	++cptr;		       \
kUTRMSptr->bhg = atoi (cptr); \
kUTRMSptr->kUnitigMatchBegin = atoi (flds[4])-1; \
kUTRMSptr->kUnitigMatchEnd = atoi (flds[5]); \
kUTRMSptr->orientedReadMatchBegin = atoi (flds[6]); \
kUTRMSptr->orientedReadMatchEnd = atoi (flds[7]); \
if (kUTRMSptr->orientedReadMatchBegin < kUTRMSptr->orientedReadMatchEnd) {  \
--kUTRMSptr->orientedReadMatchBegin; \
kUTRMSptr->ori = 'F'; \
kUTRMSptr->readMatchBegin = kUTRMSptr->orientedReadMatchBegin; \
kUTRMSptr->readMatchEnd = kUTRMSptr->orientedReadMatchEnd; }    \
else { 	   \
--kUTRMSptr->orientedReadMatchEnd; \
kUTRMSptr->ori = 'R'; \
kUTRMSptr->readMatchBegin = kUTRMSptr->orientedReadMatchEnd; \
kUTRMSptr->readMatchEnd = kUTRMSptr->orientedReadMatchBegin; }    \
kUTRMSptr->kUnitigLength = atoi (flds[8]); \
kUTRMSptr->readLength = atoi (flds[9]); \
kUTRMSptr->kUnitigNumber = atoi (flds[10]); \
if (readNum % 2 == 0) ++numEvenReadMatches; \
else ++numOddReadMatches;
#endif

// #define DEBUG 814
// #define DEBUG 1541
// #define DEBUG 349
// #define DEBUG 727157
// $define DEBUG 41218
// #define DEBUG 25129
// #define DEBUG 5374

int main (int argc, char **argv)
{
    FILE *infile;
    char line[2000], cmd[400], ori, *overlapsFn = OVERLAPS_FILE;
    char *meanAndStdevByPrefixFn = MEAN_AND_STDEV_BY_PREFIX_FILE;
    char *unitigLengthsFn = UNITIG_LENGTHS_FILE;
    char *numKUnitigsFn = NUM_UNITIGS_FILE;
    char *readVsKUnitigFile = READ_VS_KUNITIG_FILE;
    int unitig1, unitig2, overlapCount = 0, itemp, itempHold;
    int unitigNum, numUnitigs, ahg, bhg, firstUnitigNum = 0;
    int i, *iptr;
    float err1, err2;
#if DEBUG
    int unitigForDebugging = DEBUG;
#endif
    int overlapLength;
    char *cptr;
//    int numNewlyAddedMaximalUnitigs;
    int numFilenamesSet = 0;
    int numFlds;

    maxTotAllowableMissingOnEnds = 2;
    minOverlapLength = 40;
    reportPaths = 0;
//    outfile = fopen ("/localraid/tri/out.gaps5", "w");
    outfile = stdout;

    maxDiffInsertSizesForPrinting = 5;
    for (i=1; i<argc; i++) {
	 if (argv[i][0] == '-') {
	      if (strcmp (argv[i], "-max-diff-insert-sizes-for-printing") == 0) {
		   ++i;
		   maxDiffInsertSizesForPrinting = atoi (argv[i]);
	      }
	      else if (strcmp (argv[i], "-report-paths") == 0)
		   reportPaths = 1;
	      else if (strcmp (argv[i], "-min-overlap-length") == 0) {
		   ++i;
		   minOverlapLength = atoi (argv[i]);
	      }
	      else if (strcmp (argv[i], "-mean-and-stdev-by-prefix-file") == 0) {
		   ++i;
		   meanAndStdevByPrefixFn = argv[i];
	      }
	      else if (strcmp (argv[i], "-unitig-lengths-file") == 0) {
		   ++i;
		   unitigLengthsFn = argv[i];
	      }
	      else if (strcmp (argv[i], "-num-kunitigs-file") == 0) {
		   ++i;
		   numKUnitigsFn = argv[i];
	      }
	      else if (strcmp (argv[i], "-overlaps-file") == 0) {
		   ++i;
		   overlapsFn = argv[i];
	      }
	      else { // We need to allow -h later; for now we just exit
		   fprintf (stderr, "Unrecognized flag %s. Bye.\n", argv[i]);
		   return (-1);
	      }
	 }
	 else if (numFilenamesSet == 0) {
	      readVsKUnitigFile = argv[i];
	      ++numFilenamesSet; }
    }
#if DEBUG
    if (argc > 3) unitigForDebugging = atoi (argv[3]); // Must fix
#endif

    rdPrefix[2] = rdPrefixHold[2] = 0;
    strcpy (fileName, meanAndStdevByPrefixFn);
    infile = Fopen (fileName, "r");
    while (fgets (line, 2000, infile)) {
	 getFldsFromLine(line);
	 mean[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[1]);
	 stdev[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[2]);
    }
    fclose (infile);

    mateUnitig1ori = 'F'; mateUnitig2ori = 'R';
// Get the number of unitigs
    strcpy (fileName, numKUnitigsFn);
    numUnitigs = getInt (fileName);
    mallocOrDie (startOverlapByUnitig, numUnitigs + 1 + firstUnitigNum, int);
    mallocOrDie (startOverlapIndexByUnitig2, numUnitigs + 1 + firstUnitigNum, int);

    mallocOrDie (unitigLengths, numUnitigs + firstUnitigNum, int);
    // Here we read in the unitig lengths, allowing for either type of length
    // format
    strcpy (fileName, unitigLengthsFn);
    infile = Fopen (fileName, "r");
    if (! fgets (line, 2000, infile)) {
	 fprintf (stderr, "File %s is of length 0 (or can't be read). Bye!\n", fileName);
	 exit (1); }
    numFlds = getFldsFromLine (line);
    rewind (infile);
    if (numFlds == 1) {
	 int retCode; // For the stupid new compiler
	 for (i = 0, iptr = unitigLengths + firstUnitigNum; i < numUnitigs; i++, iptr++)
	      retCode = fscanf (infile, "%d\n", iptr);
    }
    else {
	 while (fgets (line, 2000, infile)) {
	      getFldsFromLine (line);
	      unitigLengths[atoi(flds[0])] = atoi (flds[1]);
	 }
    }
    fclose (infile);

    mallocOrDie (isUnitigMaximal, numUnitigs + firstUnitigNum, unsigned char);

// Set up space to keep the overlaps file
    sprintf (cmd, "zcat -f %s", overlapsFn);
    infile = Popen (cmd, "r");
    while (fgets (line, 2000, infile))
    {
	sscanf (line, "%d %d ", &unitig1, &unitig2);
	if (unitig1 >= unitig2) 
	     continue;
	overlapCount++;
	startOverlapByUnitig[unitig1]++;
	startOverlapByUnitig[unitig2]++;
    }
    mallocOrDie (overlapData, 2 * overlapCount, struct overlapDataStruct);
    mallocOrDie (unitig2OverlapIndex, 2 * overlapCount, int);
    for (unitigNum = 1; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
    for (unitigNum = 0; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	 startOverlapIndexByUnitig2[unitigNum] = startOverlapByUnitig[unitigNum];

// Set up the priority queue
    initializeEmptyHeap (priorityQueue, struct unitigLocStruct,
	  unitigLocStructCompare);

// Set up the RB trees
    initializeEmptyTrees (treeArr, numUnitigs + 1, dataArr,
	  struct abbrevUnitigLocStruct, abbrevLocStructCompForSort,
	  abbrevLocStructCompForSearch);
    mallocOrDie (treeReinitList, numUnitigs + 1, int);
    
// Set up the RB tree for the final paths
    initializeEmptyTreesWithDataSize (treeArr2, 1, dataArr2,
          struct unitigPathPrintStruct, 40000, unitigPathPrintStructComp,
	  unitigPathPrintStructComp);

// Unitig in the overlaps file
    pclose (infile);
    infile = Popen (cmd, "r");
    while (fgets (line, 2000, infile))
    {
	sscanf (line, "%d %d %c %d %d %f %f ", &unitig1, &unitig2, &ori, &ahg, &bhg,
		&err1, &err2);
	if (unitig1 >= unitig2)
	     continue;
	startOverlapByUnitig[unitig1]--;
	itemp = startOverlapByUnitig[unitig1];
	overlapData[itemp].unitig1 = unitig1;
	overlapData[itemp].unitig2 = unitig2;
	overlapData[itemp].ori = ori;
	overlapData[itemp].ahg = ahg;
	overlapData[itemp].bhg = bhg;
	unitig1 = overlapData[itemp].unitig1;
	if (overlapData[itemp].ahg >= 0)
	     overlapLength = unitigLengths[overlapData[itemp].unitig1] - overlapData[itemp].ahg;
	else
	     overlapLength = unitigLengths[overlapData[itemp].unitig2] + overlapData[itemp].ahg;
	startOverlapIndexByUnitig2[unitig2]--;
	unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig2]] = itemp;
	startOverlapByUnitig[unitig2]--;
	itempHold = itemp;
	itemp = startOverlapByUnitig[unitig2];
	overlapData[itemp].unitig1 = unitig2;
	overlapData[itemp].unitig2 = unitig1;
	overlapData[itemp].ori = ori;
	startOverlapIndexByUnitig2[unitig1]--;
	unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig1]] = itemp;
	if (ori == 'N')
	{
	    overlapData[itemp].ahg = -ahg;
	    overlapData[itemp].bhg = -bhg;
	}
	else
	{
	    overlapData[itemp].ahg = bhg;
	    overlapData[itemp].bhg = ahg;
	}
    }
    pclose (infile);

    sprintf (cmd, "zcat -f %s", readVsKUnitigFile);
    infile = Popen (cmd, "r");
    if (! fgets (line, 2000, infile))
	 return (1); // A critical file doesn't exist
    // Load the appropriate stuff
    numFlds = getFldsFromLine(line);
    cptr = flds[numFlds-1];
    rdPrefixHold[0] = cptr[0];
    rdPrefixHold[1] = cptr[1];
    cptr += 2;
    readNum = readNumHold = atoll (cptr);
    numEvenReadMatches = numOddReadMatches = 0;
    updateMatchRecords;

    while (fgets (line, 2000, infile)) {
	 numFlds = getFldsFromLine(line);
	 cptr = flds[numFlds-1];
	 rdPrefix[0] = cptr[0];
	 rdPrefix[1] = cptr[1];
	 cptr += 2;
	 readNum = atof (cptr);
	 if ((strcmp (rdPrefix, rdPrefixHold) == 0) &&
	     readNum == readNumHold) {
	      // load more data
	      updateMatchRecords;
	      continue;
	 }
	 if ((strcmp (rdPrefix, rdPrefixHold) != 0) ||
	     (readNum != readNumHold+1) ||
	     (readNum % 2 == 0)) {
	      // Get the super-read for the insert we just finished reading
	      getSuperReadsForInsert();
	      // Set up and load the new data
	      numEvenReadMatches = numOddReadMatches = 0;
	      updateMatchRecords;
	      // Update what the old data is
	      strcpy (rdPrefixHold, rdPrefix);
	      readNumHold = readNum;
	      continue;
	 }
	 // If we get here we've gotten to the second read of a mate pair
	 // load the data
	 updateMatchRecords;
	 // hold the updated read info
	 readNumHold = readNum;
    }
    fclose (infile);
    // Output the stuff for the old pair
    getSuperReadsForInsert();

    return (0);
}

int joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev)
{
     
     int lastOffsetToTest = 6000, lastOffsetToTestIfNotMate2, maxOffsetToAllow;
     int j;
     struct unitigLocStruct unitigLocVal;
     struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
     int maxNodes;
     void *vptr;
#if DEBUG
     struct unitigLocStruct *RLPtr;
#endif
     int unitig1, unitig2;
     char ori; 
     int offset;
     int elementIndex;
     int lastPositionOverlapped;
     int overlapLength;
     int ahg, bhg;


     lastOffsetToTest = insertLengthMean+4*insertLengthStdev;
     // The following assumes that all the overlaps are of length
     // minOverlapLength
     lastOffsetToTestIfNotMate2 = lastOffsetToTest - (unitigLengths[mateUnitig2]-minOverlapLength);
     // Adjust overlaps for mateUnitig1 if ori not 'F' and for
     //    mateUnitig2 if ori not 'R' to what they would be if they had
     //    the desired orientation. Note that mateUnitig1 is only a
     //    source and mateUnitig2 is only a sink
     unitigLocVal.unitig2 = mateUnitig1;
     unitigLocVal.frontEdgeOffset = unitigLengths[mateUnitig1];
     unitigLocVal.ori = mateUnitig1ori;
     numTreesUsed = 0;
     abbrevUnitigLocVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
     abbrevUnitigLocVal.ori = unitigLocVal.ori;
     // We may want to change the following
     abbrevUnitigLocVal.lastOffsetAtOverlap =
	  abbrevUnitigLocVal.frontEdgeOffset;
     if (treeArr[mateUnitig1].root == TREE_NIL)
     {
	  treeReinitList[numTreesUsed] = mateUnitig1;
	  numTreesUsed++;
     }
     abbrevUnitigLocVal.pathNum = 0;
     RBTreeInsertElement (treeArr + mateUnitig1, &abbrevUnitigLocVal);
#if 0
     printf ("Inserting in the RB tree at %d: fEO = %d, lOAO = %d, pN = %u ori = %c\n", mateUnitig1, abbrevUnitigLocVal.frontEdgeOffset, abbrevUnitigLocVal.lastOffsetAtOverlap, abbrevUnitigLocVal.pathNum, abbrevUnitigLocVal.ori);
#endif     
     heapInsert (&priorityQueue, &unitigLocVal);
     maxNodes = 1;
//     printf ("Got to 30\n");
     while (!heapIsEmpty (&priorityQueue))
     {
	  if (priorityQueue.heapSize > maxNodes)
	       maxNodes = priorityQueue.heapSize;
#if DEBUG
	  printf ("Starting new offset: "); fflush (stdout);
	  for (j = 0; j < priorityQueue.heapSize; j++)
	  {
	       setHeapValPtr (vptr, &priorityQueue, j);
	       RLPtr = vptr;
	       printf ("%d ", RLPtr->frontEdgeOffset); fflush (stdout);
	  }
	  printf ("\n"); fflush (stdout);
#endif
	  heapExtractRoot (&priorityQueue, &unitigLocVal);
	  unitig1 = unitigLocVal.unitig2;
#if DEBUG
	  printf ("unitig1 = %d, ", unitig1); fflush (stdout);
#endif
	  ori = unitigLocVal.ori;
	  offset = unitigLocVal.frontEdgeOffset;
	  abbrevUnitigLocVal.frontEdgeOffset = offset;
#if DEBUG
	  printf ("offset = %d, ", offset); fflush (stdout);
#endif
	  abbrevUnitigLocVal.ori = ori;
	  elementIndex =
	       treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
	  assert (elementIndex != TREE_NIL);
	  setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	  abbRLPtr = vptr;
	  lastPositionOverlapped = abbRLPtr->lastOffsetAtOverlap;
#if DEBUG
	  printf ("offset = %d, lastPosOvlppd = %d, ori = %c, ", offset,
		  lastPositionOverlapped, ori);  fflush (stdout);
#endif
//	  printf ("Got to 40\n");
	  for (j = startOverlapByUnitig[unitig1];
	       j < startOverlapByUnitig[unitig1 + 1]; j++)
	  {
#if DEBUG
	       printf ("Starting an overlap:");
#endif
	       unitig2 = overlapData[j].unitig2;
	       if (overlapData[j].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[j].ahg;
	       else
		    overlapLength = unitigLengths[unitig2] + overlapData[j].ahg;
	       if (overlapLength > unitigLengths[unitig1])
		    overlapLength = unitigLengths[unitig1];
	       if (overlapLength > unitigLengths[unitig2])
		    overlapLength = unitigLengths[unitig2];
#if DEBUG
	       printf ("; ovl len = %d", overlapLength);
#endif
#if 0		
	       if (overlapLength < minOverlapLength) {
#if DEBUG
		    printf (" Ovl too short\n");
#endif
		    continue;
	       }
#endif
#if DEBUG
	       printf ("; unitig2 = %d, ori = %c, ahg = %d, bhg = %d\n", unitig2,
		       overlapData[j].ori, overlapData[j].ahg,
		       overlapData[j].bhg);
#endif
//	       printf ("Got to 50, unitig1 = %d, unitig2 = %d\n", unitig1, unitig2);
	       if (ori == 'F')
	       {
		    if (overlapData[j].bhg <= 0) 
			 continue;
		    bhg = overlapData[j].bhg;
		    if (overlapData[j].ori == 'N')
			 abbrevUnitigLocVal.ori = 'F';
		    else
			 abbrevUnitigLocVal.ori = 'R';
		    abbrevUnitigLocVal.frontEdgeOffset = offset + bhg;
	       }
	       else
	       {
		    if (overlapData[j].ahg >= 0)
			 continue;
		    ahg = overlapData[j].ahg;
		    if (overlapData[j].ori == 'N')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
		    abbrevUnitigLocVal.frontEdgeOffset = offset - ahg;
	       }
	       abbrevUnitigLocVal.lastOffsetAtOverlap = offset;
	       // Skip if the offset is too large
//	       printf ("frontEdgeOffset = %d, lastOffsetToTest = %d\n", abbrevUnitigLocVal.frontEdgeOffset, lastOffsetToTest);
	       if (unitig2 == mateUnitig2)
		    maxOffsetToAllow = lastOffsetToTest;
	       else
		    maxOffsetToAllow = lastOffsetToTestIfNotMate2;
	       if (abbrevUnitigLocVal.frontEdgeOffset > maxOffsetToAllow)
		    continue;
#if DEBUG
	       printf ("cur front = %d\n", abbrevUnitigLocVal.frontEdgeOffset);
#endif
	       // Skip if abbrevUnitigLocVal alunitigy seen for unitig2
	       elementIndex =
		    treeFindElement (treeArr + unitig2, &abbrevUnitigLocVal);
	       setTreeValPtr (vptr, treeArr + unitig2, elementIndex);
	       abbRLPtr = vptr;
	       abbRLPtr->lastOffsetAtOverlap = offset;
	       
//	       printf ("Got to 60\n");
	       if (elementIndex != TREE_NIL)
		    continue;
#if DEBUG
	       printf ("Adding to the tree\n");
#endif
	       // Insert this value in the priority queue
	       unitigLocVal.unitig2 = unitig2;
	       unitigLocVal.frontEdgeOffset = abbrevUnitigLocVal.frontEdgeOffset;
	       unitigLocVal.ori = abbrevUnitigLocVal.ori;
//	       printf ("Got to 70, unitig = %d\n", unitigLocVal.unitig2);
	       {
		    heapInsert (&priorityQueue, &unitigLocVal); }
	       if (treeArr[unitig2].root == TREE_NIL)
	       {
		    // If unitig2's tree never seen before
		    //    Add unitig2 to list of trees to reinit
		    treeReinitList[numTreesUsed] = unitig2;
		    numTreesUsed++;
	       }
	       // Add offset to list for tree
	       abbrevUnitigLocVal.pathNum = 0;
	       RBTreeInsertElement (treeArr + unitig2, &abbrevUnitigLocVal);
#if DEBUG
//		if ((unitig2 == mateUnitig2) && (abbrevUnitigLocVal.ori == 'R'))
	       if (unitig2 == mateUnitig2)
		    printf ("Adding distance %d for the (rev oriented) mate, unitig %d\n",
			    abbrevUnitigLocVal.frontEdgeOffset, mateUnitig2);
#endif		
	       //   Make sure the root of the tree is updated (if needed)
	  }			// End of going through overlaps for unitig
	  if (maxNodes > MAX_NODES_ALLOWED)
	       break;
     }			// Ends heapIsEmpty line
     // Do output for the mate pair
     //      for (j=0; j<=numUnitigs; j++)
     curPathNum = 0;
     approxNumPaths = 0;
//	       printf ("Got to 80\n");
     if (treeArr[mateUnitig2].root != TREE_NIL)
     {
	  treeSize = 0;
	  inOrderTreeWalk (treeArr + mateUnitig2, treeArr[mateUnitig2].root,
			   (void (*)()) funcToGetTreeSize);
#ifdef KILLED111115
	  printf ("treeSize = %d\n", treeSize);
#ifndef NO_OUTPUT
	  fprintf (outfile, "mateUnitig2 = %d\n", mateUnitig2); // This prints
#endif
	  fprintf (outfile, "maxNodes = %d\n", maxNodes); // This prints
#endif
	  priorityQueue.compare = (void *) unitigLocStructCompareReversed;
	  // This is where the main print statement is
	  if (maxNodes <= MAX_NODES_ALLOWED)
	       inOrderTreeWalk (treeArr + mateUnitig2, treeArr[mateUnitig2].root,
//				  (void (*)(void *)) printIfGood);
				(void (*)(void *)) completePathPrint);
	  priorityQueue.compare = (void *) unitigLocStructCompare;
     }
     priorityQueue.heapSize = 0;
     for (j = 0; j < numTreesUsed; j++)
	  treeArr[treeReinitList[j]].root = TREE_NIL;
     numTreesUsed = 0;
     dataArr.arraySize = 0;

#ifdef KILLED111115
     printf ("Approx num paths returned = %d\n", approxNumPaths);
#endif
     return (approxNumPaths);
}

// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd)
{
     int thisUnitig, otherUnitig, overlapLength;
     otherUnitig = overlapData[maxOvlIndex].unitig2;
     thisUnitig = overlapData[j].unitig2;
     if (whichEnd == FRONT_END)
	  overlapLength = unitigLengths[thisUnitig] - overlapData[j].ahg;
     else
	  overlapLength = unitigLengths[thisUnitig] + overlapData[j].bhg;
     if (overlapLength < maxOvlLen) return (-1);
     if (overlapLength > maxOvlLen) return (overlapLength);
     if (unitigLengths[otherUnitig] > unitigLengths[thisUnitig]) return (-1);
     if (unitigLengths[otherUnitig] < unitigLengths[thisUnitig]) return (overlapLength);
     if (otherUnitig < thisUnitig) return (-1);
     if (otherUnitig > thisUnitig) return (overlapLength);
     if (overlapData[maxOvlIndex].ori == 'N') return (-1);
     return (overlapLength);
}

int findOtherOverlapIndex (int ovlIndex1)
{
     int unitig1, unitig2, itemp;
     unitig1 = overlapData[ovlIndex1].unitig1;
     unitig2 = overlapData[ovlIndex1].unitig2;
     for (itemp=startOverlapByUnitig[unitig2]; itemp<startOverlapByUnitig[unitig2+1]; itemp++) {
	  if (overlapData[itemp].unitig2 != unitig1) continue;
	  if (overlapData[itemp].ori != overlapData[ovlIndex1].ori) continue;
	  if (overlapData[itemp].ori == 'N') {
	       if (overlapData[itemp].ahg + overlapData[ovlIndex1].ahg != 0)
		    continue;
	  }
	  else {
	       if (overlapData[itemp].ahg != overlapData[ovlIndex1].bhg)
		    continue;
	  }
	  // If we get here they agree
	  break;
     }
     return (itemp);
}

void printIfGood (struct abbrevUnitigLocStruct *ptr)
{
     int val;
     if (ptr->ori == 'R') {
	  val = ptr->frontEdgeOffset;
	  if (mateUnitig1ori != 'F') val -= unitigLengths[mateUnitig1];
	  if (mateUnitig2ori != 'R') val -= unitigLengths[mateUnitig2];
#ifndef NO_OUTPUT
	  fprintf (outfile, "%d\n", val);
#endif
     }
}

void printPathNode (struct unitigPathPrintStruct *ptr)
{
     int beginOffset, endOffset;
     if (ptr->ori == 'F') {
	  endOffset = ptr->frontEdgeOffset;
	  beginOffset = endOffset - unitigLengths[ptr->unitig1];
     }
     else {
	  beginOffset = ptr->frontEdgeOffset;
	  endOffset = beginOffset - unitigLengths[ptr->unitig1];
     }
//     fprintf (outfile, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", ptr->unitig1, ptr->frontEdgeOffset, ptr->ori, beginOffset, endOffset, ptr->numOverlapsIn, ptr->numOverlapsOut);
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = ptr->unitig1;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = ptr->frontEdgeOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = ptr->ori;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].beginOffset = beginOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].endOffset = endOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = ptr->numOverlapsIn;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = ptr->numOverlapsOut;
     ++numUnitigPathPrintRecsOnPath;
     if (ptr->numOverlapsIn == 0)
	  approxNumPaths += ptr->numOverlapsOut;
     else if (ptr->numOverlapsOut > 0)
	  approxNumPaths += (ptr->numOverlapsOut - 1);
}

void completePathPrint (struct abbrevUnitigLocStruct *ptr)
{
     struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
     struct unitigLocStruct unitigLocVal;
     struct unitigPathPrintStruct unitigPathPrintVal, *rppvPtr1, *rppvPtr2;
     char ori;
#ifdef NEW_STUFF
     char tempOri1, tempOri2;
#endif
     int isSpecialCase, finalOffset, minConnectingOffset, i, index;
     int unitig1, unitig2, offset, elementIndex1, elementIndex, overlapLength;
     int minConnectingUnitig=0, minConnectingOverlapIndex;
     char minConnectingOri=' ';
     void *vptr;
     double numStdevsFromMean;
     // In the following we assume we move from left to right when moving from
     // mateUnitig1 to mateUnitig2.
     // First we want to find out if mateUnitig2 has an overlap which ends to
     // the left of the right-most base of mateUnitig2. If not, we put the
     // unitig which overlaps mateUnitig2 and extends the least to the left as
     // the last node of the priority queue (but we have to adjust for the 
     // missing connection between this unitig and mateUnitig2 later.)
     // Otherwise, we put mateUnitig2 on the priority queue and start.
     // Note that at this point mateUnitig2 have reverse orientation; even if
     // it was specified as forward orientation in the pair, it has been changed
     // to reverse orientation. Similarly, mateUnitig1 has forward orientation.
     ++curPathNum;
#if 0
     printf ("curPathNum = %d, nodePathNum = %d; ", curPathNum, ptr->pathNum);
     printf ("frontEdgeOffset = %d, lastOffsetAtOverlap = %d, ori = %c\n", ptr->frontEdgeOffset, ptr->lastOffsetAtOverlap, ptr->ori);
#endif
     ori = ptr->ori;
     if (ori != mateUnitig2ori) return;
//     printf ("In rtn completePathPrint\n");
     isSpecialCase = 1;
     finalOffset = ptr->frontEdgeOffset;
     minConnectingOffset = finalOffset + 1000000;
     numStdevsFromMean = (finalOffset - insertLengthMeanBetweenKUnisForInsertGlobal)/insertLengthStdevGlobal;
#ifdef KILLED111115
     fprintf (outfile, "%d %f\n", finalOffset, numStdevsFromMean);
#endif
     for (i=startOverlapIndexByUnitig2[mateUnitig2]; i<startOverlapIndexByUnitig2[mateUnitig2+1]; i++) {
	  index = unitig2OverlapIndex[i];
	  unitig1 = overlapData[index].unitig1;
	  if (overlapData[index].ahg >= 0)
	       overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	  else
	       overlapLength = unitigLengths[mateUnitig2] + overlapData[index].ahg;
	  if (overlapLength > unitigLengths[unitig1])
	       overlapLength = unitigLengths[unitig1];
	  if (overlapLength > unitigLengths[mateUnitig2])
	       overlapLength = unitigLengths[mateUnitig2];
#if 0
	  if (overlapLength < minOverlapLength)
	       continue;
#endif
	  if (overlapData[index].ori == 'N')
	       abbrevUnitigLocVal.ori = ori;
	  else {
	       if (ori == 'F') abbrevUnitigLocVal.ori = 'R';
	       else abbrevUnitigLocVal.ori = 'F';
	  }
	  if (abbrevUnitigLocVal.ori == 'F')
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset - overlapData[index].bhg;
	  else
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset + overlapData[index].ahg;
//	  printf ("We are at 1\n");
	  elementIndex = treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
	  if (elementIndex == TREE_NIL) continue;
//	  printf ("We are at 2\n");
	  setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	  abbRLPtr = vptr;
#if 0
	  printf ("frontEdgeOffset = %d\n", abbRLPtr->frontEdgeOffset);
#endif
	  if (abbRLPtr->frontEdgeOffset < finalOffset) {
	       isSpecialCase = 0;
	       break;
	  }
	  else if (abbRLPtr->frontEdgeOffset < minConnectingOffset) {
	       minConnectingOffset = abbRLPtr->frontEdgeOffset;
	       minConnectingUnitig = unitig1;
	       minConnectingOri = abbrevUnitigLocVal.ori;
	       minConnectingOverlapIndex = index;
	  }
     }
     if (isSpecialCase) {
//	  fprintf (stderr, "We shouldn't get here\n");
	  fprintf (stdout, "We shouldn't get here\n");
	  unitigLocVal.unitig2 = minConnectingUnitig;
	  unitigLocVal.frontEdgeOffset = minConnectingOffset;
	  unitigLocVal.ori = minConnectingOri;
	  // ..and for the path
	  unitigPathPrintVal.unitig1 = mateUnitig2;
	  unitigPathPrintVal.numOverlapsIn = 1;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = 'R'; // Forced; may be adjusted later
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
	  unitigPathPrintVal.unitig1 = minConnectingUnitig;
	  unitigPathPrintVal.frontEdgeOffset = minConnectingOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 1;
	  unitigPathPrintVal.ori = minConnectingOri;
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
     }
     else {
	  unitigLocVal.unitig2 = mateUnitig2;
	  unitigLocVal.frontEdgeOffset = finalOffset;
	  unitigLocVal.ori = mateUnitig2ori;
	  unitigPathPrintVal.unitig1 = mateUnitig2;
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = mateUnitig2ori;
#if 0
	  printf ("Inserting unitig1 = %d, offset = %d, ori = %c\n", unitigPathPrintVal.unitig1, unitigPathPrintVal.frontEdgeOffset, unitigPathPrintVal.ori);
#endif
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
     }
#if 0
     printf ("isSpecialCase = %d, unitigLocVal = %d, %d, %c\n", isSpecialCase, mateUnitig2, finalOffset, unitigLocVal.ori);
#endif
     numUnitigConnectionsForPathData = 0;
     heapInsert (&priorityQueue, &unitigLocVal);
     while (!heapIsEmpty (&priorityQueue)) {
	  heapExtractRoot (&priorityQueue, &unitigLocVal);
	  unitig2 = unitigLocVal.unitig2;
	  offset = unitigLocVal.frontEdgeOffset;
	  ori = unitigLocVal.ori;
	  unitigPathPrintVal.unitig1 = unitig2;
	  unitigPathPrintVal.frontEdgeOffset = offset;
	  unitigPathPrintVal.ori = ori;
	  elementIndex1 = treeFindElement (treeArr2, &unitigPathPrintVal);
	  setTreeValPtr (vptr, treeArr2, elementIndex1);
	  rppvPtr1 = vptr;
#if 0
	  printf ("unitig2 = %d, offset = %d, ori = %c; elementIndex1 = %d\n", unitig2, offset, ori, elementIndex1);
#endif
	  for (i=startOverlapIndexByUnitig2[unitig2]; i<startOverlapIndexByUnitig2[unitig2+1]; i++) {
	       index = unitig2OverlapIndex[i];
	       unitig1 = overlapData[index].unitig1;
#if 0
	       printf ("unitig1 = %d\n", unitig1);
#endif
	       if (overlapData[index].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	       else 
		    overlapLength = unitigLengths[unitig2] + overlapData[index].ahg;
	       if (overlapLength > unitigLengths[unitig1])
		    overlapLength = unitigLengths[unitig1];
	       if (overlapLength > unitigLengths[unitig2])
		    overlapLength = unitigLengths[unitig2];
#if 0
	       if (overlapLength < minOverlapLength)
		    continue;
#endif
	       if (overlapData[index].ori == 'N')
		    abbrevUnitigLocVal.ori = ori;
	       else {
		    if (ori == 'F')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
	       }
	       if (abbrevUnitigLocVal.ori == 'F')
		    abbrevUnitigLocVal.frontEdgeOffset = offset - overlapData[index].bhg;
	       else
		    abbrevUnitigLocVal.frontEdgeOffset = offset + overlapData[index].ahg;
	       elementIndex = treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
//	       printf ("Got to 21, unitig1 = %d\n", unitig1);
//	       printf ("fEO = %d, lOAO = %d, pN = %u ori = %c\n", abbrevUnitigLocVal.frontEdgeOffset, abbrevUnitigLocVal.lastOffsetAtOverlap, abbrevUnitigLocVal.pathNum, abbrevUnitigLocVal.ori);
	       if (elementIndex == TREE_NIL) continue;
//	       printf ("Got to 215\n");
	       setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	       abbRLPtr = vptr;
//	       printf ("Got to 22\n");
	       if (abbRLPtr->frontEdgeOffset >= offset) continue;
//	       printf ("Got to 23, abbRLVpathNum = %d, curPathNum = %d\n", abbRLPtr->pathNum, curPathNum);
	       // It hasn't been seen in the retrace, so put on the queue
	       if (abbRLPtr->pathNum < curPathNum) {
#if 0
		    printf ("Adding node: unitig2 = %d, offset = %d, ori = %c\n", unitigLocVal.unitig2, unitigLocVal.frontEdgeOffset, unitigLocVal.ori);
#endif
		    abbRLPtr->pathNum = curPathNum;
		    unitigLocVal.unitig2 = unitig1;
		    unitigLocVal.frontEdgeOffset = abbRLPtr->frontEdgeOffset;
		    unitigLocVal.ori = abbRLPtr->ori;
		    heapInsert (&priorityQueue, &unitigLocVal);
		    unitigPathPrintVal.unitig1 = unitig1;
		    unitigPathPrintVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
		    unitigPathPrintVal.ori = unitigLocVal.ori;
		    unitigPathPrintVal.numOverlapsIn = 0;
		    unitigPathPrintVal.numOverlapsOut = 0;
		    RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
	       }
		    
	       unitigPathPrintVal.unitig1 = unitig1;
	       unitigPathPrintVal.frontEdgeOffset = abbRLPtr->frontEdgeOffset;
	       unitigPathPrintVal.ori = abbRLPtr->ori;
	       elementIndex = treeFindElement (treeArr2, &unitigPathPrintVal);
	       setTreeValPtr (vptr, treeArr2, elementIndex);
	       rppvPtr2 = vptr;
	       ++(rppvPtr2->numOverlapsOut);
	       // The following must be recalced in case the array had to be
	       // moved due to needing more space
	       setTreeValPtr (vptr, treeArr2, elementIndex1);
	       rppvPtr1 = vptr;
	       ++(rppvPtr1->numOverlapsIn);
#ifdef NEW_STUFF
	       tempOri1 = rppvPtr1->ori;
	       if (rppvPtr1->unitig1 == mateUnitig1) tempOri1 = mateUnitig1ori;
	       if (rppvPtr1->unitig1 == mateUnitig2) tempOri1 = mateUnitig2ori;
	       tempOri2 = rppvPtr2->ori;
	       if (rppvPtr2->unitig1 == mateUnitig1) tempOri2 = mateUnitig1ori;
	       if (rppvPtr2->unitig1 == mateUnitig2) tempOri2 = mateUnitig2ori;
//	       printf ("Node (%d, %d, %c) -> (%d, %d, %c)\n", rppvPtr2->unitig1, rppvPtr2->frontEdgeOffset, tempOri2, rppvPtr1->unitig1, rppvPtr1->frontEdgeOffset, tempOri1);
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig1 = rppvPtr2->unitig1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig2 = rppvPtr1->unitig1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset1 = rppvPtr2->frontEdgeOffset;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset2 = rppvPtr1->frontEdgeOffset;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori1 = tempOri2;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori2 = tempOri1;
	       ++numUnitigConnectionsForPathData;
#endif
//	       printf ("Got to 24\n");
	       
	  }
     }
#ifdef KILLED111115
     for (i=0; i<numUnitigConnectionsForPathData; i++)
	  printf ("Node (%d, %d, %c) -> (%d, %d, %c)\n", unitigConnectionsForPathData[i].unitig1, unitigConnectionsForPathData[i].frontEdgeOffset1, unitigConnectionsForPathData[i].ori1, unitigConnectionsForPathData[i].unitig2, unitigConnectionsForPathData[i].frontEdgeOffset2, unitigConnectionsForPathData[i].ori2);
#endif
#if 0
     printf ("tree root = %d\n", treeArr2[0].root);
#endif
     numUnitigPathPrintRecsOnPath = 0;
     if ((treeSize <= maxDiffInsertSizesForPrinting) && (reportPaths))
	  inOrderTreeWalk (treeArr2, treeArr2[0].root,
			   (void (*)(void *)) printPathNode);
#ifdef KILLED111115
     for (i=0; i<numUnitigPathPrintRecsOnPath; i++)
	  fprintf (outfile, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].frontEdgeOffset, augmentedUnitigPathPrintData[i].ori, augmentedUnitigPathPrintData[i].beginOffset, augmentedUnitigPathPrintData[i].endOffset, augmentedUnitigPathPrintData[i].numOverlapsIn, augmentedUnitigPathPrintData[i].numOverlapsOut);
#endif
     if (approxNumPaths == 1) {
	  int isReversed, superReadLength;
	  // the following uses augmentedUnitigPathPrintData
	  superReadLength = getSuperReadLength ();
	  isReversed = setSuperReadNameFromAugmentedPath ();
	  printf ("%s%lld %s ", rdPrefixHold, readNumHold-1, superReadName);
	  if (! isReversed)
	       printf ("%d F\n", lengthAdjustment1);
	  else
	       printf ("%d R\n", superReadLength - lengthAdjustment1);
	  printf ("%s%lld %s ", rdPrefixHold, readNumHold, superReadName);
	  if (! isReversed)
	       printf ("%d R\n", superReadLength - lengthAdjustment2);
	  else
	       printf ("%d F\n", lengthAdjustment2);
     }
#if 0
     printf ("final offset = %d, arraySize = %d\n", finalOffset, dataArr2.arraySize);
#endif
     treeArr2[0].root = TREE_NIL;
     dataArr2.arraySize = 0;
}

int setSuperReadNameFromAugmentedPath (void)
{
     int isReversed=0, i;
     char *cptr;
     for (i=0; i<numUnitigPathPrintRecsOnPath/2; i++) {
	  if (augmentedUnitigPathPrintData[i].unitig1 != augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1) {
	       if (augmentedUnitigPathPrintData[i].unitig1 < augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (augmentedUnitigPathPrintData[i].ori == augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].ori) {
	       if (augmentedUnitigPathPrintData[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
     }
     cptr = superReadName;
     if (isReversed == 0) {
	  sprintf (cptr, "%d%c", augmentedUnitigPathPrintData[0].unitig1, augmentedUnitigPathPrintData[0].ori);
	  cptr += strlen (cptr);
	  for (i=1; i<numUnitigPathPrintRecsOnPath; i++) {
	       sprintf (cptr, "_%d_%d%c", minOverlapLength, augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].ori);
	       cptr += strlen (cptr);
	  }
     }
     else {
	  sprintf (cptr, "%d%c", augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].unitig1, (augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].ori == 'F') ? 'R' : 'F');
	  cptr += strlen (cptr);
	  for (i=numUnitigPathPrintRecsOnPath-2; i>=0; i--) {
	       sprintf (cptr, "_%d_%d%c", minOverlapLength, augmentedUnitigPathPrintData[i].unitig1, (augmentedUnitigPathPrintData[i].ori == 'F') ? 'R' : 'F');
	       cptr += strlen (cptr);
	  }
     }
     return (isReversed);
}

int getSuperReadLength(void)
{
     int totLen, i;
     totLen = unitigLengths[augmentedUnitigPathPrintData[0].unitig1];
     for (i=1; i<numUnitigPathPrintRecsOnPath; i++)
	  totLen += (unitigLengths[augmentedUnitigPathPrintData[i].unitig1] - minOverlapLength);

     return (totLen);
}

void funcToGetTreeSize (void *ptr)
{
     ++treeSize;
}

int unitigLocStructCompare (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2)
{
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (1);
    if (ptr1->unitig2 > ptr2->unitig2)
	return (-1);
    if (ptr1->unitig2 < ptr2->unitig2)
	return (1);
    if (ptr1->ori > ptr2->ori)
	return (-1);
    if (ptr1->ori < ptr2->ori)
	return (1);
    return (0);
}

int unitigLocStructCompareReversed (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2)
{
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (1);
    if (ptr1->unitig2 < ptr2->unitig2)
	return (-1);
    if (ptr1->unitig2 > ptr2->unitig2)
	return (1);
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    return (0);
}

int abbrevLocStructCompForSort (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2)
{
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (1);
    return (0);
}

int abbrevLocStructCompForSearch (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2)
{
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    if (ptr1->frontEdgeOffset <
	  ptr2->frontEdgeOffset - DEFAULT_MAX_OFFSET_CONSIDERED_SAME)
	return (-1);
    if (ptr1->frontEdgeOffset >
	  ptr2->frontEdgeOffset + DEFAULT_MAX_OFFSET_CONSIDERED_SAME)
	return (1);
    return (0);
}

int unitigPathPrintStructComp (struct unitigPathPrintStruct *ptr1,
      struct unitigPathPrintStruct *ptr2)
{
     if (ptr1->unitig1 == mateUnitig2) {
	  if (ptr1->unitig1 == ptr2->unitig1)
	       return (0);
	  else
	       return(1);
     }
     if (ptr2->unitig1 == mateUnitig2) return(-1);
     if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset) return (-1);
     if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset) return (1);
     if (ptr1->unitig1 < ptr2->unitig1) return (-1);
     if (ptr1->unitig1 > ptr2->unitig1) return (1);
     if (ptr1->ori < ptr2->ori) return (-1);
     if (ptr1->ori > ptr2->ori) return (1);
     return (0);
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

FILE *Popen (const char *fn, const char *mode)
{
    FILE *result;
    result = popen (fn, mode);
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

int getInt (char *fname)
{
    FILE *infile;
    int tval;

    infile = Fopen (fname, "r");
    if (! fscanf (infile, "%d\n", &tval)) {
	 fprintf (stderr, "Couldn't read file %s. Bye!\n", fname);
	 exit (1);
    }
    fclose (infile);
    return (tval);
}

int getFldsFromLine (char *cptrHold)
{
     int numFlds=0, state = 0;
     char *cptr;

     for (cptr=cptrHold; *cptr; cptr++) {
          if (isspace (*cptr)) { state = 0;*cptr = 0; }
          else {
               if (state == 1) continue;
               flds[numFlds] = cptr;
               ++numFlds;
               state = 1;
          }
     }
     return (numFlds);
}

void findSingleReadSuperReads(char *readName)
{
     long long tempInt;
     char *cptr=readName+2;
     int countOfMatchingKUnitigs, offsetOfReadInSuperRead;
     int minReadOffset, maxReadOffset, minReadOffsetSeen, maxReadOffsetSeen;
     int i, j, recNumToUse=0;
     int isReversed=0;

     tempInt = atoll(cptr);
#ifdef KILLED111115
     printf ("findSingleReadSuperReads\n");
#endif
     if (tempInt % 2 == 0) {
	  countOfMatchingKUnitigs = numEvenReadMatches;
	  kUTRMSptr = &(evenReadMatchStructs[0]);
     }
     else {
	  countOfMatchingKUnitigs = numOddReadMatches;
	  kUTRMSptr = &(oddReadMatchStructs[0]);
     }
     
//     printf ("countOfMatchingKUnitigs = %d\n", countOfMatchingKUnitigs);
     i = 0;
     minReadOffsetSeen = kUTRMSptr[i].readMatchBegin;
     maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
     matchStructIsUsed[i] = 1;
     for (i=1; i<countOfMatchingKUnitigs; i++) {
	  matchStructIsUsed[i] = 0;
	  if (kUTRMSptr[i].readMatchEnd <= maxReadOffsetSeen)
	       continue;
	  if (kUTRMSptr[i].readMatchBegin < maxReadOffsetSeen-minOverlapLength)
	       continue; // Otherwise the k-unitigs overlap too much
	  if (kUTRMSptr[i].readMatchBegin > maxReadOffsetSeen)
	       return; // Part of the middle of the read is uncovered by k-unis
	  maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
	  matchStructIsUsed[i] = 1;
//	  printf ("Struct number %d is used, kUnitig = %d\n", i, kUTRMSptr[i].kUnitigNumber);
     }
     if (minReadOffsetSeen + (kUTRMSptr[0].readLength - maxReadOffsetSeen) > maxTotAllowableMissingOnEnds)
	  return;
     
     i=-1; j=countOfMatchingKUnitigs;
     isReversed = 0;
     while (1) {
	  ++i; --j;
	  while (!matchStructIsUsed[i])
	       ++i;
	  while (!matchStructIsUsed[j])
	       --j;
	  if (kUTRMSptr[i].kUnitigNumber != kUTRMSptr[j].kUnitigNumber) {
	       if (kUTRMSptr[i].kUnitigNumber < kUTRMSptr[j].kUnitigNumber)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (kUTRMSptr[i].ori == kUTRMSptr[j].ori) {
	       if (kUTRMSptr[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (j<=i)
	       break;
     }
     cptr = superReadName;
     if (isReversed == 0) {
	  for (i=0; 1; i++)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf (cptr, "%d%c", kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	  cptr += strlen(cptr);
	  maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  for (++i; i<countOfMatchingKUnitigs; i++) {
	       if (! matchStructIsUsed[i])
		    continue;
	       sprintf (cptr, "_%d_%d%c", maxReadOffset-kUTRMSptr[i].readMatchBegin, kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	       cptr += strlen(cptr);
	       maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchBegin - kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber] - kUTRMSptr[recNumToUse].kUnitigMatchEnd - kUTRMSptr[i].readMatchBegin;
//	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, kUTRMSptr[recNumToUse].ori);
	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, 'F');	  
     }
     else { // The k-unitigs are reversed from those reported
	  for (i=countOfMatchingKUnitigs-1; 1; i--)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf (cptr, "%d%c", kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	  cptr += strlen(cptr);
	  minReadOffset = kUTRMSptr[i].readMatchBegin;
	  for (--i; i>=0; i--) {
	       if (! matchStructIsUsed[i])
		    continue;
	       sprintf (cptr, "_%d_%d%c", kUTRMSptr[i].readMatchEnd-minReadOffset, kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	       cptr += strlen(cptr);
	       minReadOffset = kUTRMSptr[i].readMatchBegin;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = (unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber]-kUTRMSptr[recNumToUse].kUnitigMatchBegin) + kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchEnd + kUTRMSptr[recNumToUse].readMatchBegin;
	  // The k-unitigs are reversed from those reported
//	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, (kUTRMSptr[recNumToUse].ori == 'F') ? 'R' : 'F');
	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, 'R');
     }
//     printf ("At 50\n");
}

void getSuperReadsForInsert (void)
{
     char readNameSpace[200];
     int insertLengthMean;

     // Output the stuff for the old pair
#ifdef KILLED111115
     printf ("%s%lld %d %d\n", rdPrefixHold, readNumHold, numEvenReadMatches, numOddReadMatches);
#endif
     sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
     if ((numEvenReadMatches == 0) || (numOddReadMatches == 0)) {
	  findSingleReadSuperReads(readNameSpace);
	  return; }
     // If we get here both the even read and the odd read have
     // matches to k-unitigs
     // The next takes care of the case where both the source and
     // destination k-unitig are the same (We don't join in this case)
     mateUnitig1 = evenReadMatchStructs[0].kUnitigNumber;
     mateUnitig2 = oddReadMatchStructs[0].kUnitigNumber;
     if (mateUnitig1 == mateUnitig2) {
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads (readNameSpace);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads (readNameSpace);
	  return;
     }
     if ((evenReadMatchStructs[0].readMatchBegin + oddReadMatchStructs[0].readMatchBegin <= maxTotAllowableMissingOnEnds)) {
	  mateUnitig1ori = evenReadMatchStructs[0].ori;
	  if (oddReadMatchStructs[0].ori == 'F')
	       mateUnitig2ori = 'R';
	  else
	       mateUnitig2ori = 'F';
	  if (mateUnitig1ori == 'F')
	       lengthAdjustment1 = evenReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment1 = - evenReadMatchStructs[0].bhg;
	  if (mateUnitig2ori == 'R')
	       lengthAdjustment2 = oddReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment2 = - oddReadMatchStructs[0].bhg;
	  insertLengthMean = mean[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]] + (lengthAdjustment1 + lengthAdjustment2);
//	  printf ("lA1 = %d; lA2 = %d\n", lengthAdjustment1, lengthAdjustment2);
		   
#ifdef KILLED111115
	  printf ("joinKUnitigsFromMates for pair %s%lld %s%lld using mean %d\n", rdPrefixHold, readNumHold-1, rdPrefixHold, readNumHold, insertLengthMean);
#endif
	  // The following to check what we're doing
	  insertLengthMeanBetweenKUnisForInsertGlobal = insertLengthMean;
	  insertLengthStdevGlobal = stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]];
	  if (joinKUnitigsFromMates (insertLengthMean, stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]]) == 1)
	       return;
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads(readNameSpace);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads(readNameSpace);
     }
     return;
}

