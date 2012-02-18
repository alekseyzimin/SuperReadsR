#!/usr/bin/perl
#
# This exec takes a (set of) input read files (in fasta format) and a
# file of input k-unitigs (specified with the switch -kunitigsfile) and outputs
# the set of super-reads for these reads (in fasta format).
#
# The args are the input fasta files as well as (optionally) the directory
# where you want the work to occur. If the directory is not specified, the
# work is done in the current directory.
# If the work directory doesn't exist then it is created.
# 
# The flags are as follows:
# -l merLen : the length of the k-mer to use for the calculations (31)
# -s tableSize : the size of the table when running jellyfish (2,000,000,000)
# -t numProcessors : the number of processors to run jellyfish and create_k_unitigs (16)
# -kunitigsfile filename : a user-given k-unitigs file; otherwise we calculate
# -mean-and-stdev-by-prefix-file filename : a file giving mate info about each
#                      library. Each line is the 2-letter prefix for the reads
#                      in the library followed by its mean and stdev. This
#                      file is mandatory unless -jumplibraryreads is specified
# -mkudisr numBaseDiffs : max base diffs between overlapping k-unitigs in super-reads (0)
# -minreadsinsuperread minReads : super-reads containing fewer than numReads
#                                reads will be eliminated (2)
# -merged-unitig-data-prefix prefix : the prefix for the filenames relating to
#                      merged unitig data. We assume that the k-unitig sequence
#                      is in  'prefix'.fasta, and the map file from orig to
#                      merged k-unitigs is in 'prefix'.map.
# --stopAfter target : Stop the run after one of the following "target" names:
#               createLengthStatisticsFiles
#               createKUnitigHashTable
#               addMissingMates
#               findReadKUnitigMatches
#               createLengthStatisticsForMergedKUnitigsFiles
#               createKUnitigMaxOverlaps
#               joinKUnitigs
#               getSuperReadInsertCounts
#               createFastaSuperReadSequences
#               reduceSuperReads
#               createFinalReadPlacementFile
#               createFinalSuperReadFastaSequences
# -noclean : don't clean up the files afterwards
# -mikedebug : don't kill off intermediate results
# -jumplibraryreads : we are generating for jump-library reads; a k-unitigs
#                                 file must be specified
# -h : help 
use File::Basename;
use Cwd;
$exeDir = dirname ($0);
$pwd = cwd;
if ($exeDir !~ /^\//) {
    $exeDir = "$pwd/$exeDir"; }

$noReduce=0;
&processArgs;
$merLenMinus1 = $merLen - 1;
$maxHashFillFactor = .8;

$successFile = "$workingDirectory/superReads.success";
unlink ($successFile) if (-e $successFile);
# The following is set to 1 when the first success file for a step is missing
$mustRun = 0;
if (! -d $workingDirectory) {
    $cmd = "mkdir $workingDirectory";
    print "$cmd\n"; system ($cmd); }
# We now require that a k-unitigs file was passed on the command line
if ($kUnitigsFile !~ /^\//) {
    $kUnitigsFile = "$pwd/$kUnitigsFile"; }
$jellyfishKUnitigDataPrefix = "$workingDirectory/organismMerCountsForKUnitigs";
$jellyfishKUnitigHashFile = $jellyfishKUnitigDataPrefix . "_0";
$kUnitigLengthsFile = "$workingDirectory/kUnitigLengths.txt";
# The following stores the actual number of k-unitigs
$numKUnitigsFile = "$workingDirectory/numKUnitigs.txt";
# The following stores the largest k-unitig number (+1)
$maxKUnitigNumberFile = "$workingDirectory/maxKUnitigNumber.txt";
$totBasesInKUnitigsFile = "$workingDirectory/totBasesInKUnitigs.txt";
if ($mergedUnitigDataPrefix) {
    $mergedUnitigInputKUnitigsFile = $mergedUnitigDataPrefix . ".fasta";
    $mergedUnitigInputKUnitigMappingFile = $mergedUnitigDataPrefix . ".map";
    &runCommandAndExitIfBad ("", $mergedUnitigInputKUnitigsFile, 1, "mergedKUnitigFastaFileExists");
    &runCommandAndExitIfBad ("", $mergedUnitigInputKUnitigMappingFile, 1, "mergedKUnitigMapFileExists");
    $mergedKUnitigLengthsFile = "$workingDirectory/mergedKUnitigs.kUnitigLengths.txt";
    $mergedNumKUnitigsFile = "$workingDirectory/mergedKUnitigs.numKUnitigs.txt";
    $mergedMaxKUnitigNumberFile = "$workingDirectory/mergedKUnitigs.maxKUnitigNumber.txt";
    $mergedTotBasesInKUnitigsFile = "$workingDirectory/mergedKUnitigs.totBasesInKUnitigs.txt"; }
else {
    $mergedUnitigInputKUnitigsFile = $kUnitigsFile;
    $mergedKUnitigLengthsFile = $kUnitigLengthsFile;
    $mergedNumKUnitigsFile = $numKUnitigsFile;
    $mergedMaxKUnitigNumberFile = $maxKUnitigNumberFile;
    $mergedTotBasesInKUnitigsFile = $totBasesInKUnitigsFile; }
    
$totReadFile = "$workingDirectory/inputReads.fasta";
if ($#fastaFiles == 0) {
    $totReadFile = $fastaFiles[0]; }
$readsAfterAddingMissingMates = "$workingDirectory/inputreads.fa";
$numReadsFile = "$workingDirectory/numReads.txt";
$prefixForOverlapsBetweenKUnitigs = "$workingDirectory/overlap";
$kUnitigOverlapsFile = "${prefixForOverlapsBetweenKUnitigs}.overlaps";

$joinerOutputPrefix = "$workingDirectory/readPositionsInSuperReads";
$myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.txt";
$myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly";
$sequenceCreationErrorFile = "$workingDirectory/createFastaSuperReadSequences.errors.txt";
# $myProgOutput2 = "$workingDirectory/readPlacementsInSuperReads.postMateMerge.read.superRead.offset.ori.txt";
# $myProgOutput3 = "$workingDirectory/superReadCounts.count.superRead.txt";
$finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.fasta";
$finalReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt";
# $reducedReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.reduced.read.superRead.offset.ori.txt";

if ($jumpLibraryReads) {
    $myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.jumpLibrary.txt";
    $myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly.jumpLibrary";
    $myProgOutput8 = "$workingDirectory/superReadSequences.jumpLibrary.bothSidesExtended.fasta";
    $myProgOutput12 = "$workingDirectory/readPlacementsInSuperReads.forJumpLibraryWBothSidesExtended.read.superRead.offset.ori.txt";
    $finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.jumpLibrary.fasta";
}

# goto startHere;


&cleanUpFailureDirectories;

# In addition to obvious output file, this also generates the files
# numKUnitigs.txt, maxKUnitigNumber.txt, and totBasesInKUnitigs.txt in
# $workingDirectory
$cmd = "cat $kUnitigsFile | $exeDir/getLengthStatisticsForKUnitigsFile.perl $workingDirectory > $kUnitigLengthsFile";
&runCommandAndExitIfBad ($cmd, $kUnitigLengthsFile, 1, "createLengthStatisticsFiles", $totBasesInKUnitigsFile, $numKUnitigsFile, $maxKUnitigNumberFile, $kUnitigLengthsFile);

$minSizeNeededForTable = &reportMinJellyfishTableSizeForKUnitigs;
 redoKUnitigsJellyfish:
    $cmd = "jellyfish count -m $merLen -r -o $jellyfishKUnitigDataPrefix -c 6 -p 126 --both-strands -s $minSizeNeededForTable -t $numProcessors $kUnitigsFile";
&runCommandAndExitIfBad ($cmd, $jellyfishKUnitigHashFile, 1, "createKUnitigHashTable", "$workingDirectory/organismMerCountsForKUnitigs_0");

$tableResizeFactor = &returnTableResizeAmount ($jellyfishKUnitigDataPrefix, $jellyfishKUnitigHashFile);
if ($tableResizeFactor > 1) {
    $tableSize *= 2;
    print "Resizing the table to $tableSize for the k-unitig jellyfish run\n";
    goto redoKUnitigsJellyfish; }

$cmd = "cat $totReadFile | $exeDir/add_missing_mates.pl >  $readsAfterAddingMissingMates";
&runCommandAndExitIfBad ($cmd, $readsAfterAddingMissingMates, 1, "addMissingMates", $readsAfterAddingMissingMates);

$cmd = "$exeDir/findMatchesBetweenKUnitigsAndReads $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1_1prefix $kUnitigsFile $maxKUnitigNumberFile  $readsAfterAddingMissingMates";
&runCommandAndExitIfBad ($cmd, $myProgOutput1_1prefix . "*", 1, "findReadKUnitigMatches", "$workingDirectory/newTestOutput.nucmerLinesOnly_*");
if (! $mikedebug) { &killFiles ($jellyfishKUnitigHashFile, $readsAfterAddingMissingMates); }

if ($jumpLibraryReads) {
    goto jumpLibraryCalculations; }

# In addition to obvious output file, this also generates the files
# mergedKUnitigs.numKUnitigs.txt, mergedKUnitigs.maxKUnitigNumber.txt, and mergedKUnitigs.totBasesInKUnitigs.txt in
# $workingDirectory
if ($mergedUnitigDataPrefix) {
    $cmd = "cat $mergedUnitigInputKUnitigsFile | $exeDir/getLengthStatisticsForKUnitigsFile.perl -output-prefix mergedKUnitigs $workingDirectory > $mergedKUnitigLengthsFile";
    &runCommandAndExitIfBad ($cmd, $mergedKUnitigLengthsFile, 1, "createLengthStatisticsForMergedKUnitigsFiles", $mergedTotBasesInKUnitigsFile, $mergedNumKUnitigsFile, $mergedMaxKUnitigNumberFile, $mergedKUnitigLengthsFile); }
else { # The following is so we stop here if we are using --stopAfter createLengthStatisticsForMergedKUnitigsFiles but don't use the merged k-unitigs
    &runCommandAndExitIfBad ("", "", 0, "createLengthStatisticsForMergedKUnitigsFiles"); }

open (FILE, $mergedMaxKUnitigNumberFile); $maxKUnitigNumber = <FILE>; chomp ($maxKUnitigNumber); close (FILE);
$cmd = "$exeDir/createKUnitigMaxOverlaps $mergedUnitigInputKUnitigsFile -kmervalue $merLen -largest-kunitig-number ".(int($maxKUnitigNumber)+1)." $prefixForOverlapsBetweenKUnitigs";
&runCommandAndExitIfBad($cmd, $kUnitigOverlapsFile, 1, "createKUnitigMaxOverlaps", $kUnitigOverlapsFile, "$workingDirectory/overlap.coords");

# Do the shooting method here
if ($mergedUnitigDataPrefix) {
    $mergedUnitigDataFileStr = "--kunitigs-translation-file $mergedUnitigInputKUnitigMappingFile"; }
$cmd = "$exeDir/joinKUnitigs_v3 --mean-and-stdev-by-prefix-file $meanAndStdevByPrefixFile --unitig-lengths-file $mergedKUnitigLengthsFile --num-kunitigs-file $mergedMaxKUnitigNumberFile --overlaps-file $kUnitigOverlapsFile --min-overlap-length $merLenMinus1 --prefix $joinerOutputPrefix $mergedUnitigDataFileStr --num-file-names $numProcessors $myProgOutput1_1prefix";
&runCommandAndExitIfBad ($cmd, $joinerOutputPrefix . "*", 1, "joinKUnitigs", "$workingDirectory/readPositionsInSuperReads_*");

$cmd= "$exeDir/getSuperReadInsertCountsFromReadPlacementFileTwoPasses -n `cat $numKUnitigsFile` -o $workingDirectory/superReadCounts.all ${joinerOutputPrefix}_*";
&runCommandAndExitIfBad ($cmd,"$workingDirectory/superReadCounts.all", 1, "getSuperReadInsertCounts", "$workingDirectory/superReadCounts.all");

if($noReduce==0) {
    if ($mergedUnitigDataPrefix) {
	$mergedUnitigDataFileStr = "-maxunitignumberfile $mergedMaxKUnitigNumberFile"; }
    $cmd = "cat $workingDirectory/superReadCounts.all | $exeDir/createFastaSuperReadSequences $workingDirectory /dev/fd/0 -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread $minReadsInSuperRead $mergedUnitigDataFileStr -good-sr-filename $workingDirectory/superReadNames.txt -kunitigsfile $mergedUnitigInputKUnitigsFile 2> $sequenceCreationErrorFile | tee $finalSuperReadSequenceFile.all | perl -ane 'BEGIN{my \$seq_length=0}{if(\$F[0] =~ /^>/){if(\$seq_length>0){print \$seq_length,\"\\n\";} print substr(\$F[0],1),\" \";\$seq_length=0;}else{\$seq_length+=length(\$F[0]);}}END{if(\$seq_length>0){print \$seq_length,\"\\n\";}}' | sort -nrk2,2 -S 40% -T ./ > $workingDirectory/sr_sizes.tmp";
    &runCommandAndExitIfBad ($cmd,"$workingDirectory/sr_sizes.tmp", 1, "createFastaSuperReadSequences", "$workingDirectory/superReadSequences.fasta.all", "$workingDirectory/superReadNames.txt", "$workingDirectory/sr_sizes.tmp");

    $cmd = "cat $workingDirectory/sr_sizes.tmp| $exeDir/reduce_sr $maxKUnitigNumber  > $workingDirectory/reduce.tmp";
    &runCommandAndExitIfBad ($cmd,"$workingDirectory/reduce.tmp", 1, "reduceSuperReads", "$workingDirectory/reduce.tmp", "$workingDirectory/createFastaSuperReadSequences.errors.txt");

    $cmd = "cat ${joinerOutputPrefix}_* | $exeDir/eliminateBadSuperReadsUsingList /dev/fd/0 $workingDirectory/superReadNames.txt | perl -e '{open(FILE,\$ARGV[0]);while(\$line=<FILE>){chomp(\$line);\@F=split(\" \",\$line);\$sr{\$F[0]}=\$F[1]} while(\$line=<STDIN>){chomp(\$line);\@l=split(\" \",\$line);if(defined(\$sr{\$l[1]})){print \"\$l[0] \",\$sr{\$l[1]},\" \$l[2] \$l[3]\\n\";}else{print \"\$line\\n\";}}}' $workingDirectory/reduce.tmp >  $finalReadPlacementFile"; 
    &runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1, "createFinalReadPlacementFile", "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt");

    $cmd = "awk '{print \$1}' $workingDirectory/reduce.tmp > $workingDirectory/sr_to_eliminate.tmp; $exeDir/extractreads_not.pl $workingDirectory/sr_to_eliminate.tmp $finalSuperReadSequenceFile.all 1 > $finalSuperReadSequenceFile";
    &runCommandAndExitIfBad ($cmd,$finalSuperReadSequenceFile, 1, "createFinalSuperReadFastaSequences", "$workingDirectory/superReadSequences.fasta", "$workingDirectory/sr_to_eliminate.tmp"); }
else {
    $cmd = "cat $workingDirectory/superReadCounts.all | $exeDir/createFastaSuperReadSequences $workingDirectory /dev/fd/0 -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread $minReadsInSuperRead -good-sr-filename $workingDirectory/superReadNames.txt -kunitigsfile $mergedUnitigInputKUnitigsFile 2> $sequenceCreationErrorFile > $finalSuperReadSequenceFile";
    &runCommandAndExitIfBad ($cmd,$finalSuperReadSequenceFile, 1, "createFastaSuperReadSequences", "$workingDirectory/superReadSequences.fasta", "$workingDirectory/superReadNames.txt");

    $cmd = "cat ${joinerOutputPrefix}_* | $exeDir/eliminateBadSuperReadsUsingList /dev/fd/0 $workingDirectory/superReadNames.txt > $finalReadPlacementFile";
    &runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1, "createFinalReadPlacementFile", "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt", "$workingDirectory/createFastaSuperReadSequences.errors.txt");
}

$cmd = "touch $successFile";
system ($cmd);

exit (0);

 jumpLibraryCalculations:
# Must run createFastaSuperReadSequences here
    
$cmd = "$exeDir/outputSuperReadSeqForJumpLibrary.perl $myProgOutput8 $myProgOutput12 > $finalSuperReadSequenceFile";
&runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1, "createFastaSuperReadSequences");

$cmd = "touch $successFile";
system ($cmd);

sub processArgs
{
    $tableSize = 2000000000;
    $merLen = 31;
    $numProcessors = 16;
    $minReadsInSuperRead = 2;
    $seqDiffMax = 0;
    $help = 0;
    if ($#ARGV < 0) {
	$help = 1; }
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "-l") {
	    ++$i;
	    $merLen = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-s") {
	    ++$i;
	    $tableSize = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-t") {
	    ++$i;
	    $numProcessors = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-kunitigsfile") {
	    ++$i;
	    $kUnitigsFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mkudisr") {
	    ++$i;
	    $seqDiffMax = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-minreadsinsuperread") {
	    ++$i;
	    $minReadsInSuperRead = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-merged-unitig-data-prefix") {
	    ++$i;
	    $mergedUnitigDataPrefix = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mean-and-stdev-by-prefix-file") {
	    ++$i;
	    $meanAndStdevByPrefixFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "--stopAfter") {
	    ++$i;
	    $stopAfter = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-jumplibraryreads") {
	    $jumpLibraryReads = 1;
	    next; }
	elsif ($ARGV[$i] eq "-noreduce") {
            $noReduce = 1;
            next; }
	elsif ($ARGV[$i] eq "-h") {
	    $help = 1;
	    next; }
	elsif ($ARGV[$i] eq "-mikedebug") {
	    $mikedebug = 1;
	    next; }
	elsif (-f $ARGV[$i]) {
	    push (@fastaFiles, $ARGV[$i]);
	    next; }
	else {
	    if (! $workingDirectory) {
		$workingDirectory = $ARGV[$i]; }
	    else {
		print STDERR "Working directory was specified as \"$workingDirectory\", now specified as ",$ARGV[$i],".\nThis message could also occur if an input read file doesn't exist.\n";
		$help = 1; } } }
    if ($#fastaFiles < 0) {
	$help = 1; }
    if (! $kUnitigsFile) {
	print STDERR "A k-unitigs file must be supplied\n";
	$help = 1; }
    if ((! $meanAndStdevByPrefixFile) && (! $jumpLibraryReads)) {
	print STDERR "A file specifying mean and stdev of each library must be given (using the \n  '-mean-and-stdev-by-prefix-file' switch) if -jumplibraryreads is not set.\n";
	$help = 1; }
    if ($help) {
	&giveUsageAndExit; }
}

sub giveUsageAndExit
{
    open (FILE, $0);
    $line = <FILE>;
    while ($line = <FILE>) {
	chomp ($line);
	last unless ($line =~ /^\#/);
	substr ($line, 0, 2) = "";
	print "$line\n"; }
    exit (0);
}

sub returnTableResizeAmount
{
    my ($dataPrefix, $hashFile) = @_;
    
    for ($i=1; 1; $i++) {
	$fn = $dataPrefix . "_$i";
	last unless (-e $fn);
	unlink ($fn); }
    unlink ($hashFile) if ($i > 1);
    return ($i);
}

sub getNumLines
{
    my ($fn) = @_;
    my ($cmd, $line, @flds);
    
    $cmd = "wc -l $fn |";
    open (CRAZYFILE, $cmd);
    $line = <CRAZYFILE>;
    close (CRAZYFILE);
    @flds = split (" ", $line);
    return ($flds[0]);
}

sub reportMinJellyfishTableSizeForKUnitigs
{
    my ($numKMersInKUnitigs, $numKUnitigs, $minSizeNeededForTable);
    
    open (FILE, $totBasesInKUnitigsFile);
    $numKMersInKUnitigs = <FILE>;
    close (FILE);
    chomp ($numKMersInKUnitigs);
    
    open (FILE, $numKUnitigsFile);
    $numKUnitigs = <FILE>;
    chomp ($numKUnitigs);
    close (FILE);
    $numKMersInKUnitigs -= ($numKUnitigs * ($merLenMinus1));
    $minSizeNeededForTable = int ($numKMersInKUnitigs/$maxHashFillFactor + .5);
    return ($minSizeNeededForTable);
}

sub killFiles
{
    my (@filesToKill) = @_;
    my ($file);
    
    for (@filesToKill) {
	$file = $_;
	if (-e $file) {
	    unlink ($file); } }
}

# If localCmd is set, it captures the return code and makes
# sure it ran properly. Otherwise it exits.
# If one sets the fileName then we assume it must exist
# With minSize one can set the minimum output file size
sub runCommandAndExitIfBad
{
    my ($localCmd, $fileName, $minSize, $stepName, @filesCreated) = @_;
    my ($retCode, $exitValue, $sz);
    my ($totSize, $tempFilename, $cmd);
    my ($successFile, $failDir);
    
#    sleep (5); # For testing only
    $failDir = $workingDirectory . "/" . $stepName . ".Failed";
    if (-e $failDir) {
	$cmd = "rm -rf $failDir"; print "$cmd\n"; system ($cmd); }
    $successFile = $workingDirectory . "/" . $stepName . ".success";
    if (-e $successFile) {
	if ($mustRun) {
	    unlink ($successFile); }
	else {
	    print STDERR "Step $stepName already completed. Continuing.\n";
	    goto successfulRun; } }
    else {
	$mustRun = 1; }
    
    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
        system ("time $localCmd");
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
	    goto failingRun; }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
	    goto failingRun; }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
		$retCode = $exitValue;
		goto failingRun; }
        }
    }
    goto successfulRun unless ($fileName =~ /\S/);
    if ($fileName =~ /\*/) {
	goto multipleFiles; }
    if (! -e $fileName) {
        print STDERR "Output file \"$fileName\" doesn't exist. Bye!\n";
	$retCode = 1;
	goto failingRun; }
    $sz = -s $fileName;
    if ($sz < $minSize) {
        print STDERR "Output file \"$fileName\" is of size $sz, must be at least of size $minSize. Bye!\n";
	$retCode = 1;
	goto failingRun; }
    goto successfulRun;
  multipleFiles:
    $cmd = "ls $fileName |";
    $totSize = 0;
    open (CMD, $cmd);
    while ($tempFilename = <CMD>) {
	chomp ($tempFilename);
	$sz = -s $tempFilename;
	$totSize += $sz; }
    close (CMD);
    if ($totSize < $minSize) {
	print STDERR "The combined output files from \"$fileName\" have a total size of $totSize, must be at least of size $minSize. Bye!\n";
	$retCode = 1;
	goto failingRun; }
  successfulRun:
    if (! -e $successFile) {
	$cmd = "touch $successFile"; print "$cmd\n"; system ($cmd); }
    if ($stopAfter eq $stepName) {
	print STDERR "Stopping after step ${stepName}. Bye!\n";
	exit (0); }
    return;
    
  failingRun:
    $outdir = "$workingDirectory/${stepName}.Failed";
    mkdir ($outdir);
    for (@filesCreated) {
	$cmd = "mv $_ $outdir";
	print "$cmd\n"; system ($cmd); }
    exit ($retCode);
}

sub cleanUpFailureDirectories
{
    my ($tfile, $totFile, $cmd);

    if (! -e $workingDirectory) {
	return; }
    opendir (DIR, $workingDirectory);
    while ($tfile = readdir (DIR)) {
	next unless ($tfile =~ /\.Failed$/);
	$totFile = "$workingDirectory/$tfile";
	next unless (-d $totFile);
	$cmd = "rm -rf $totFile";
	print "$cmd\n"; system ($cmd); }
    closedir (DIR);
}


    
