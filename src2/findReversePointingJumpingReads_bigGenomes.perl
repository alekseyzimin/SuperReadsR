#!/usr/bin/env perl
use Cwd;
use File::Basename;
use File::Copy;
$exeDir = dirname ($0);
$localReadsFile = "localReadsFile";
$joiningEndPairsOrig = "readsToJoinReads.fasta";
$meanAndStdevFn = "meanAndStdevByPrefix.sj.txt";
$meanAndStdevFile = returnAbsolutePath ($meanAndStdevFn);
&processArgs;
$globalIsGood = 1;
$isGood = returnIfAllFilesExist (@readsFiles, @readsForKUnitigsFiles);
if (! $isGood) {
    $globalIsGood = 0; }
if ($#readsFiles < 0) {
    print STDERR "You must enter some reads files.\n";
    $globalIsGood = 0; }
if ($#readsForKUnitigsFiles < 0) {
    if ($dirForKUnitigs !~ /\S/) {
	print STDERR "You must enter some files containing reads to generate k-unitigs or a directory which has k-unitigs files we can use.\n";
	$globalIsGood = 0; }
    else {
	if (! -d $dirForKUnitigs) {
	    print STDERR "The directory specified to get the k-unitigs, \'$dirForKUnitigs\' doesn't exist or isn't a directory.\n";
	    $globalIsGood = 0; }
    }
}
if (! -e $meanAndStdevFn) {
    print STDERR "\'$meanAndStdevFn\' doesn't exist.\n";
    $globalIsGood = 0; }


for ($i=0; $i<=$#readsFiles; ++$i) {
    $readsFiles[$i] = returnAbsolutePath ($readsFiles[$i]); }
for ($i=0; $i<=$#readsForKUnitigsFiles; ++$i) {
    $readsForKUnitigsFiles[$i] = returnAbsolutePath ($readsForKUnitigsFiles[$i]); }
if (! -e $dirToChangeTo) {
    $cmd = "mkdir $dirToChangeTo";
    runCommandAndExitIfBad ($cmd); }

if ($dirToChangeTo) {
    chdir ($dirToChangeTo); }

$cmd = "zcat -f @readsFiles > $joiningEndPairsOrig";
runCommandAndExitIfBad ($cmd);

$tempVal = $minKMerLen - 1;
$joiningEndPairs = $joiningEndPairsOrig . ".$tempVal";
$cmd = "ln -s $joiningEndPairsOrig $joiningEndPairs";
runCommandAndExitIfBad ($cmd);

$totInputSize = getReadFileSize (@readsForKUnitigsFiles);
if (! $jellyfishHashSize) {
    $jellyfishHashSize = $totInputSize; }
$readPrefix = &getReadPrefix ($joiningEndPairs);
$joiningEndPairNamesFile = "readsToExclude.txt";
unlink ($joiningEndPairNamesFile);

&runMainLoop;

sub returnIfAllFilesExist
{
    my (@files) = @_;
    my ($isGood);
    $isGood = 1;

    for (@files) {
	$file = $_;
	if (! -e $file) {
	    print STDERR "File \'$file\' doesn't exist\n";
	    $isGood = 0; }
	elsif (-d $file) {
	    print STDERR "File \'$file\' is a directory\n";
	    $isGood = 0; }
    }

    return ($isGood);
}

sub runMainLoop
{
    my ($k, $suffix, $cmd, $directoryPassed);
    my ($fn1, $fn2);
    my ($localCode, $tempLine);
    my ($outFn, $sz, $sz2, $minContinuation);

    $directoryPassed=0; 
    # Here's the main loop (k-mer values going up)
    for ($k = $minKMerLen; $k<=$maxKMerLen; $k+=$kmerStepSize) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	$tempKUnitigsFile = "k_unitigs_${suffix}.fa";
#	$minContinuation = int ($k/2);
	$minContinuation = $k-1;
	if (-e $tempKUnitigsFile) {
	    goto afterKUnitigCreation; }
	if ($dirForKUnitigs =~ /\S/) {
	    $tfile = "$dirForKUnitigs/$tempKUnitigsFile";
	    if ((-e $tfile) && (-s $tfile > 0)) {
		$cmd = "ln -s $tfile $tempKUnitigsFile";
		runCommandAndExitIfBad ($cmd);
		goto afterKUnitigCreation; }
	}
	$cmd = "$exeDir/create_k_unitigs_large_k2 -c $minContinuation -t $numThreads -m $k -n $jellyfishHashSize -l $k @readsForKUnitigsFiles  |  grep --text -v '^>' | perl -ane '{\$seq=\$F[0]; \$F[0]=~tr/ACTGacgt/TGACtgac/;\$revseq=reverse(\$F[0]); \$h{(\$seq ge \$revseq)?\$seq:\$revseq}=1;}END{\$n=0;foreach \$k(keys \%h){print \">\",\$n++,\" length:\",length(\$k),\"\\n\$k\\n\"}}' > $tempKUnitigsFile";
	if (runCommandAndReturnIfBad ($cmd)) {
	    last; }
	if (-s $tempKUnitigsFile == 0) {
	    print STDERR "create_k_unitigs_large_k2 for k = $k failed.\n";
	    last; }

      afterKUnitigCreation:
	$cmd = "\\rm -rf out.$suffix"; system ($cmd);
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);

	if (! -e $meanAndStdevFn) {
	    $cmd = "cp $meanAndStdevFile $meanAndStdevFn";
	    runCommandAndExitIfBad ($cmd);
	}

	if ($k == $minKMerLen) {
	    $lastKMerSize = $k-1; }
	else {
	    $lastKMerSize = $k - $kmerStepSize; }
	$inputEndPairs = $joiningEndPairsOrig . ".$lastKMerSize";
	$outputEndPairs = $joiningEndPairsOrig . ".$k";

#	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file $meanAndStdevFn -num-stdevs-allowed $numStdevsAllowed -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file $meanAndStdevFn -num-stdevs-allowed $numStdevsAllowed -closegaps -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -join-aggressive $joinAggressive -maxnodes $maxNodes -mkudisr 0 --stopAfter joinKUnitigs work_${suffix} $inputEndPairs 1>>out.$suffix 2>>out.$suffix";
	if (runCommandAndReturnIfBad ($cmd)) {
	    last; }

	$readPositionFile = "work_$suffix/readPositionsInSuperReads";
	$readNamesForNextPassFile = "readNamesForNextPass.${k}.txt";
	$cmd = "$exeDir/extractJoinableAndNextPassReadsFromJoinKUnitigs.perl $readPositionFile $joiningEndPairNamesFile $readNamesForNextPassFile";
	runCommandAndExitIfBad ($cmd);

	undef %readsForNextPass;
	$sz = -s $readNamesForNextPassFile;
	last if ($sz == 0);
	open (FILE, $readNamesForNextPassFile);
	while ($line = <FILE>) {
	    chomp ($line);
	    $readsForNextPass{$line} = 1; }
	close (FILE);
	unlink ($readNamesForNextPassFile);

	open (FILE, $inputEndPairs);
	open (OUTFILE, ">$outputEndPairs");
	while ($line = <FILE>) {
	    $line2 = <FILE>;
	    chomp ($line);
	    ($rd) = ($line =~ /^.(\S+)/);
	    if ($readsForNextPass{$rd}) {
		print OUTFILE "$line\n$line2"; }
	}
	close (FILE); close (OUTFILE);
	if (! $keepKUnitigs) {
	    unlink ($tempKUnitigsFile); }
	if (! $keepDirectories) {
	    $cmd = "\\rm -rf work_$suffix"; system ($cmd); }
    }
    return;
    
    $passingKMerFile = "passingKMer.txt";
    open (FILE, ">$passingKMerFile");
    if ($directoryPassed) { 
        print FILE "$k\n"; }
    else {
	print FILE "11\n"; }
    close(FILE);
    $passingReadsFile = "passingReadsFile.txt";
    $passingReadsFileHold = "passingReadsFile.orig.txt";
    $passingReadsFileSize = -s $passingReadsFile;
    if ($passingReadsFileSize > 0) {
	undef %isPassing;
	move ($passingReadsFile, $passingReadsFileHold);
	open (FILE, $passingReadsFileHold);
	while ($line = <FILE>) {
	    chomp ($line);
	    $isPassing{$line} = 1; }
	close (FILE);
	open (OUTFILE, ">$passingReadsFile");
	open (FILE, $joiningEndPairs);
	while ($line = <FILE>) {
	    chomp ($line);
	    @flds = split (" ", $line);
	    ($name) = ($flds[0] =~ /^.(.+)$/);
	    if ($isPassing{$name}) {
		print OUTFILE "$flds[1]\n"; }
	}
	close (FILE);
    }
}

sub getReadFileSize
{
    my (@files) = @_;
    my ($totSize, $file);
    $totSize = 0;
    for (@files) {
	$file = $_;
	$totSize += (-s $file); }
    return ($totSize);
}

sub runCommandAndExitIfBad
{
    my ($localCmd) = @_;
    my ($retCode, $exitValue);
    
    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
	system ($localCmd);
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
  successfulRun:
    return (0);
    
  failingRun:
    exit ($retCode);
}

sub runCommandAndReturnIfBad
{
    my ($localCmd) = @_;
    my ($retCode, $exitValue);
    
    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
	$localCmd = "time $localCmd"; # For debugging 11/23/12
	system ($localCmd);
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
	    goto failingRun2; }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
	    goto failingRun2; }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
		$retCode = $exitValue;
		goto failingRun2; }
        }
    }
  successfulRun2:
    return (0);
    
  failingRun2:
    return ($retCode);
}

sub processArgs
{
    my ($arg, $cmd, @kmerLens, $i);
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 65;
    $minKMerLen = 17;
    $numThreads = 1;
    $maxNodes = 2000;
    $numStdevsAllowed = 5;
    $joinAggressive = 1;
    $keepDirectories = 0;
    $kmerStepSize = 1;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if (($arg eq "-t") || ($arg eq "--num-threads")) {
	    ++$i;
	    $numThreads = $ARGV[$i];
	    next; }
	if ($arg eq "--num-stdevs-allowed") {
	    ++$i;
	    $numStdevsAllowed = $ARGV[$i];
	    next; }
	if ($arg eq "--kmer-step-size") {
	    ++$i;
	    $kmerStepSize = $ARGV[$i];
	    next; }
	if ($arg eq "--jellyfish-hash-size") {
	    ++$i;
	    $jellyfishHashSize = $ARGV[$i];
	    next; }
#        if ($ARGV[$i] eq "--join-aggressive") {
#            ++$i;
#            $joinAggressive = $ARGV[$i];
#            next; }
	if ($arg eq "--min-kmer-len") {
	    ++$i;
	    $minKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--keep-kunitigs") {
	    $keepKUnitigs = 1; }
	if ($arg eq "--keep-directories") {
	    $keepDirectories = 1; }
	if ($arg eq "--max-kmer-len") {
	    ++$i;
	    $maxKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--reads-file") {
	    ++$i;
	    push (@readsFiles, $ARGV[$i]);
	    next; }
	if ($arg eq "--reads-for-kunitigs-file") {
	    ++$i;
	    push (@readsForKUnitigsFiles, $ARGV[$i]);
	    next; }
	if ($arg eq "--dir-for-kunitigs") {
	    ++$i;
	    $dirForKUnitigs = $ARGV[$i];
	    $dirForKUnitigs = returnAbsolutePath ($dirForKUnitigs);
	    next; }
#	if ($arg eq "--kunitig-continuation-number") {
#	    ++$i;
#	    $kUnitigContinuationNumber = $ARGV[$i];
#	    next; }
	if ($arg eq "--maxnodes") {
	    ++$i;
	    $maxNodes = $ARGV[$i];
	    next; }
	if ($arg eq "--dir-to-change-to") {
	    ++$i;
	    $dirToChangeTo = $ARGV[$i];
	    next; }
	if (-f $arg) {
	    push (@readsFiles, $arg);
	    next; }
	if ($arg =~ /^\d+/) {
	    push (@kmerLens, $arg);
	    next; }
    }
    @kmerLens = sort byNum @kmerLens;
    if ($#kmerLens >= 0) {
	$maxKMerLen = $kmerLens[0]; }
    if ($#kmerLens >= 1) {
	$minKMerLen = $kmerLens[1]; }
    if ($#kmerLens >= 2) {
	$kUnitigContinuationNumber = $kmerLens[2]; }
}

sub byNum
{
    return ($b <=> $a);
}

sub determineIfMultipleJoinsAreRun
{
    my ($cmd, $tline, @flds);
    $cmd = "wc -l $joiningEndPairs";
    open (FILE, "$cmd |"); $tline = <FILE>; close (FILE);
    @flds = split (" ", $tline);
    if ($flds[0] > 4) {
	return (1); }
    else {
	return (0); }
}

sub getReadPrefix
{
    my ($inputFile) = @_;
    my ($tempLine, $prefix);

    open (FILE, $inputFile);
    $tempLine = <FILE>;
    close (FILE);
    ($prefix) = ($tempLine =~ /^.(..)/); # The first skips the '>' at the beginning of the line
    return ($prefix);
}

sub returnAbsolutePath
{
    my ($file) = @_;
    if (! $cwd) {
	$cwd = cwd; }
    if ($file !~ /^\//) {
        $file = "$cwd/$file"; }
    return ($file);
}
