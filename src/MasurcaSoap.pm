package MasurcaSoap;

use POSIX qw(floor);
use MasurcaConf qw(fail);

my $SOAP_dir = "SOAP_assembly";
my $SOAP_CONF = "soap_config";

sub round { floor($_[0] + 0.5); }

sub SOAPconfig {
  my (%config) = @_;

  open(my $io, ">", $SOAP_CONF) or fail("Can't create SOAP configuration file '$SOAP_CONF': $!");
  print $io <<"EOS";
# Configuration file autogenerated by MaSuRSOAP.
max_rd_len=65536 
[LIB]
avg_ins=200
reverse_seq=0
asm_flags=1
rank=1
f=../work1/superReadSequences.fasta.all
EOS
#special case -- only one library
if(scalar(@{$config{PE_INFO}})==1){
my $lib=${$config{PE_INFO}}[0];
print $io <<"EOS";
[LIB]
avg_ins=@$lib[1]
reverse_seq=0
asm_flags=2
rank=1
map_len=63
p=../pe.cor.fa
EOS
}else{

  foreach my $lib (@{$config{PE_INFO}}) {
    my $mean = abs(@$lib[1]);
    my $prefix = @$lib[0];
print $io <<"EOS";
[LIB]
avg_ins=$mean
reverse_seq=0
asm_flags=2
rank=1
map_len=63
p=../$prefix.cor.fa
EOS
}
}

#special case -- only one library
if(scalar(@{$config{JUMP_INFO}})==1){
my $lib = ${$config{JUMP_INFO}}[0];
my $mean = abs(@$lib[1]);
my $rev_seq = @$lib[1] < 0 ? "0" : "1"; 
print $io <<"EOS";
[LIB]
avg_ins=$mean
reverse_seq=$rev_seq
asm_flags=2
rank=2
map_len=63
p=../sj.cor.clean2.fa
EOS
}else{
  my %ranks;
  foreach my $lib (@{$config{JUMP_INFO}}) {
    my $mean = abs(@$lib[1]);
    push(@{$ranks{round($mean / 2000)}}, $lib);
  }
  my @okeys = sort { $a <=> $b } keys(%ranks);
  for(my $i = 0; $i < @okeys; $i++) {
    my $rank = $i + 2;
    foreach my $lib (@{$ranks{$okeys[$i]}}) {
      my ($name, $mean, $stdev, $f1, $f2) = @$lib;
      my $map_len=$mean>0?51:35;
      my $file = "../" . $name . ".cor.clean.fa";
      my $rev_seq = $mean < 0 ? "0" : "1";
      $mean=abs($mean);
      print $io <<"EOS";
[LIB]
avg_ins=$mean
reverse_seq=$rev_seq
asm_flags=2
rank=$rank
map_len=$map_len
p=$file
EOS
     }
    }
  }
}

sub runSOAP {
  my ($out, $reads_file, %config) = @_;
  my $cmdline_jump="splitFileByPrefix.pl " . join(" ", map { $$_[0] } @{$config{JUMP_INFO}});
  my $cmdline_pe="splitFileByPrefix.pl " . join(" ", map { $$_[0] } @{$config{PE_INFO}});
 

print $out "log 'SOAPdenovo'\n";
print $out "$cmdline_pe < pe.cor.fa\n" unless(scalar(@{$config{PE_INFO}})==1); 
print $out "$cmdline_jump < sj.cor.clean2.fa\n" unless(scalar(@{$config{JUMP_INFO}})==1);

print $out <<"EOS";
mkdir -p $SOAP_dir
( cd $SOAP_dir
  [ \$KMER -le 63 ] && cmd=SOAPdenovo-63mer || cmd=SOAPdenovo-127mer
  \$cmd all -F -u -w -p $config{NUM_THREADS} -D 0 -d 0 -K \$KMER -k 35 -R -o asm -s ../$SOAP_CONF 1>../SOAPdenovo.err 2>\&1
)
[ -e "$SOAP_dir/asm.scafSeq" ] || fail SOAPdenovo failed, Check SOAPdenovo.err for problems.
EOS
  
my $reads_argument= join(" ", map { "--reads-file '$_'" } @$reads_file);
print $out <<"EOS";
log 'Gap closing'
closeGapsInScaffFastaFile.perl --split 1 --max-reads-in-memory 1000000000 -s $config{JF_SIZE} --scaffold-fasta-file  $SOAP_dir/asm.scafSeq $reads_argument --output-directory SOAP_gapclose --min-kmer-len 19 --max-kmer-len \$((\$PE_AVG_READ_LENGTH-5)) --num-threads $config{NUM_THREADS} --contig-length-for-joining \$((\$PE_AVG_READ_LENGTH-1)) --contig-length-for-fishing 200 --reduce-read-set-kmer-size 25 1>gapClose.err 2>&1
[ -e "SOAP_gapclose/genome.ctg.fasta" ] || fail Gap close failed, you can still use pre-gap close scaffold in file $SOAP_dir/asm.scafSeq. Check gapClose.err for problems.
log 'Rescaffolding'
(cd $SOAP_dir
finalFusion -K \$KMER -g asm2 -c ../SOAP_gapclose/genome.scf.fasta -D >> ../SOAPdenovo.err
  [ \$KMER -le 63 ] && cmd=SOAPdenovo-63mer || cmd=SOAPdenovo-127mer
  \$cmd map -s ../$SOAP_CONF -g asm2 -p $config{NUM_THREADS} -k 35  1>>../SOAPdenovo.err 2>\&1
  \$cmd scaff -g asm2 -p $config{NUM_THREADS} -w -u   1>>../SOAPdenovo.err 2>\&1
  restore_ns.pl asm2.contig asm2.contigPosInscaff > asm2.scafSeq2
)
[ -e "$SOAP_dir/asm2.scafSeq2" ] || fail SOAPdenovo failed, Check SOAPdenovo.err for problems.
log Assembly success. Output sequence is in SOAP_assembly/asm2.scafSeq2
EOS
}

1;
