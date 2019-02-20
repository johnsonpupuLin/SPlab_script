#!/usr/bin/perl
use Getopt::Long;
$|=1;
$help_desk='
================================= proTRAC ====================================
VERSION: 2.0.5                                   LAST MODIFIED: 15. April 2015

Please cite:
Rosenkranz D, Zischler H. proTRAC - a software for probabilistic piRNA cluster
detection, visualization and analysis. 2012. BMC Bioinformatics 13:5.

and (for proTRAC 2.0 and later):
Rosenkranz D, Rudloff S, Bastuck K, Ketting RF, Zischler H. Tupaia small RNAs
provide insights into function and evolution of RNAi-based transposon defense
in mammals. 2015. RNA 21(5):1-12.

Contact:
David Rosenkranz
Institute of Anthropology, small RNA group
Johannes Gutenberg University Mainz
email: rosenkrd@uni-mainz.de

You can find the latest proTRAC version at:
http://sourceforge.net/projects/protrac/files
http://www.smallRNAgroup-mainz.de/software
==============================================================================
                              +++ OPTIONS +++
 ftp  = floating point number
 int  = integer
 0..1 = floating point number between 0 and 1.

 -genome [file]          Name of the file that contains the genomic sequence
                         that was used for mapping the sequence reads.
 -map [file]             Name of the file that contains mapped reads in ELAND3
                         format. Use SeqMap with option /output_all_matches or
                         sRNAmapper to create an appropriate file.
 -help                   Show this information.
 -h                      Show this information.
 -repeatmasker [file]    Name of the file that contains the RepeatMasker
                         annotation. Make sure that the names for the
                         chromosomes/scaffolds are identical in your Repeat-
                         Masker and genome file.
 -geneset [file]         Name of the file that contains gene annotation (GTF-
                         file from Ensembl database). Make sure that the names
                         for the chromosomes/scaffolds are identical in your
                         GTF- and genome file.
 -swsize [int]           Size of the sliding window (default=5000)
 -swincr [int]           Increament of the sliding window (default=1000)
 -nh                     Normalize hit count by total number of genomic hits
                         for that sequence.
 -nr                     Normalize hit count by number of reads for that
                         sequence. Read number must be the FASTA header in
                         the file that was used for mapping.
 -rpm                    Normalize the hit count with the total number of
                         mapped sequence reads (reads per million).
                         This will make the values comparable accross
                         different piRNA libraries and or organisms.
 -dens [fpt]             Define an absolute minimum number of (normalized) 
                         read counts per kb.
 -pdens [0..1]           Define a p-value for minimum number of (normalized)
                         read counts per kb. A p-value of 0.01 means that the
                         the (normalized) read counts in a sliding window must
                         belong to the top 1% of all sliding windows.
 -est [int]              Use that option together with -pdens. Estimate the 
                         required minimum number of (normalized) read counts
                         in a sliding windows based on n random 1kb sliding
                         windows (faster). Without that option proTRAC will
                         scan the map file and calculate the required minimum
                         number of (normalized) read counts in a sliding
                         window based on the observed distribution. We
                         recommend not to use this option.
 -pisize [0..1]          Fraction of (normalized) read counts that have
                         the typical piRNA size (default=0.75).
 -pimin [int]            Define the minimum length of a piRNA (default=24).
 -pimax [int]            Define the maximum length of a piRNA (default=30).
 -1Tor10A [0..1]         Fraction of (normalized) read counts that have 1T
                         (1U) or 10A (default=0.75).
 -1Tand10A [0..1]        If the fraction of (normalized) read counts with 1T
                         (1U) OR 10A is below the value defined by -1Tor10A,
                         accept the sliding window if BOTH the 1T (1U) and the
                         10A fraction reach this value (default=0.5).
 -distr [int-int]        To avoid false positive piRNA cluster annotation of
                         loci with one or few mapped sequences represented by
                         exceptionally many reads. You can e.g. type
                         \'-distr 10-75\' which means that the TOP 10% of
                         mapped sequences account for max. 75% of the piRNA
                         clusters (normalized) read counts. Otherwise the
                         locus is rejected.
 -clsize [int]           Set the minimum size for a piRNA cluster (default=
                         5000).
 -clhits [int]           Minimum number of sequence hit loci per piRNA cluster
                         (default=0).
 -clhitsn [ftp]          Minimum number of normalized sequence read counts per
                         piRNA cluster (default=0).
 -clstrand [0..1]        Fraction of (normalized) read counts that map to the
                         main strand (default=0.75).
 -clsplit [0..1]         Minimum fraction of (normalized) read counts on the
                         smaller arm of a bi-directional piRNA cluster.
                         Otherwise the cluster will be annotated as
                         mono-directional (default=0.1).
 -noimage                Do not output image files for each piRNA cluster.
 -notable                Do not output a summary table.
 -nofaspi                Do not output a FASTA file comprising mapped piRNAs
                         for each cluster.
 -nofascl                Do not output a FASTA file comprising all piRNA
                         cluster sequences.
 -nomotif                Do not search for transcripotion factor binding
                         sites.
 -flank [int]            Include n bp of flanking sequence in output files.
 -pti                    Output a file that contains information on mapped
                         sequence reads in a tab-delimited table that
                         comprises sequence, reads, genomic hits e.g:
                         TGGGCACGCAAATTCGAGTATCG   12   4


';


GetOptions
	(
	"help"=>\$print_help,# print help/info text
	"h"=>\$print_help,# print help/info text
	"genome=s"=>\$genome,# genome in fasta format
	"map=s"=>\$map,# eland output from sRNAmapper or SeqMap using the option /output_all_matches
	"repeatmasker=s"=>\$RMannotation,# output from RepeatMasker
	"geneset=s"=>\$GeneSet,# Gene Set (GTF-file) from ensembl database
	"swsize=i"=>\$sliding_window_size,
	"swincr=i"=>\$sliding_window_increament,
	"nh"=>\$normalize_by_hit_number,# normalize hits by number of genomic hits
	"nr"=>\$normalize_by_read_number,# normalize hits by number of reads
	"rpm"=>\$normalize_by_total_number_of_mapped_reads,# normalize hits by number of reads
	"dens=f"=>\$threshold_density_absolute,# set user-defined minimum number of normalized hits per kb
	"pdens=f"=>\$threshold_density_pvalue,# p value for hit density
	"est=i"=>\$accuracy,# number of resampling steps when calculating minimum density based on expected hit density
	"pisize=f"=>\$threshold_size,# fraction of normalized reads with typical piRNA size
	"pimin=i"=>\$min_piRNAsize,# minimum size for a typical piRNA [nt]
	"pimax=i"=>\$max_piRNAsize,# maximum size for a typical piRNA [nt]
	"1Tor10A=f"=>\$threshold_T1A10,# fraction of normalized reads with 1T or 10A (either/or)
	"1Tand10A=f"=>\$threshold_T1A10_both,# fraction of normalized reads with 1T or 10A (both)
	"distr=s"=>\$avoid_peaks,# avoid false positive piRNA cluster prediction due to peaks of few mapped sequences with many reads (x-y which means: x% of top sequences account for max. y% of reads. Default: 5-50);
	"clsize=i"=>\$threshold_clustersize,# minimum size of a piRNA cluster [bp]
	"clhits=i"=>\$threshold_hits,# minimum number of sequence hits per cluster
	"clhitsn=f"=>\$threshold_hits_normalized,# minimum number of mapped reads (normalized) per cluster
	"clstrand=f"=>\$threshold_strandbias,# fraction of normalized reads on mainstrand
	"clsplit=f"=>\$bidirectionality_split,# minimum fraction of normalized hits per cluster arm (for bi-directional clusters)
	"noimage"=>\$image_files,# make an image file for each cluster [0=no/1=yes]
	"notable"=>\$results_table,# make a summary table with data for each cluster
	"nofaspi"=>\$fasta_piRNA_files,# make a fasta file for each cluster comprising mapped sequences [0=no/1=yes]
	"nofascl"=>\$fasta_cluster_file,# make one fasta file for comprising all cluster sequences [0=no/1=yes]
	"nomotif"=>\$search_bindingsites,# search sequence motifs in piRNA clusters (e.g. binding sites for transcription factors like A-MYB)
	"flank=i"=>\$flank_size,# flanking sequence [bp]
	"pti"=>\$save_info,# save information from map file and genome file in ~.pTi for faster future computation
	);

###   CHECK COMMAND LINE OPTIONS   ###
if($print_help){print$help_desk;exit;}
unless(-e$map){print"Could not open map file: $map.\n$!\n";exit;}
unless(-e$genome){print"Could not open genome file: $genome.\n$!\n";exit;}
if($RMannotation){unless(-e$RMannotation){print"Could not open RepeatMasker file: $RMannotation.\n$!\nIgnore and continue...\n";}}else{$RMannotation="n.a."}
if($GeneSet){unless(-e$GeneSet){print"Could not open GeneSet file: $GeneSet.\n$!\nIgnore and continue...\n";}}else{$GeneSet="n.a."}
unless($sliding_window_size){$sliding_window_size=5000;}
unless($sliding_window_increament){$sliding_window_increament=1000;}
unless($normalize_by_read_number){$normalize_by_read_number=0;$normalization_in_words="normalization: 1";}else{$normalization_in_words="normalization: 1*(reads)";}
unless($normalize_by_hit_number){$normalize_by_hit_number=0;}else{$normalization_in_words.="/(genomic hits)";}
unless($normalize_by_total_number_of_mapped_reads){$normalize_by_total_number_of_mapped_reads=0;}else{$normalization_in_words.="/((total mapped reads)*10^6)";}$normalization_in_words=~s/1$/no/;
if($accuracy){$exp_or_obs=1;unless($threshold_density_pvalue){$threshold_density_pvalue=0.05;}}
if($threshold_density_absolute&&$threshold_density_pvalue){print"ERROR: Do not use options -dens and -pdens together!\n";exit;}
elsif($threshold_density_absolute){$exp_or_obs=0;$options_info="Minimum hit density: $threshold_density_absolute hits/kb";}
elsif($threshold_density_pvalue&&$accuracy){$exp_or_obs=1;$threshold_density_absolute=0;$options_info="Significant (p<=$threshold_density_pvalue) hit density will be estimated based\non $accuracy random 1kb sliding windows.";}
#else{$exp_or_obs=2;$threshold_density_pvalue=0.01;$threshold_density_absolute=0;$options_info="Significant (p<=$threshold_density_pvalue) hit density will be calculated based\non observed hit distribution.";} #kw_mod
else{$exp_or_obs=2;unless($threshold_density_pvalue){$threshold_density_pvalue=0.01;}$threshold_density_absolute=0;$options_info="Significant (p<=$threshold_density_pvalue) hit density will be calculated based\non observed hit distribution.";} #kw_mod
unless($threshold_size){$threshold_size=0.75}
unless($min_piRNAsize){$min_piRNAsize=24;}
unless($max_piRNAsize){$max_piRNAsize=30;}
unless($threshold_T1A10){$threshold_T1A10=0.75;}
unless($threshold_T1A10_both){$threshold_T1A10_both=0.5;}
unless($avoid_peaks)
	{
	$avoid_peaks[0]=5;
	$avoid_peaks[1]=50;
	}
else
	{
	if($avoid_peaks=~/^\d+-\d+$/)
		{
		@avoid_peaks=split('-',$avoid_peaks);
		}
	else
		{
		print"\nERROR: Option -distr requires format: [int]-[int] (e.g. 5-50)";
		exit;
		}
	}
unless($threshold_clustersize){$threshold_clustersize=5000;}
unless($threshold_hits){$threshold_hits=0;}
unless($threshold_hits_normalized){$threshold_hits_normalized=0;}
unless($threshold_strandbias){$threshold_strandbias=0.75;}
unless($bidirectionality_split){$bidirectionality_split=0.1;}
unless($image_files){$image_files=1;}else{$image_files=0;}
unless($results_table){$results_table=1;}else{$results_table=0;}
unless($fasta_piRNA_files){$fasta_piRNA_files=1;}else{$fasta_piRNA_files=0;}
unless($fasta_cluster_file){$fasta_cluster_file=1;}else{$fasta_cluster_file=0;}
unless($search_bindingsites){$search_bindingsites=1;}else{$search_bindingsites=0;}
unless($flank_size){$flank_size=0;}
unless($save_info){$save_info=0;}

$proTRAC_runinfo=
"================================= proTRAC ====================================
VERSION: 2.0.5                                    LAST MODIFIED: 15. April 2015

Please cite:
Rosenkranz D, Zischler H. proTRAC - a software for probabilistic piRNA cluster
detection, visualization and analysis. 2012. BMC Bioinformatics 13:5.

Contact:
David Rosenkranz
Institute of Anthropology, small RNA goup
Johannes Gutenberg University Mainz
email: rosenkrd\@uni-mainz.de

You can find the latest proTRAC version at:
http://sourceforge.net/projects/protrac/files
http://www.smallRNAgroup-mainz.de/software
==============================================================================

PARAMETERS:
Map file: ...............$map
Genome file: ............$genome
RepeatMasker annotation: $RMannotation
GeneSet:.................$GeneSet

$options_info

Sliding window size: ........................................ $sliding_window_size bp
Sliding window increament: .................................. $sliding_window_increament bp
Normalize each hit by number of genomic hits: ............... $normalize_by_hit_number [0=no/1=yes]
Normalize each hit by number of sequence reads: ............. $normalize_by_read_number [0=no/1=yes]
Normalize values (-> per million mapped reads): ............. $normalize_by_total_number_of_mapped_reads [0=no/1=yes]
Min. fraction of hits with 1T(U) or 10A: .................... $threshold_T1A10
Alternatively: Min. fraction of hits with 1T(U) and 10A: .... $threshold_T1A10_both
Min. fraction of hits with typical piRNA length: ............ $threshold_size
Typical piRNA length: ....................................... $min_piRNAsize-$max_piRNAsize nt
Min. size of a piRNA cluster: ............................... $threshold_clustersize bp.
Min. number of hits (absolute): ............................. $threshold_hits
Min. number of hits (normalized): ........................... $threshold_hits_normalized
Min. fraction of hits on the mainstrand: .................... $threshold_strandbias
Top fraction of mapped sequences (in terms of read counts): . $avoid_peaks[0]%
Top fraction accounts for max. n% of sequence reads: ........ $avoid_peaks[1]%
Min. fraction of hits on each arm of a bidirectional cluster: $bidirectionality_split
Output image file for each cluster: ......................... $image_files [0=no/1=yes]
Output a summary table: ..................................... $results_table [0=no/1=yes]
Output a FASTA file for each cluster (piRNA sequences): ..... $fasta_piRNA_files [0=no/1=yes]
Output a FASTA file comprising cluster sequences: ........... $fasta_cluster_file [0=no/1=yes]
Search DNA motifs in clusters: .............................. $search_bindingsites [0=no/1=yes]
Output flanking sequences: +/- .............................. $flank_size bp
Output ~.pTi file: .......................................... $save_info [0=no/1=yes]
==============================================================================\n\n\n";
print"\n\n$proTRAC_runinfo";

# create folder for results
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$mon++;
$year+=1900;
$out_folder="proTRAC_$map\_$year\y$mon\m$mday\d$hour\h$min\m$sec\s";
if(-d $out_folder)
	{
	$folder_id=1;
	while(1)
		{
		$folder_id++;
		$out_folder="proTRAC_$year\y$mon\m$mday\d$hour\h$min\m$sec\s_$folder_id";
		unless(-d $out_folder)
			{
			mkdir("$out_folder")||die print"Could not create results folder.\n$!\n";
			last;
			}
		}
	}
else
	{
	mkdir("$out_folder")||die print"Could not create results folder.\n$!\n";
	}


# read genome
print"read genome ($genome)...";
%list_titles_genome=();
%scaffolds_to_analyze=();
open(GENOME,$genome)||die print"Could not open genome file.\n$!\n";
$genome_size=0;
while(<GENOME>)
	{
	$_=~s/\s*$//;
	if($_!~/^>/)
		{
		$genome{$title}.=$_;
		$_=~s/[NnXx-]//g;
		$gap_size+=length$&;
		$genome_size+=length$_;
		}
	else
		{
		# make a list of scaffolds/locations that are large enough to comprise a piRNA cluster
		$length_check=length$genome{$title};
		$_=~s/>//;
		$list_titles_genome{$_}=1;
		$title=$_;
		}
	}
close GENOME;
print" done.\ngenome size (without gaps): $genome_size bp\ngaps (N/X/-): $gap_size bp\n\n";

# read map file or info file
$read_info=$map.".pTi";
if(-e$read_info)
	{
	print"read information from $read_info...";
	open(INFO,$read_info);
	while(<INFO>)
		{
		$_=~s/\s*$//;
		@d=split("\t",$_);
		$total_reads+=$d[1];
		$genomic_hits+=$d[2];
		
		if($normalize_by_hit_number==0)
			{
			$hits_per_seq{$d[0]}=1;
			}
		else
			{
			$hits_per_seq{$d[0]}=$d[2];
			}
		
		if($normalize_by_read_number==0)
			{
			$reads_per_seq{$d[0]}=1;
			}
		else
			{
			$reads_per_seq{$d[0]}=$d[1];
			}
		}
	close INFO;
	}
else
	{
	print"read map file ($map)...";
	open(MAP,$map)||die print$!;
	while(<MAP>)
		{
		if($_!~/^trans_id/&&$_!~/^#/)
			{
			@d=split("\t",$_);
			$genomic_hits++;
			
			#unless($total_reads{$d[4]}) #kw_mod
			#	{ #kw_mod
				# $total_reads{$d[4]}=1; #kw_mod
				# $total_reads+=$d[3]; #kw_mod
			#	} #kw_mod
			$total_reads{$d[4]}+=$d[3]; #kw_mod
			$total_reads+=$d[3]; #kw_mod
			
			if($normalize_by_hit_number==1)
				{
				$hits_per_seq{$d[4]}++;
				}
			else
				{
				$hits_per_seq{$d[4]}=1;
				}
			
			if($normalize_by_read_number==1)
				{
				# $reads_per_seq{$d[4]}=$d[3]; #kw_mod
				$reads_per_seq{$d[4]}+=$d[3]; #kw_mod
				}
			else
				{
				$reads_per_seq{$d[4]}=1;
				}
			}
		}
	close MAP;
	}
$non_ID=keys%hits_per_seq;
print" done.\nmapped reads: $total_reads\nnon-identical sequences: $non_ID\nhit loci: $genomic_hits\n\n";

# save read counts and hit counts to .pTi-file
if($save_info==1)
	{
	unless(-e$read_info)
		{
		$save_info_out=$map.".pTi";
		open(SAVEINFO,">$save_info_out")||die print$!;
		foreach(keys%hits_per_seq)
			{
			print SAVEINFO"$_\t$reads_per_seq{$_}\t$hits_per_seq{$_}\n";
			}
		close SAVEINFO;
		}
	else
		{
		print"INFORMATION: $read_info already exists. Will not overwrite it.\n";
		}
	}

# read RepeatMasker output file
$RM_index=0;
%list_titles_RM=();
if(-e$RMannotation&&$image_files==1)
	{
	$RM_index=1;
	print"read RepeatMasker annotation $RMannotation...";
	open(RM,$RMannotation);
	while(<RM>)
		{
		if($_=~s/^ *\d+ +//)
			{
			@d=split(/ +/,$_);
			$s=int($d[4]/1000000);
			$loc="$d[3]#$s";
			$list_titles_RM{$d[3]}=1;
			$repeats_on_loc{$loc}++;
			$div_to_cons=$d[0]+$d[1]+$d[2];
			$RM{$loc}{$repeats_on_loc{$loc}}=
				[
				$d[4],
				$d[5],
				$d[8],
				$d[7],
				$div_to_cons
				];
			}
		}
	close RM;
	print" done.\n";
	# check if RM titles are consistent with genome titles
	foreach$rmtitle(keys%list_titles_RM)
		{
		unless($list_titles_genome{$rmtitle})
			{
			print"Location $rmtitle in RepeatMasker file is not part of $genome!\n";
			}
		}
	undef%list_titles_RM;
	}


# read Gene Set (GTF) file
$GTF_index=0;
if(-e$GeneSet&&$image_files==1)
	{
	$GTF_index=1;
	print"read Gene Set $GeneSet...";
	open(GTF,$GeneSet);
	while(<GTF>)
		{
		if($_!~/^#/)
			{
			@d=split("\t",$_);
			if($d[2]=~/exon/||$d[2]=~/UTR/)
				{
				$s=int($d[3]/1000000);
				$loc="$d[0]#$s";
				$list_titles_GTF{$d[0]}=1;
				$genes_on_loc{$loc}++;
				$gtf_info_text="";
				
				# extract information from gtf annotation
				$d[8]=~s/gene_id "[^"]+"//;
				$gene_id=$&;
				$gene_id=~s/^[^"]+"//;
				$gene_id=~s/"//;
				
				$d[8]=~s/gene_name "[^"]+"//;
				$gene_name=$&;
				$gene_name=~s/^[^"]+"//;
				$gene_name=~s/"//;
				
				$d[8]=~s/transcript_id "[^"]+"//;
				$transcript_id=$&;
				$transcript_id=~s/\d+//;
				$transcript_id=$&;
				
				$gtf_info_text="$gene_name ($gene_id) Tr:$transcript_id";
				if($d[8]=~s/exon_number "[^"]+"//)
					{
					$exon_number=$&;
					$exon_number=~s/\d+//;
					$gtf_info_text.=" Ex:$&";
					}
				elsif($d[2]=~/UTR/)
					{
					$gtf_info_text.=" UTR";
					}
				
				$GeneSet{$loc}{$genes_on_loc{$loc}}=
					[
					$d[3],
					$d[4],
					$gtf_info_text,
					$d[6]
					];
				}
			}
		}
	close GTF;
	print" done.\n";
	# check if GTF titles are consistent with genome titles
	foreach$gtftitle(keys%list_titles_GTF)
		{
		unless($list_titles_genome{$gtftitle})
			{
			print"Location $gtftitle in GeneSet file is not part of $genome!\n";
			}
		}
	undef%list_titles_RM;
	}
undef%list_titles_genome;


# load sequence motifs
if($search_bindingsites==1)
	{
	#$binding_motifs{'TBP'}='TATAAA';
	#$binding_motifs{'TBP rc'}='TTTATA';
	#$binding_motifs{'RFX-palindrome'}='CCTAGG';				# RFX palindrome from Horvath et al. 2004
	$binding_motifs{'RFX4_1'}='C.T[AG]GCAAC';			# top k-mer in UniProbe database
	$binding_motifs{'RFX4_1 rc'}='GTTGC[TC]A.G';		# top k-mer in UniProbe database
	$binding_motifs{'RFX4_2'}='C.T[AG]G[TA]TAC';		# top k-mer in UniProbe database
	$binding_motifs{'RFX4_2 rc'}='GTA[AT]C[TC]A.G';	# top k-mer in UniProbe database
	$binding_motifs{'Mybl1_1'}='AACCGTTA';			# top k-mer in UniProbe database
	$binding_motifs{'Mybl1_1 rc'}='TAACGGTT';			# top k-mer in UniProbe database
	$binding_motifs{'A-MYB'}='[TA]G[GA]CAGTTGG';		# motif from Li et al. Mol Cell. 2013 50(1): 67–81
	$binding_motifs{'A-MYB rc'}='CCAACTG[CT]C[AT]';	# motif from Li et al. Mol Cell. 2013 50(1): 67–81
	}


if($normalize_by_read_number==0)
	{
	$total_reads=$non_ID;
	}
if($exp_or_obs==1)
	{
	# estimate significant hit density (normalized) per kb...
	print"estimating significant (p<=$threshold_density_pvalue) hit density...";
	$exp_hitdensity=$non_ID/$genome_size*1000*($total_reads/$non_ID);
	@values=values%hits_per_seq;
	foreach(1..$accuracy)
		{
		$i=0;
		foreach(1..1000)
			{
			if((rand)<=$genomic_hits/$genome_size)
				{
				$i+=$non_ID/$genomic_hits*($total_reads/$non_ID);
				}
			}
		push(@i,$i);
		}
	@i=sort{$b<=>$a}@i;
	foreach(1..$accuracy*$threshold_density_pvalue)
		{
		$last_value=shift@i;
		if($last_value>0)
			{
			$rescue_value=$last_value;
			}
		}
	$sig_hitdensity=shift@i;
	if($sig_hitdensity==0)
		{
		$sig_hitdensity=$rescue_value;
		}
	
	undef@i;
	print" done.\nexpectation: $exp_hitdensity hits/kb*\nsignificant density: $sig_hitdensity hits/kb*\n*=normalized\n\n";
	}

elsif($exp_or_obs==2)
	{
	# or check density distribution
	print"\ncheck genomic distribution of hit density:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
	open(MAP,$map)||die print"Could not open map file: $map.\n$!\n";
	if($^O=~/MSWin/)
		{
		$seek=1; # -1 is necessary for windows systems because newline is represented by 2 characters.
		}
	else
		{
		$seek=0; # on linux or mac systems a newline is represented by 1 character.
		}
	$processed_genomic_hits=0;
	$printed_dots=0;
	@hit_densities=();
	@sw=('0')x$sliding_window_size;
	while(<MAP>)
		{
		$processed_genomic_hits++;
		if($processed_genomic_hits>=$genomic_hits/71)
			{
			$processed_genomic_hits=0;
			print".";
			$printed_dots++;
			}
		if($_!~/^trans_id/&&$_!~/^#/)
			{
			@d=split("\t",$_);
			
			# change of chromosome/scaffold/reference
			if($d[0]ne$prev_scaff)
				{
				# process previous sliding window 
				$sw_value=0;
				foreach$pos(0..$sliding_window_size-1)
					{
					$sw_value+=$sw[$pos];
					}
				push(@hit_densities,$sw_value);
				$sw_start=1;
				@sw=('0')x$sliding_window_size;
				
				if($^O=~/MSWin/)
					{
					seek(MAP,-(length$_)-1,1); # -1 is necessary for windows systems because newline is represented by 2 characters.
					}
				else
					{
					seek(MAP,-length$_,1); # on linux or mac systems a newline is represented by 1 character.
					}
				
				$processed_genomic_hits=$processed_genomic_hits-1;
				}
			
			# hit is within the sliding window
			elsif($d[1]>=$sw_start&&$d[1]<$sw_start+$sliding_window_size)
				{
				# use standard normalized hit counts
				if(@d==7)
					{
					$sw[$d[1]-$sw_start]+=(1/$hits_per_seq{$d[4]})*$reads_per_seq{$d[4]};
					}
				# use weighted hit counts as calculated by map-arrange.pl
				elsif(@d==9)
					{
					$sw[$d[1]-$sw_start]+=$d[8];
					}
				}
			
			# process previous sliding window and move sliding window
			else
				{
				$sw_value=0;
				foreach$pos(0..$sliding_window_size-1)
					{
					$sw_value+=$sw[$pos];
					}
				push(@hit_densities,$sw_value);
				splice(@sw,0,$sliding_window_increament);
				$sw_start+=$sliding_window_increament;
				seek(MAP,-(length$_)-$seek,1); # -1 is necessary for windows systems because newline is represented by 2 characters.
				$processed_genomic_hits=$processed_genomic_hits-1;
				}
			$prev_scaff=$d[0];
			}
		}
	close MAP;
	@hit_densities=sort{$b<=>$a}@hit_densities;
	$sig_hitdensity=$hit_densities[abs($threshold_density_pvalue*@hit_densities)];
	
	# check if significant hit density = 0
	while($sig_hitdensity==0)
		{
		$sig_hitdensity=pop@hit_densities;
		last if(@hit_densities==0);
		}
	if($sig_hitdensity==0)
		{
		die print"\nERROR: Hit density is 0 for all sliding windows!\nCheck map file.\n\n";
		}
	foreach($printed_dots..70)
		{
		print".";
		}
	print"\nsignificant density: $sig_hitdensity hits/kb*\n*=normalized\n\n";
	undef@hit_densities;
	}
else
	{
	$sig_hitdensity=$threshold_density_absolute;
	}

# read map file with sliding window
$prev_trans_id="";
$sw_id=1-($sliding_window_increament/$sliding_window_size);
@sw=();
print"search and tag loci with significant hit density:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
open(MAP,$map)||die print$!;
$processed_genomic_hits=0;
$printed_dots=0;
while(<MAP>)
	{
	$processed_genomic_hits++;
	if($processed_genomic_hits>=$genomic_hits/71)
		{
		$processed_genomic_hits=0;
		print".";
		$printed_dots++;
		}
	
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		@d=split("\t",$_);
		
		# check if hit is on the same chromosome/scaffold
		if($d[0]ne$prev_trans_id)
			{
			$sw_id=1-($sliding_window_increament/$sliding_window_size);
			if(@sw>0)
				{
				calc_and_check();
				}
			@sw=();
			}
		
		# add hit to array if coordinate is in the sliding window
		if($d[1]<=$sw_id*$sliding_window_size)
			{
			push(@sw,$_);
			}
		
		elsif(@sw>0)
			{
			calc_and_check();
			sub calc_and_check
				{
				# calculate hit density in sliding window
				$sw_hits=0;
				%sw_coord=();
				$count=0;
				foreach$hit(@sw)
					{
					@hit=split("\t",$hit);
					$sw_coord{"$hit[0]-$hit[1]"}=1;
					
					# use standard normalized hit counts
					if(@hit==7)
						{
						$sw_hits+=(1/$hits_per_seq{$hit[4]})*$reads_per_seq{$hit[4]};
						}
					# use weighted hit counts as calculated by map-arrange.pl
					elsif(@hit==9)
						{
						$sw_hits+=$hit[8];
						}
					}
				
				# check hit density
				if(($sw_hits/$sliding_window_size*1000)>=$sig_hitdensity&&($sw_hits/$sliding_window_size*1000)>=$threshold_density_absolute)
					{
					foreach$key(keys%sw_coord)
						{
						$all_coord{$key}=1;
						}
					}
				undef%sw_coord;
				}
			
			# move sliding window
			move_sw();
			sub move_sw
				{
				$sw_increment_steps=0;
				while($d[1]>$sw_id*$sliding_window_size)
					{
					$sw_id+=($sliding_window_increament/$sliding_window_size);
					$sw_increment_steps++;
					}
				if($sw_increment_steps==1)
					{
					@new_sw=();
					foreach$hit(@sw)
						{
						@hit=split("\t",$hit);
						if($hit[1]>($sw_id-1)*$sliding_window_size)
							{
							push(@new_sw,$hit);
							}
						}
					@sw=@new_sw;
					undef@new_sw;
					}
				else
					{
					@sw=();
					}
				push(@sw,$_);
				}
			}
		else
			{
			move_sw();
			}
		$prev_trans_id=$d[0];
		}
	}
close MAP;
foreach($printed_dots..70)
	{
	print".";
	}
print"\n\n";

# print to output files
if($image_files==1)
	{
	eval"use GD; 1" or $module_unavailable=1;
	if($module_unavailable==1)
		{
		print"\nWARNING: GD module not found on this machine.\nYou need to install GD in order to output piRNA cluster image files.\n\n";
		$image_files=0;
		}
	else
		{
		# load module again because loading module in eval will cause problems with GD::Image::string and GD::Font
		use GD;
		}
	}
$id=0;
$stat=0;
open(MAP,$map)||die print$!;

if($fasta_cluster_file==1)
	{
	open(CLUSTER_FASTA,">$out_folder/clusters.fasta");
	}
if($results_table==1)
	{
	if($exp_or_obs==0)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb (user-defined): $threshold_density_absolute";
		}
	elsif($exp_or_obs==1)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb based on theoretically expected hit density. p=$threshold_density_pvalue -> $sig_hitdensity hits/kb\n1kb-window resampling steps to calculate p: $accuracy";
		}
	elsif($exp_or_obs==2)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb based on observed hit density. p=$threshold_density_pvalue -> $sig_hitdensity hits/kb";
		}
	open(RESULTS_TEXT,">$out_folder/results.table");
	print RESULTS_TEXT$proTRAC_runinfo;
	}

$total_size_of_all_clusters=0;
%total_clustered_sequences=();
@flank_buffer=();
$prev_cl_coordinate=0;
$prev_hit_location=0;
while(<MAP>)
	{
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$_=~s/\s*$//;
		@hit_data=split("\t",$_);
		
		# check if chromosome/scaffold changed and last chromosome/scaffold hit is inside a cluster candidate
		if($hit_data[0]ne$prev_hit_location&&$stat==1)
			{
			cluster_candidate_check();
			$stat=0;
			}
		$prev_hit_location=$hit_data[0];
		
		
		# store 5' hits if flank_size > 0
		if($flank_size>0)
			{
			push(@flank_buffer,$_);
			if(@flank_buffer>$flank_size*10)
				{
				shift@flank_buffer;
				}
			}
		
		
		if($all_coord{"$hit_data[0]-$hit_data[1]"})
			{
			# check if the gap between tagged hits is bigger than 10kb -> split
			if($hit_data[1]-$prev_cl_coordinate>=10000&&$stat==1)
				{
				cluster_candidate_check();
				$stat=0;
				}
			$prev_cl_coordinate=$hit_data[1];
			
			if($stat==0)
				{
				$hits_abs=0;
				$hits_norm=0;
				$cl_T1=0;
				$cl_A10=0;
				$cl_size=0;
				$id++;
				$start=$hit_data[1];
				@shape=();
				open(OUT,">$out_folder/$id.fasta");
				
				# check if there are 5' flanking hits if flank_size > 0
				if($flank_size>0)
					{
					pop@flank_buffer;
					foreach$flank_hit(@flank_buffer)
						{
						@d_flank=split("\t",$flank_hit);
						if($d_flank[0]eq$hit_data[0]&&$d_flank[1]>=$hit_data[1]-$flank_size)
							{
							unless($d_flank[8])
								{
								$d_flank[8]=$reads_per_seq{$d_flank[4]}/$hits_per_seq{$d_flank[4]};
								}
							print OUT ">$d_flank[0]\t$d_flank[1]\t$reads_per_seq{$d_flank[4]}\t$hits_per_seq{$d_flank[4]}\t$d_flank[8]\t$d_flank[6]\tFLANK\n$d_flank[4]\n";
							}
						}
					}
				}
			$hits_abs++;
			
			if(@hit_data==7)
				{
				$hits_norm+=(1/$hits_per_seq{$hit_data[4]})*$reads_per_seq{$hit_data[4]};
				}
			elsif(@hit_data==9)
				{
				$hits_norm+=$hit_data[8]
				}
			
			if($hit_data[4]=~/^T/)
				{
				if(@hit_data==7)
					{
					$cl_T1+=(1/$hits_per_seq{$hit_data[4]})*$reads_per_seq{$hit_data[4]};
					}
				elsif(@hit_data==9)
					{
					$cl_T1+=$hit_data[8]
					}
				}
			if($hit_data[4]=~/^.{9}A/)
				{
				if(@hit_data==7)
					{
					$cl_A10+=(1/$hits_per_seq{$hit_data[4]})*$reads_per_seq{$hit_data[4]};
					}
				elsif(@hit_data==9)
					{
					$cl_A10+=$hit_data[8]
					}
				}
			if(length$hit_data[4]>=$min_piRNAsize&&length$hit_data[4]<=$max_piRNAsize)
				{
				if(@hit_data==7)
					{
					$cl_size+=(1/$hits_per_seq{$hit_data[4]})*$reads_per_seq{$hit_data[4]};
					}
				elsif(@hit_data==9)
					{
					$cl_size+=$hit_data[8]
					}
				}
			$cl_location=$hit_data[0];
			$end=($hit_data[1]+length$hit_data[4])-1;
			$seqs_in_candidate{$hit_data[4]}=1;
			unless($hit_data[8])
				{
				$hit_data[8]=$reads_per_seq{$hit_data[4]}/$hits_per_seq{$hit_data[4]};
				}
			print OUT ">$hit_data[0]\t$hit_data[1]\t$reads_per_seq{$hit_data[4]}\t$hits_per_seq{$hit_data[4]}\t$hit_data[8]\t$hit_data[6]\n$hit_data[4]\n";
			if(@hit_data==7)
				{
				push(@shape,$reads_per_seq{$hit_data[4]}/$hits_per_seq{$hit_data[4]});
				}
			elsif(@hit_data==9)
				{
				push(@shape,$hit_data[8]);
				}
			
			
			$stat=1;
			}
		
		# check if there are 3' flanking hits if flank_size > 0
		elsif($hit_data[0]eq$cl_location&&$hit_data[1]-$flank_size<=$end)
			{
			unless($hit_data[8])
				{
				$hit_data[8]=$reads_per_seq{$hit_data[4]}/$hits_per_seq{$hit_data[4]};
				}
			print OUT ">$hit_data[0]\t$hit_data[1]\t$reads_per_seq{$hit_data[4]}\t$hits_per_seq{$hit_data[4]}\t$hit_data[8]\t$hit_data[6]\tFLANK\n$hit_data[4]\n";
			}
		
		else
			{
			if($stat==1)
				{
				cluster_candidate_check();
				sub cluster_candidate_check
					{
					close OUT;
					
					# check which value is lower, 1T or 10A | compare to $threshold_T1A10_both
					if($cl_A10<$cl_T1)
						{
						$low1T10A=$cl_A10/$hits_norm;
						}
					else
						{
						$low1T10A=$cl_T1/$hits_norm;
						}
					
					# calculate hit count for top fraction
					@shape=sort{$b<=>$a}@shape;
					$top_hits=0;
					foreach(0..int((@shape*($avoid_peaks[0]/100))-1))
						{
						$top_hits+=$shape[$_];
						}
					
					# check whether cluster results from one or few peaks
					if($top_hits>$hits_norm*($avoid_peaks[1]/100))
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# check cluster size
					elsif($end-$start+1<$threshold_clustersize)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# check number of hits (absolute)
					elsif($hits_abs<$threshold_hits)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# check number of hits (normalized)
					elsif($hits_norm<$threshold_hits_normalized)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# check number of hits (normalized) with typical piRNA size
					elsif($cl_size/$hits_norm<$threshold_size)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# check number of hits (normalized) with 1T or 10A
					elsif($cl_T1/$hits_norm<$threshold_T1A10&&$cl_A10/$hits_norm<$threshold_T1A10&&$low1T10A<$threshold_T1A10_both)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					# calculate strand bias
					else
						{
						open(CHECKSTRAND,"$out_folder/$id.fasta");
						@values=();
						@coordinates_for_split=();
						$plus_strand=0;
						while(<CHECKSTRAND>)
							{
							chomp$_;
							if($_=~/^>/&&$_!~/FLANK/)
								{
								@head=split("\t",$_);
								push(@coordinates_for_split,$head[1]);
								
								if($head[5]eq'+')
									{
									$plus_strand+=$head[4];
									push(@values,$head[4]);
									}
								else
									{
									push(@values,($head[4])*-1);
									}
								}
							}
						close CHECKSTRAND;
						
						# calculate monodirectionality
						if($plus_strand/$hits_norm>=0.5)
							{
							$directionality="mono:plus";
							$strandbias=$plus_strand/$hits_norm;
							}
						else
							{
							$directionality="mono:minus";
							$strandbias=($hits_norm-$plus_strand)/$hits_norm;
							}
						
						# calculate bidirectionality
						@left_arm=();
						@left_fraction=();
						$fraction=0;
						$plus_strand=0;
						foreach$value(@values)
							{
							$replicate=$value;
							$fraction+=abs$replicate;
							push(@left_fraction,$fraction);
							if($value>0)
								{
								$plus_strand+=$value;
								}
							push(@left_arm,$plus_strand);
							}
							
						@values=reverse@values;
						
						@right_arm=();
						@right_fraction=();
						$fraction=0;
						$minus_strand=0;
						foreach$value(@values)
							{
							$replicate=$value;
							$fraction+=abs$replicate;
							push(@right_fraction,$fraction);
							if($value<0)
								{
								$minus_strand=$minus_strand-$value;
								}
							push(@right_arm,$minus_strand);
							}
						
						$best_split=-1;
						foreach$split_pos(0..@values-2)
							{
							# both arms with a minimum of normalized hits?
							if($left_fraction[$split_pos]/$hits_norm>=$bidirectionality_split&&$right_fraction[-2-$split_pos]/$hits_norm>=$bidirectionality_split)
								{
								# left arm with proper strans bias?
								if($left_arm[$split_pos]/$left_fraction[$split_pos]>=$threshold_strandbias||$left_arm[$split_pos]/$left_fraction[$split_pos]<=1-$threshold_strandbias)
									{
									# right arm with proper strans bias?
									if($right_arm[-2-$split_pos]/$right_fraction[-2-$split_pos]>=$threshold_strandbias||$right_arm[-2-$split_pos]/$right_fraction[-2-$split_pos]<=1-$threshold_strandbias)
										{
										$split_value=$left_arm[$split_pos]+$right_arm[-2-$split_pos];
										if($split_value/$hits_norm>0.5)
											{
											$plus_minus="plus-minus";
											$split_value=$split_value/$hits_norm;
											}
										else
											{
											$plus_minus="minus-plus";
											$split_value=($hits_norm-$split_value)/$hits_norm;
											}
										if($split_value>$strandbias)
											{
											$best_split=$split_pos;
											$strandbias=$split_value;
											$directionality="bi:$plus_minus";
											}
										}
									}
								}
							}
						
						# check strand bias
						if($strandbias<$threshold_strandbias)
							{
							$directionality="dual"; #kw_mod
							#unlink"$out_folder/$id.fasta"; #kw_mod
							#$id=$id-1; #kw_mod
							}
						if(0) {} #kw_mod	
						else
							{
							print"Predicted Cluster: $cl_location $start-$end directionality: $directionality\n";
							if($fasta_cluster_file==1)
								{
								$cluster_sequence=substr($genome{$cl_location},$start-1-$flank_size,$end-$start+1+($flank_size*2));
								print CLUSTER_FASTA"$cl_location $start-$end (+-$flank_size bp) directionality: $directionality\n$cluster_sequence\n";
								}
							foreach$cl_seq(keys%seqs_in_candidate)
								{
								$total_clustered_sequences{$cl_seq}=1;
								}
							if($image_files==1||$results_table==1)
								{
								$picture_height=720;
								
								# check number of RepeatMasker elements / GTF elements and make the image larger (if required)
								if($RM_index==1||$GTF_index==1)
									{
									$count_RMelements=0;
									$count_GTFelements=0;
									$s=int($start/1000000);
									$loc="$cl_location#$s";
									countRMelements();
									sub countRMelements
										{
										$repeat_id=0;
										$repeat_in_cl=0;
										while(1)
											{
											$repeat_id++;
											if($RM{$loc}{$repeat_id})
												{
												$ref=$RM{$loc}{$repeat_id};
												@rm=@$ref;
												# repeat inside piRNA cluster
												if($rm[0]>=$start&&$rm[1]<=$end)
													{
													$count_RMelements++;
													}
												# repeat overlaps the 5' end of the piRNA cluster
												elsif($rm[0]<=$start&&$rm[1]>=$start)
													{
													$count_RMelements++;
													}
												# repeat overlaps the 3' end of the piRNA cluster
												elsif($rm[0]<=$end&&$rm[1]>=$end)
													{
													$count_RMelements++;
													}
												}
											else
												{
												last;
												}
											}
										}
									countGTFelements();
									sub countGTFelements
										{
										$gene_id=0;
										$gene_in_cl=0;
										while(1)
											{
											$gene_id++;
											if($GeneSet{$loc}{$gene_id})
												{
												$ref=$GeneSet{$loc}{$gene_id};
												@gtf=@$ref;
												# gene inside piRNA cluster
												if($gtf[0]>=$start&&$gtf[1]<=$end)
													{
													$count_GTFelements++;
													}
												# gene overlaps the 5' end of the piRNA cluster
												elsif($gtf[0]<=$start&&$gtf[1]>=$start)
													{
													$count_GTFelements++;
													}
												# gene overlaps the 3' end of the piRNA cluster
												elsif($gtf[0]<=$end&&$gtf[1]>=$end)
													{
													$count_GTFelements++;
													}
												}
											else
												{
												last;
												}
											}
										}
									
									if(int($start/1000000)<int($end/1000000))
										{
										$s=int($end/1000000);
										$loc="$cl_location#$s";
										countRMelements();
										countGTFelements();
										}
									if($count_RMelements>17||$count_GTFelements>17)
										{
										if($count_RMelements>$count_GTFelements)
											{
											$enlarge_image=$count_RMelements;
											}
										else
											{
											$enlarge_image=$count_GTFelements;
											}
										$enlarge_image=$enlarge_image-17;
										$picture_height+=$enlarge_image*10;
										}
									}
								
								$image=new GD::Image 1100,$picture_height;
								
								$white=$image->colorAllocate(255,255,255);
								$black=$image->colorAllocate(0,0,0);
								$red=$image->colorAllocate(240,90,90);
								$light_red=$image->colorAllocate(255,200,200);
								$blue=$image->colorAllocate(90,150,240);
								$light_blue=$image->colorAllocate(200,220,255);
								$gray=$image->colorAllocate(200,200,200);
								$dark_gray=$image->colorAllocate(150,150,150);
								$light_gray=$image->colorAllocate(240,240,240);
								$color[0]=$image->colorAllocate(0,150,0);
								$color[1]=$image->colorAllocate(0,200,0);
								$color[2]=$image->colorAllocate(100,200,0);
								$color[3]=$image->colorAllocate(200,200,0);
								$color[4]=$image->colorAllocate(200,140,0);
								$color[5]=$image->colorAllocate(200,70,0);
								$color[6]=$image->colorAllocate(200,0,0);
								
								$image->string(gdMediumBoldFont,10,5,"Cluster $id",$black);
								$image->string(gdSmallFont,100,35,"hits ($normalization_in_words)",$black);
								$image->string(gdSmallFont,100,285,"hits (color code refers to genomic hits of the read)",$black);
								$image->string(gdTinyFont,50,140,"+ strand",$black);
								$image->string(gdTinyFont,50,152,"- strand",$black);
								$image->string(gdTinyFont,50,390,"+ strand",$black);
								$image->string(gdTinyFont,50,402,"- strand",$black);
								
								# draw RepeatMasker annotation
								if($RM_index==1)
									{
									$image->string(gdTinyFont,30,259,"RepeatMasker:",$black);
									$image->line(99,263,600,263,$black);
									
									# colors encode identity of the repeat to its consensus
									$rm_blue[0]=$image->colorAllocate(0,0,210); # 100-98% identity
									$rm_blue[1]=$image->colorAllocate(0,0,230); # <98-95% identity
									$rm_blue[2]=$image->colorAllocate(0,0,250); # <95-90% identity
									$rm_blue[3]=$image->colorAllocate(60,60,250); # <90-85% identity
									$rm_blue[4]=$image->colorAllocate(100,100,250); # <85-80% identity
									$rm_blue[5]=$image->colorAllocate(140,140,250); # <80-75% identity
									$rm_blue[6]=$image->colorAllocate(180,180,250); # <75-70% identity
									$rm_blue[7]=$image->colorAllocate(210,210,250); # <70% identity
									$rm_red[0]=$image->colorAllocate(210,0,0); # 100-98% identity
									$rm_red[1]=$image->colorAllocate(230,0,0); # <98-95% identity
									$rm_red[2]=$image->colorAllocate(250,0,0); # <95-90% identity
									$rm_red[3]=$image->colorAllocate(250,60,60); # <90-85% identity
									$rm_red[4]=$image->colorAllocate(250,100,100); # <85-80% identity
									$rm_red[5]=$image->colorAllocate(250,140,140); # <80-75% identity
									$rm_red[6]=$image->colorAllocate(250,180,180); # <75-70% identity
									$rm_red[7]=$image->colorAllocate(250,210,210); # <70% identity
									$image->string(gdSmallFont,410,530,"RepeatMasker annotation:",$black);
									$image->string(gdTinyFont,411,548,"+",$white);
									$image->string(gdTinyFont,426,548,"-",$white);
									$image->rectangle(410,548,420,558,$rm_blue[0]);$image->fill(411,549,$rm_blue[0]);$image->rectangle(425,548,435,558,$rm_red[0]);$image->fill(426,549,$rm_red[0]);$image->string(gdTinyFont,440,549,"100-98% identity",$black);
									$image->rectangle(410,558,420,568,$rm_blue[1]);$image->fill(411,559,$rm_blue[1]);$image->rectangle(425,558,435,568,$rm_red[1]);$image->fill(426,559,$rm_red[1]);$image->string(gdTinyFont,440,559,"<98-95% identity",$black);
									$image->rectangle(410,568,420,578,$rm_blue[2]);$image->fill(411,569,$rm_blue[2]);$image->rectangle(425,568,435,578,$rm_red[2]);$image->fill(426,569,$rm_red[2]);$image->string(gdTinyFont,440,569,"<95-90% identity",$black);
									$image->rectangle(410,578,420,588,$rm_blue[3]);$image->fill(411,579,$rm_blue[3]);$image->rectangle(425,578,435,588,$rm_red[3]);$image->fill(426,579,$rm_red[3]);$image->string(gdTinyFont,440,579,"<90-85% identity",$black);
									$image->rectangle(410,588,420,598,$rm_blue[4]);$image->fill(411,589,$rm_blue[4]);$image->rectangle(425,588,435,598,$rm_red[4]);$image->fill(426,589,$rm_red[4]);$image->string(gdTinyFont,440,589,"<85-80% identity",$black);
									$image->rectangle(410,598,420,608,$rm_blue[5]);$image->fill(411,599,$rm_blue[5]);$image->rectangle(425,598,435,608,$rm_red[5]);$image->fill(426,599,$rm_red[5]);$image->string(gdTinyFont,440,599,"<80-75% identity",$black);
									$image->rectangle(410,608,420,618,$rm_blue[6]);$image->fill(411,609,$rm_blue[6]);$image->rectangle(425,608,435,618,$rm_red[6]);$image->fill(426,609,$rm_red[6]);$image->string(gdTinyFont,440,609,"<75-70% identity",$black);
									$image->rectangle(410,618,420,628,$rm_blue[7]);$image->fill(411,619,$rm_blue[7]);$image->rectangle(425,618,435,628,$rm_red[7]);$image->fill(426,619,$rm_red[7]);$image->string(gdTinyFont,440,619,"<70% identity",$black);
									
									$s=int($start/1000000);
									$loc="$cl_location#$s";
									$repeat_id=0;
									$repeat_in_cl=0;
									draw_RM();
									sub draw_RM
										{
										while(1)
											{
											$repeat_id++;
											if($RM{$loc}{$repeat_id})
												{
												$ref=$RM{$loc}{$repeat_id};
												@rm=@$ref;
												$cluster_plus_flank=$end-$start+1+($flank_size*2);
												# repeat inside piRNA cluster
												if($rm[0]>=$start&&$rm[1]<=$end)
													{
													$rm_start=$rm[0];
													$rm_end=$rm[1];
													draw_RMelement();
													sub draw_RMelement
														{
														$repeat_in_cl++;
														if($rm[4]<2){$rm_color_code=0;}
														elsif($rm[4]<5){$rm_color_code=1;}
														elsif($rm[4]<10){$rm_color_code=2;}
														elsif($rm[4]<15){$rm_color_code=3;}
														elsif($rm[4]<20){$rm_color_code=4;}
														elsif($rm[4]<25){$rm_color_code=5;}
														elsif($rm[4]<30){$rm_color_code=6;}
														else{$rm_color_code=7;}
														
														foreach$pos($rm_start..$rm_end)
															{
															last if($pos>$end);
															$relative_position=$pos-$start+1+$flank_size;
															if($rm[3]=~/\+/)
																{
																$image->line(($relative_position/$cluster_plus_flank)*500+100,262,($relative_position/$cluster_plus_flank)*500+100,255,$rm_blue[$rm_color_code]);
																}
															elsif($rm[3]=~/C/)
																{
																$image->line(($relative_position/$cluster_plus_flank)*500+100,264,($relative_position/$cluster_plus_flank)*500+100,271,$rm_red[$rm_color_code]);
																}
															}
														if($rm[3]=~/\+/)
															{
															$image->string(gdTinyFont,(($rm_start-$start+1+$flank_size)/$cluster_plus_flank)*500+100,255,$repeat_in_cl,$black);
															}
														elsif($rm[3]=~/C/)
															{
															$image->string(gdTinyFont,(($rm_start-$start+1+$flank_size)/$cluster_plus_flank)*500+100,264,$repeat_in_cl,$black);
															}
														$image->string(gdTinyFont,550,(539+($repeat_in_cl*10)),"$repeat_in_cl: $rm[2] $rm[0]-$rm[1]",$black);
														}
													}
												# repeat overlaps the 5' end of the piRNA cluster
												elsif($rm[0]<=$start&&$rm[1]>=$start)
													{
													$rm_start=$start;
													$rm_end=$rm[1];
													draw_RMelement();
													}
												# repeat overlaps the 3' end of the piRNA cluster
												elsif($rm[0]<=$end&&$rm[1]>=$end)
													{
													$rm_start=$rm[0];
													$rm_end=$end;
													draw_RMelement();
													}
												}
											else
												{
												last;
												}
											}
										}
									if(int($start/1000000)<int($end/1000000))
										{
										$s=int($end/1000000);
										$loc="$cl_location#$s";
										draw_RM();
										}
									}
								
								# draw Gene Set (GTF) information
								if($GTF_index==1)
									{
									$image->string(gdTinyFont,50,273,"Gene Set:",$black);
									$image->line(99,277,600,277,$black);
									$s=int($start/1000000);
									$loc="$cl_location#$s";
									$gene_id=0;
									$gene_in_cl=0;
									draw_GTF();
									$image->string(gdSmallFont,760,530,"Gene Set:",$black);
									$image->rectangle(820,530,830,540,$blue);$image->fill(821,531,$blue);
									$image->rectangle(840,530,850,540,$red);$image->fill(841,531,$red);
									$image->string(gdTinyFont,822,532,"+",$black);
									$image->string(gdTinyFont,843,531,"-",$black);
									sub draw_GTF
										{
										while(1)
											{
											$gene_id++;
											if($GeneSet{$loc}{$gene_id})
												{
												$ref=$GeneSet{$loc}{$gene_id};
												@gtf=@$ref;
												$cluster_plus_flank=$end-$start+1+($flank_size*2);
												# gene inside piRNA cluster
												if($gtf[0]>=$start&&$gtf[1]<=$end)
													{
													$gtf_start=$gtf[0];
													$gtf_end=$gtf[1];
													draw_GTFelement();
													sub draw_GTFelement
														{
														$gene_in_cl++;
														if($gtf[3]=~/\+/)
															{
															$gtf_color=$blue;
															}
														elsif($gtf[3]=~/-/)
															{
															$gtf_color=$red;
															}
														foreach$pos($gtf_start..$gtf_end)
															{
															last if($pos>$end);
															$relative_position=$pos-$start+1+$flank_size;
															$image->line(($relative_position/$cluster_plus_flank)*500+100,281,($relative_position/$cluster_plus_flank)*500+100,273,$gtf_color);
															}
														$image->string(gdTinyFont,(($gtf_start-$start+1+$flank_size)/$cluster_plus_flank)*500+100,274,$gene_in_cl,$black);
														$image->string(gdTinyFont,760,(539+($gene_in_cl*10)),"$gene_in_cl: $gtf[2] $gtf[0]-$gtf[1]",$black);
														}
													}
												# gene overlaps the 5' end of the piRNA cluster
												elsif($gtf[0]<=$start&&$gtf[1]>=$start)
													{
													$gtf_start=$start;
													$gtf_end=$gtf[1];
													draw_GTFelement();
													}
												# gene overlaps the 3' end of the piRNA cluster
												elsif($gtf[0]<=$end&&$gtf[1]>=$end)
													{
													$gtf_start=$gtf[0];
													$gtf_end=$end;
													draw_GTFelement();
													}
												}
											else
												{
												last;
												}
											}
										}
									if(int($start/1000000)<int($end/1000000))
										{
										$s=int($end/1000000);
										$loc="$cl_location#$s";
										draw_GTF();
										}
									}
								
								$print_1T=(int((($cl_T1/$hits_norm)*1000)+0.5))/10;
								$print_10A=(int((($cl_A10/$hits_norm)*1000)+0.5))/10;
								$print_size=(int((($cl_size/$hits_norm)*1000)+0.5))/10;
								$print_strandbias=(int(($strandbias*1000)+0.5))/10;
								$print_clustersize=$end-$start+1;
								$print_density=($hits_norm/$print_clustersize)*1000;
								$total_size_of_all_clusters+=$print_clustersize;
								$cl_location_for_image=$cl_location;
								
								# calculate values per million mapped reads (rpm)
								if($normalize_by_total_number_of_mapped_reads==1)
									{
									$hits_norm=($hits_norm/$total_reads)*1000000;
									$print_density=($print_density/$total_reads)*1000000;
									}
								
								if(length$cl_location>39)
									{
									$cl_location_for_image=substr($cl_location,0,34);
									$cl_location_for_image.='[...]';
									}
								$image->string(gdSmallFont,100,530,"Location: $cl_location_for_image",$black);
								$image->string(gdSmallFont,100,545,"Coordinates: $start-$end",$black);
								$image->string(gdSmallFont,100,560,"Size [bp]: $print_clustersize",$black);
								$image->string(gdSmallFont,100,575,"Hits (absolute, non-identical): $hits_abs",$black);
								$image->string(gdSmallFont,100,590,"Hits (normalized): $hits_norm",$black);
								$image->string(gdSmallFont,100,605,"Hits (normalized) per kb: $print_density",$black);
								$image->string(gdSmallFont,100,620,"Normalized hits with 1T: $print_1T%",$black);
								$image->string(gdSmallFont,100,635,"Normalized hits with 10A: $print_10A%",$black);
								$image->string(gdSmallFont,100,650,"Normalized hits with length $min_piRNAsize-$max_piRNAsize nt: $print_size%",$black);
								$image->string(gdSmallFont,100,665,"Normalized hits on the main strand(s): $print_strandbias%",$black);
								if($best_split>-1)
									{
									$split1=$coordinates_for_split[$best_split];
									$split2=$coordinates_for_split[$best_split+1];
									$directionality.=" (split between $split1 and $split2)";
									}
								$image->string(gdSmallFont,100,680,"Predicted directionality: $directionality",$black);
								if($results_table==1)
									{
									print RESULTS_TEXT"Cluster $id\tLocation: $cl_location\tCoordinates: $start-$end\tSize [bp]: $print_clustersize\tHits (absolute): $hits_abs\tHits (normalized): $hits_norm\tHits (normalized) per kb: $print_density\tNormalized hits with 1T: $print_1T%\tNormalized hits with 10A: $print_10A%\tNormalized hits $min_piRNAsize-$max_piRNAsize nt: $print_size%\tNormalized hits on the main strand(s): $print_strandbias%\tPredicted directionality: $directionality";
									}
								
								$image->rectangle(100,50,600,251,$light_gray);
								$image->rectangle(100,300,600,501,$light_gray);
								$image->fill(101,51,$light_gray);
								$image->fill(101,301,$light_gray);
								$image->line(99,150,600,150,$black);
								$image->line(99,400,600,400,$black);
								$image->line(99,50,99,251,$black);
								$image->line(99,300,99,501,$black);
								$image->line(96,50,103,50,$black);
								$image->line(96,251,103,251,$black);
								
								$image->rectangle(10,301,20,310,$color[0]);$image->fill(11,302,$color[0]);$image->string(gdTinyFont,25,302,"1 hit",$black);
								$image->rectangle(10,311,20,320,$color[1]);$image->fill(11,312,$color[1]);$image->string(gdTinyFont,25,312,"2-5 hits",$black);
								$image->rectangle(10,321,20,330,$color[2]);$image->fill(11,322,$color[2]);$image->string(gdTinyFont,25,322,"6-10 hits",$black);
								$image->rectangle(10,331,20,340,$color[3]);$image->fill(11,332,$color[3]);$image->string(gdTinyFont,25,332,"11-20 hits",$black);
								$image->rectangle(10,341,20,350,$color[4]);$image->fill(11,342,$color[4]);$image->string(gdTinyFont,25,342,"21-50 hits",$black);
								$image->rectangle(10,351,20,360,$color[5]);$image->fill(11,352,$color[5]);$image->string(gdTinyFont,25,352,"51-100 hits",$black);
								$image->rectangle(10,361,20,370,$color[6]);$image->fill(11,362,$color[6]);$image->string(gdTinyFont,25,362,"> 100 hits",$black);
								
								$cluster_plus_flank=$end-$start+1+($flank_size*2);
								if($best_split>-1)
									{
									$split_coordinate=(($split1+$split2)/2)-$start+$flank_size;
									$image->line(($split_coordinate/$cluster_plus_flank)*500+100,50,($split_coordinate/$cluster_plus_flank)*500+100,251,$gray);
									$image->line(($split_coordinate/$cluster_plus_flank)*500+100,300,($split_coordinate/$cluster_plus_flank)*500+100,501,$gray);
									
									if($directionality=~/plus-minus/)
										{
										$image->fill(101,51,$light_blue);
										$image->fill(599,250,$light_red);
										$image->fill(101,301,$light_blue);
										$image->fill(599,500,$light_red);
										}
									elsif($directionality=~/minus-plus/)
										{
										$image->fill(599,51,$light_blue);
										$image->fill(101,250,$light_red);
										$image->fill(599,301,$light_blue);
										$image->fill(101,500,$light_red);
										}
									}
								elsif($directionality eq"mono:plus")
									{
									$image->fill(101,51,$light_blue);
									$image->fill(101,301,$light_blue);
									}
								elsif($directionality eq"mono:minus")
									{
									$image->fill(104,250,$light_red);
									$image->fill(104,500,$light_red);
									}
									
								if($flank_size>0)
									{
									$image->line(($flank_size/$cluster_plus_flank*500)+99,50,($flank_size/$cluster_plus_flank*500)+99,251,$black);
									$image->line(600-($flank_size/$cluster_plus_flank*500),50,600-($flank_size/$cluster_plus_flank*500),251,$black);
									
									$image->line(($flank_size/$cluster_plus_flank*500)+99,300,($flank_size/$cluster_plus_flank*500)+99,501,$black);
									$image->line(600-($flank_size/$cluster_plus_flank*500),300,600-($flank_size/$cluster_plus_flank*500),501,$black);
									
									$image->fill(599,51,$dark_gray);
									$image->fill(599,250,$dark_gray);
									$image->fill(599,301,$dark_gray);
									$image->fill(599,500,$dark_gray);
									
									$image->fill(101,51,$dark_gray);
									$image->fill(101,250,$dark_gray);
									$image->fill(101,301,$dark_gray);
									$image->fill(101,500,$dark_gray);
									}
								
								open(IN,"$out_folder/$id.fasta");
								%transcription_plus=();
								%transcription_minus=();
								$extreme=0;
								while(<IN>)
									{
									if($_=~/^>/)
										{
										@d=split("\t",$_);
										if($d[5]=~/\+/)
											{
											if($normalize_by_total_number_of_mapped_reads==0)
												{
												foreach$pos($d[1]..($d[1]+length$d[1])-1)
													{
													$transcription_plus{$pos}+=$d[4];
													if($transcription_plus{$pos}>$extreme)
														{
														$extreme=$transcription_plus{$pos};
														}
													}
												}
											# calculate values per million mapped reads (rpm)
											else
												{
												foreach$pos($d[1]..($d[1]+length$d[1])-1)
													{
													$transcription_plus{$pos}+=($d[4]/$total_reads)*1000000;
													if($transcription_plus{$pos}>$extreme)
														{
														$extreme=$transcription_plus{$pos};
														}
													}
												}
												
											}
										else
											{
											if($normalize_by_total_number_of_mapped_reads==0)
												{
												foreach$pos($d[1]..($d[1]+length$d[1])-1)
													{
													$transcription_minus{$pos}+=$d[4];
													if($transcription_minus{$pos}>$extreme)
														{
														$extreme=$transcription_minus{$pos};
														}
													}
												}
											# calculate values per million mapped reads (rpm)
											else
												{
												foreach$pos($d[1]..($d[1]+length$d[1])-1)
													{
													$transcription_minus{$pos}+=(($d[4])/$total_reads)*1000000;
													if($transcription_minus{$pos}>$extreme)
														{
														$extreme=$transcription_minus{$pos};
														}
													}
												}
											}
										
										}
									}
								close IN;
								$print_extreme=(int(($extreme*100)+0.5))/100;
								$length_extremestring=length$print_extreme;
								$extreme_poscorr=$length_extremestring*5;
								
								$image->string(gdTinyFont,95-$extreme_poscorr,47,$print_extreme,$black);
								$image->string(gdTinyFont,95-$extreme_poscorr,247,$print_extreme,$black);

								foreach$pos($start-$flank_size..$end+$flank_size)
									{
									$relative_position=$pos-$start+1+$flank_size;
									if($transcription_plus{$pos})
										{
										$image->line(($relative_position/$cluster_plus_flank)*500+100,149,($relative_position/$cluster_plus_flank)*500+100,int((149-(($transcription_plus{$pos}/$extreme)*100))+0.5),$blue);
										}
									if($transcription_minus{$pos})
										{
										$image->line(($relative_position/$cluster_plus_flank)*500+100,151,($relative_position/$cluster_plus_flank)*500+100,int((151+(($transcription_minus{$pos}/$extreme)*100))+0.5),$red);
										}
									}
								undef%transcription_plus;
								undef%transcription_minus;
								
								# search binding motifs
								if($search_bindingsites==1)
									{
									%sites_in_cluster=();
									$sites_in_cluster="";
									$found_motifs=0;
									$found_motifs_rc=0;
									$image->string(gdSmallFont,605,35,"binding sites (motif: position)",$black);
									foreach$motif_name(keys%binding_motifs)
										{
										$check_sequence=$cluster_sequence;
										while(1)
											{
											if($check_sequence=~/$binding_motifs{$motif_name}/)
												{
												$hit_motif=$&;
												$position=index($check_sequence,$hit_motif);
												$position++;
												$replace="";
												foreach(1..length$hit_motif)
													{
													$replace.="X";
													}
												$check_sequence=~s/$binding_motifs{$motif_name}//;
												$check_sequence=$replace.$check_sequence;
												$print_motiv_name=$motif_name;
												$print_motiv_name=~s/ rc$//;
												if($motif_name=~/ rc$/)
													{
													unless($sites_in_cluster{"$hit_motif-$position"})
														{
														$found_motifs_rc++;
														$sites_in_cluster{"$hit_motif-$position"}=1;
														$image->rectangle((($position/$cluster_plus_flank)*500)+99,149,(($position/$cluster_plus_flank)*500)+101,151,$color[0]);
														$image->line((($position/$cluster_plus_flank)*500)+100,149,(($position/$cluster_plus_flank)*500)+100,256-(10*$found_motifs_rc),$color[0]);
														$image->line((($position/$cluster_plus_flank)*500)+100,256-(10*$found_motifs_rc),602,256-(10*$found_motifs_rc),$color[0]);
														$position=$position-$flank_size;
														$image->string(gdTinyFont,605,252-(10*$found_motifs_rc),"$print_motiv_name ($hit_motif: $position)",$black);
														$sites_in_cluster.="$print_motiv_name ($hit_motif: $position) ";
														}
													}
												else
													{
													unless($sites_in_cluster{"$hit_motif-$position"})
														{
														$found_motifs++;
														$sites_in_cluster{"$hit_motif-$position"}=1;
														$image->rectangle((($position/length$cluster_sequence)*500)+99,149,(($position/length$cluster_sequence)*500)+101,151,$color[0]);
														$image->line((($position/length$cluster_sequence)*500)+100,149,(($position/length$cluster_sequence)*500)+100,44+(10*$found_motifs),$color[0]);
														$image->line((($position/length$cluster_sequence)*500)+100,44+(10*$found_motifs),602,44+(10*$found_motifs),$color[0]);
														$position=$position-$flank_size;
														$image->string(gdTinyFont,605,40+(10*$found_motifs),"$print_motiv_name ($hit_motif: $position)",$black);
														$sites_in_cluster.="$print_motiv_name ($hit_motif: $position) ";
														$sites_in_cluster{"$hit_motif-$position"}=1;
														}
													}
												}
											else
												{
												last;
												}
											}
										}
									if($results_table==1&&$found_motifs+$found_motifs_rc>0)
										{
										$sites_in_cluster=~s/ $//;
										print RESULTS_TEXT"\tBinding sites: $sites_in_cluster\n";
										}
									else
										{
										print RESULTS_TEXT"\n";
										}
									undef%sites_in_cluster;
									}
								elsif($results_table==1)
									{
									print RESULTS_TEXT"\n";
									}
								
								@redundancy_code=("1-1","2-5","6-10","11-20","21-50","51-100","101-999999999");
								foreach$redundancy(0..6)
									{
									open(IN,"$out_folder/$id.fasta");
									$color_id=6-$redundancy;
									$range=pop@redundancy_code;
									@range=split('-',$range);
									while(<IN>)
										{
										if($_=~/^>/)
											{
											@d=split("\t",$_);
											if($d[3]>=$range[0]&&$d[3]<=$range[1])
												{
												$relative_position=$d[1]-$start+1+$flank_size;
												if($d[5]=~/\+/)
													{
													$image->line(($relative_position/$cluster_plus_flank)*500+100,399,($relative_position/$cluster_plus_flank)*500+100,390,$color[$color_id]);
													}
												else
													{
													$image->line(($relative_position/$cluster_plus_flank)*500+100,401,($relative_position/$cluster_plus_flank)*500+100,410,$color[$color_id]);
													}
												}
											}
										}
									close IN;
									}
								if($image_files==1)
									{
									open(GD,">$out_folder/$id.png");
									binmode GD;
									print GD $image->png;
									close GD;
									}
								}
							if($fasta_piRNA_files==0)
								{
								unlink"$out_folder/$id.fasta";
								}
							}
						}
					undef%seqs_in_candidate;
					} 
				}
			$stat=0;
			}
		}
	}
close OUT;
close MAP;
$percent_clusters_to_genome=(int((($total_size_of_all_clusters/$genome_size)*100000)+0.5))/1000;
$total_clustered_sequences=keys%total_clustered_sequences;
$total_clustered_reads=0;
foreach(keys%total_clustered_sequences)
	{
	$total_clustered_reads+=$reads_per_seq{$_};
	}
$percent_clustered_reads=(int((($total_clustered_reads/$total_reads)*100000)+0.5))/1000;
$percent_clustered_to_all_seq=(int((($total_clustered_sequences/$non_ID)*100000)+0.5))/1000;
if($results_table==1)
	{
	print "\n$exp_or_obs_text";  #kw_mod
	print RESULTS_TEXT"\n$exp_or_obs_text\n";  #kw_mod
	print RESULTS_TEXT"Non-identical sequences: $non_ID\n"; #kw_mod
	print RESULTS_TEXT"Total Clustered Sequences: $total_clustered_sequences\n"; #kw_mod
	print RESULTS_TEXT"Total Clustered Reads: $total_clustered_reads\n"; #kw_mod
	print RESULTS_TEXT"Total Reads: $total_reads\n"; #kw_mod

	print RESULTS_TEXT"\nTotal size of $id predicted piRNA clusters: $total_size_of_all_clusters bp ($percent_clusters_to_genome%)\nNon identical sequences that can be assigned to clusters: $total_clustered_sequences ($percent_clustered_to_all_seq%)\nSequence reads that can be assigned to clusters: $total_clustered_reads ($percent_clustered_reads%)\n";
	close RESULTS_TEXT;
	}
print"\n\nNon-identical sequences: $non_ID\n"; #kw_mod
print"Total Clustered Sequences: $total_clustered_sequences\n"; #kw_mod
print"Total Clustered Reads: $total_clustered_reads\n"; #kw_mod
print"Total Reads: $total_reads\n"; #kw_mod
print"\n\nTotal size of $id predicted piRNA clusters: $total_size_of_all_clusters bp ($percent_clusters_to_genome%)\nNon identical sequences that can be assigned to clusters: $total_clustered_sequences ($percent_clustered_to_all_seq%)\nSequence reads that can be assigned to clusters: $total_clustered_reads ($percent_clustered_reads%)\n\n";
exit;
