#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin", "$FindBin::Bin/KmerGraphLib");

use Fasta_reader;

use ColorGradient;

use KmerNode;
use KmerGraph;
use Carp;

use Getopt::Long qw(:config no_ignore_case bundling);

no warnings 'recursion';

my $usage = <<_EOUSAGE_;

#####################################################################################
#
#  Required:
#
#  --fasta          fasta files containing reads (comma-separated list of files:  fileA.reads,fileB.reads,...) 
#
#  and:
#       --assemble
#    and/or
#       --graph [filename]    dot output file for viewing structure in GraphViz
#
#  Optional:
#  -K               Kmer length (default: 25)
#  -R               maximum recursive search depth from terminus for extension (default: 5) # note can exponentially increases runtime
#
#
#  -L                    minimum assembled sequence length to report (default: 100)
#  -E                    minimum Kmer seed entropy (default: 1.5;  note maximum entropy = 2 for nucleotide sequences)
#  --min_coverage        minimum kmer coverage to report kmer from reads in graph.
#
#
#  -v                    verbose
#
#  --CDS             reference CDS sequence to thread through the graph
#  --iworm           inchworm greedy path results (fasta)
#  --min_iworm_coverage   min average kmer coverage for inchworm assembly to add to the graph.
#  
#  --no_compact          do NOT compact the graph
#
#  --min_edge_threshold    default(0.1)
#  --min_leaf_node_length  default(76)
#  --min_leaf_node_avg_cov default(1.2)
#
#  -N                    max assemblies to report (default: infinity)
#
#####################################################################################

_EOUSAGE_

	;



my $help_flag;

my $fasta_file;
my $KmerLength = 25;

my $max_recurse_depth = 5;

my $min_length = 100;
my $min_cov = 1;

our $VERBOSE =0;

my $MIN_KMER_ENTROPY = 1.5;
my $graph_file;
my $assemble_flag;
my $cds_file;
my $iworm_file;
my $min_iworm_coverage = 0;
my $MIN_EDGE_THRESHOLD = 0.1;
my $MIN_LEAF_NODE_LENGTH = 76;
my $MIN_LEAF_NODE_AVG_COV = 1.2;


my $NO_COMPACT = 0;
my $max_asmbls_report = -1;

&GetOptions ( 'h' => \$help_flag,
                          
			  ## required
			  'fasta=s' => \$fasta_file,
			  'K=i' => \$KmerLength,
			  'R=i' => \$max_recurse_depth,
			  
			  ## options
			  'L=i' => \$min_length,
			  'E=f' => \$MIN_KMER_ENTROPY,
			  "min_coverage=i" => \$min_cov,
			  "graph=s" => \$graph_file,
			  "assemble" => \$assemble_flag,
			  "iworm=s" => \$iworm_file,
			  'min_iworm_coverage=i' => \$min_iworm_coverage,
			  "CDS=s" => \$cds_file,

			  "min_edge_threshold=f" => \$MIN_EDGE_THRESHOLD,
			  "min_leaf_node_length=i" => \$MIN_LEAF_NODE_LENGTH,
			  "min_leaf_node_avg_cov=f" =>  \$MIN_LEAF_NODE_AVG_COV,
			  
			  "no_compact" => \$NO_COMPACT,
			  'N=i' => \$max_asmbls_report,

			  # hidden option for debugging
			  'v' => \$VERBOSE,
			  

	);

if ($help_flag) {
	die $usage;
}

unless ($fasta_file && ($graph_file || $assemble_flag)) {
	die $usage;
}


our $SEE = $VERBOSE;


main: {
	
	my $KmerGraph = KmerGraph->new($KmerLength);

	my $fasta_reader = new Fasta_reader($fasta_file);
	
	my $seq_counter = 0;
	
	while (my $seq_obj = $fasta_reader->next()) {
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();

		$seq_counter++;
		
		print STDERR "\r-parsing seq($seq_counter)        ";
		
		$KmerGraph->add_sequence_to_graph($acc, $sequence, 1, 'black');
		
	}

	if ($cds_file) {
		my $fasta_reader = new Fasta_reader($cds_file);
		while (my $seq_obj = $fasta_reader->next()) {
			my $seq = $seq_obj->get_sequence();
			my $acc = $seq_obj->get_accession();
			print STDERR "-adding reference sequence $acc from $cds_file\n";
			$KmerGraph->add_sequence_to_graph($acc, $seq, 0, 'red');
		}
	}

	if ($iworm_file) {
		# my $seq = &get_highest_scoring_iworm_assembly($iworm_file);
		
		my $fr = new Fasta_reader($iworm_file);

		my @seqs;

		while (my $seq_obj = $fr->next()) {
			my $seq = $seq_obj->get_sequence();
			my $acc = $seq_obj->get_accession();
			
			if (length($seq) >= $min_length) {
				push (@seqs, [$acc, $seq]);
			}
			
		}
		
		my @colors;
		if (scalar @seqs > 1) {
			@colors = &ColorGradient::get_RGB_gradient(scalar @seqs);
		
			@colors = &ColorGradient::convert_RGB_hex(@colors);
		}
		else {
			@colors = ('green');
		}
		
		while (@seqs) {
			my $seq_info = shift @seqs;
			my ($acc, $seq) = @$seq_info;
			my $color = shift @colors;
			print STDERR "Adding additional seq trace: $acc, $color\n";
			$KmerGraph->add_sequence_to_graph($acc, $seq, 0, $color);
			
		}
	}
	

	
	if ($min_cov > 1) {
		print STDERR "\n-removing kmer nodes below coverage: $min_cov\n";
		
		$KmerGraph->purge_nodes_below_count($min_cov);
	}
	


	print STDERR "\nDone adding sequences\n";



	unless ($NO_COMPACT) {
		
		my $round = 0;

		my ($compacted_graph_flag, $number_dangling_nodes_pruned, $number_pruned_edges);
		
		do {


			$compacted_graph_flag = $KmerGraph->compact_graph();
						   
			$number_dangling_nodes_pruned = $KmerGraph->prune_dangling_nodes(min_leaf_node_length => $MIN_LEAF_NODE_LENGTH,
																				min_leaf_node_avg_cov => $MIN_LEAF_NODE_AVG_COV);

			$number_pruned_edges = $KmerGraph->prune_low_weight_edges(edge_weight_threshold => $MIN_EDGE_THRESHOLD);
			
			$round++;
			
			# no op
			print STDERR "Cycling through graph compaction, round: $round\n";
			print STDERR "\tnodes_pruned: $number_dangling_nodes_pruned\n"
				. "\tedges_pruned: $number_pruned_edges\n";
			
			
		} while ($compacted_graph_flag && ($number_dangling_nodes_pruned || $number_pruned_edges));
		
		
		
	}
	
	print STDERR "done.\n";
	
	
	if ($assemble_flag) {

		my @path = $KmerGraph->score_nodes();
		
		my $assembled_seq = $KmerGraph->extract_sequence_from_path(@path);
		
		unless ($assembled_seq) {
			confess "Error, no assembled sequence form path: @path";
		}
		

		my $score = &sum_node_counts(@path);
		
		my $counter = 0;
		if (length($assembled_seq) >= $min_length) {
			$counter++;
			
			my $length_normalized_score = sprintf("%.2f", $score / length($assembled_seq));
			
			print ">assembly_$counter;$score;$length_normalized_score len: " . length($assembled_seq) . "\n$assembled_seq\n";
		}
		
		
		while (! $KmerGraph->all_nodes_visited()) {
			
			if ($max_asmbls_report > 0 && $counter >= $max_asmbls_report) { 
				last;
			}
			
			my @path = $KmerGraph->extract_path_from_unvisited_node(@path);
			
			#$KmerGraph->print_path(@path);
			
			#foreach my $node (@path) {
			#	print $node->toString();
			#}
			
			my $assembled_seq = $KmerGraph->extract_sequence_from_path(@path);  # marks as visited.

			unless ($assembled_seq) {
			    confess "Error, no sequence extracted from path: @path";
			}
			
			if ( length($assembled_seq) >= $min_length) {
				$counter++;
				
				my $score = &sum_node_counts(@path);
				my $length_normalized_score = sprintf("%.2f", $score / length($assembled_seq));
				
				print ">assembly_$counter;$score;$length_normalized_score len: " . length($assembled_seq) . "\n$assembled_seq\n";
			}
		}
	}
	
	## Write graph file in dot format for GraphViz
	if ($graph_file) {
		open (my $ofh, ">$graph_file") or die "Error, cannot write graph file to $graph_file ";
		print $ofh $KmerGraph->toGraphViz( no_short_singletons => 100);
		close $ofh;
	}

		
	
	exit(0);
}


####
sub sum_node_counts {
	my @path = @_;
	
	my $sum = 0;

	foreach my $node (@path) {
		my $count = $node->get_count();
		$sum += $count;
	}

	return($sum);
}

	
####
sub compute_entropy {
	my ($string) = @_;
        
	my @chars = split(//, $string);
        
	my %char_counter;
	foreach my $char (@chars) {
		$char_counter{$char}++;
	}
        
	my $entropy = 0;
        
	my $num_chars = length($string);
	foreach my $char (keys %char_counter) {
		my $count = $char_counter{$char};
		my $prob = $count/$num_chars;
                
		my $val = $prob * log(1/$prob)/log(2);
		$entropy += $val;
	}
        
	return($entropy);
}


####
sub get_highest_scoring_iworm_assembly {
	my ($iworm_file) = @_;
	
	my @entries;

	my $fasta_reader = new Fasta_reader($iworm_file);
	while (my $seq_obj = $fasta_reader->next()) {
		my $seq = $seq_obj->get_sequence();
		my $acc = $seq_obj->get_accession();
		
		my @parts = split(/;/, $acc);
		my $cov = pop @parts;

		my $score = $cov * length($seq);
		
		push (@entries, [$score, $seq]);
	}
	
	@entries = reverse sort {$a->[0]<=>$b->[0]} @entries;
	
	my $best_entry = shift @entries;

	return($best_entry->[1]);
}
