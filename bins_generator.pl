use strict;
use Cwd;
use Getopt::Long;

my $help;
my $type;
my $mode;
my $length;
my $nbins;
my $ref;

GetOptions(
  'h|help' => \$help,
  'm|mode=s' => \$mode,
  't|type=s' => \$type,
  'l|length=i' => \$length,
  'n|nbins=s' => \$nbins,
  'r|reference=s' => \$ref,

);

if($help || !defined $type) {
  print_help();
  exit(0);
}

if($type eq "GB") {
  unless($ref) {
    warn "please provide a valide gene reference";
    print_help();
    exit(0);
  }
  
  if(!-f $ref) {
    die "error $ref: no such file or directory\n";
  }
  
  unless($nbins) {
    warn "NBINS is set as 50\n";
    $nbins = 50;
  }
  
  &GeneBins($ref, $nbins)
} elsif($type eq "GI") {
  unless($ref) {
    warn "please provide a valide gene reference";
    print_help();
    exit(0);
  }
  
  if(!-f $ref) {
    die "error $ref: no such file or directory\n";
  }
  
  &GetGenesInfo($ref);
} elsif($type eq "TB") {
  unless($ref) {
    warn "please provide a valide gene reference";
    print_help();
    exit(0);
  }
  
  if(!-f $ref) {
    die "error $ref: no such file or directory\n";
  }
  
  unless($nbins) {
    warn "NBINS is set as 50\n";
    $nbins = 50;
  }
  
  unless($length) {
    warn "LENGTH is set as +/-5000bp\n";
    $length = 5000;
  }
  
  &GetTSSBins($ref, $length, $nbins);
} elsif($type eq "RB") {
  if(@ARGV<=0 || @ARGV>=2) {
    print_help();
    exit();
  }
  
  unless($nbins) {
    warn "NBINS is set as 50\n";
    $nbins = 50;
  }
  
  unless($length) {
    warn "LENGTH is set as +/-5000bp\n";
    $length = 5000;
  }
  
  unless($mode) {
    warn "MODE is set as CENTER\n";
    $mode = "center";
  }
  
  &GetRegionBins($ARGV[0], $length, $nbins, $mode)
}

sub GeneBins
{
	my $genesFile = shift;
  
  my $pstart;
  my $len;
  my $count;
  my $n = shift;
	
	open(IN, "<$genesFile") or die "$!\n";
	<IN>; # omit header. 
  open(OUT, "| gzip -c >gene.tables.txt.gz") or die "$!\n";
  print OUT "chrom\tstart\tend\tfeature\tbin\tgene_id\n";
  
	while(<IN>)
	{
		chomp;
		my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds, $gene_id) = (split("\t"))[1..10,12];
		
		next if($chrom =~ /_random/i || $chrom =~ /chrX/i || $chrom =~ /chrY/i || $chrom =~ /chrM/i || $chrom =~ /chrUn/i);
		
		if($strand eq "+")
		{
      $pstart = $txStart - 10000 + 1; $len = 8000 / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Distal_promoter", $i, $gene_id), "\n";
      }
      
      $pstart = $txStart - 2000 + 1; $len = 2000 / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Promoter", $i, $gene_id), "\n";
      }
      
      $pstart = $txStart + 1; $len = ($cdsStart - $txStart) / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "5pUTR", $i, $gene_id), "\n";
      }
      
      $pstart = $cdsEnd + 1; $len = ($txEnd - $cdsEnd) / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "3pUTR", $i, $gene_id), "\n";
      }
  	}
  	elsif($strand eq "-")
  	{
      $pstart = $txEnd + 2000 + 1; $len = 8000 / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Distal_promoter", 20 - $i + 1, $gene_id), "\n";
      }
      
      $pstart = $txEnd + 1; $len = 2000 / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Promoter", 20 - $i + 1, $gene_id), "\n";
      }
      
      $pstart = $cdsEnd + 1; $len = ($txEnd - $cdsEnd) / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "5pUTR", 20 - $i + 1, $gene_id), "\n";
      }
      
      $pstart = $txStart + 1; $len = ($cdsStart - $txStart) / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "3pUTR", 20 - $i + 1, $gene_id), "\n";
      }
  	}
    		
		my @exStart = split(",", $exonStarts);
		my @exEnd = split(",", $exonEnds);
		for(my $i = 0; $i < $exonCount; $i++)
		{
      $pstart = $exStart[$i] + 1; $len = ($exEnd[$i] - $exStart[$i]) / $n;
      for(my $i = 1; $i <= $n; $i ++) {
        print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Exon", $i, $gene_id), "\n";
      }
      
      if($i > 0) {
        $pstart = $exEnd[$i - 1] + 1; $len = ($exStart[$i] - $exEnd[$i - 1]) / $n;
        for(my $i = 1; $i <= $n; $i ++) {
          print OUT join("\t", $chrom, int($pstart + ($i - 1) * $len), int($pstart + $i * $len - 1), "Intron", $i, $gene_id), "\n";
        }
      }
		}
	}
	close(IN);
  close(OUT);
}

sub GetGenesInfo
{
  my $annotation = 0;
  $annotation = 1 if @_ > 1;
	my $genesFile = shift;
  
	my %pomoter = ();
	my %_5UTR = ();
	my %_3UTR = ();
	my %exon = ();
	my %intron = ();
  my %intergenic = ();
  my %genes = ();
  my %distal_promoter = ();
	
	open(IN, "<$genesFile") or die "$!\n";
	<IN>; # omit header. 
	while(<IN>)
	{
		chomp;
		my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds, $gene_id) = (split("\t"))[1..10,12];
		
		next if($chrom =~ /_random/i || $chrom =~ /chrX/i || $chrom =~ /chrY/i || $chrom =~ /chrM/i || $chrom =~ /chrUn/i);
		
		if($strand eq "+")
		{
			push(@{$pomoter{$chrom}}, [$txStart - 2000 > 0 ? $txStart - 2000 : 0, $txStart - 1, $gene_id]);
      push(@{$distal_promoter{$chrom}}, [$txStart - 10000 > 0 ? $txStart - 10000 : 0, $txStart - 2001, $gene_id]);
			push(@{$_5UTR{$chrom}}, [$txStart, $cdsStart, $gene_id]);
			push(@{$_3UTR{$chrom}}, [$cdsEnd, $txEnd, $gene_id]);
  	}
  	elsif($strand eq "-")
  	{
  		push(@{$pomoter{$chrom}}, [$txEnd + 1,  $txEnd + 2000, $gene_id]);
      push(@{$distal_promoter{$chrom}}, [$txEnd + 2001,  $txEnd + 10000, $gene_id]);
			push(@{$_3UTR{$chrom}}, [$txStart, $cdsStart, $gene_id]);
			push(@{$_5UTR{$chrom}}, [$cdsEnd, $txEnd, $gene_id]);
  	}
    
    push(@{$genes{$chrom}}, [$txStart, $txEnd, $gene_id]);
		
		my @exStart = split(",", $exonStarts);
		my @exEnd = split(",", $exonEnds);
		for(my $i = 0; $i < $exonCount; $i++)
		{
			push(@{$exon{$chrom}}, [$exStart[$i], $exEnd[$i], $gene_id]);
			push(@{$intron{$chrom}}, [$exEnd[$i - 1], $exStart[$i], $gene_id]) unless $i == 0;
		}
	}
	close(IN);
	
  foreach my $key (keys %genes)
  {
    @{$genes{$key}} = sort {$a->[0] <=> $b->[0]} @{$genes{$key}};
    
    my @record = ();
    my $current_start = ${$genes{$key}}[1][0];
    my $current_end = ${$genes{$key}}[1][1];
    for(my $i=2; $i <= $#{$genes{$key}}; $i ++)
    {
      if(${$genes{$key}}[$i][0] <= $current_end)
      {
        $current_end = $current_end < ${$genes{$key}}[$i][1] ? ${$genes{$key}}[$i][1] : $current_end;
        
        if($i == $#{$genes{$key}})
        {
          push(@record, [$current_start, $current_end]);
          last;
        }
      }
      else
      {
        push(@record, [$current_start, $current_end]);
        $current_start = ${$genes{$key}}[$i][0];
        $current_end = ${$genes{$key}}[$i][1];
        
        if($i == $#{$genes{$key}})
        {
          push(@record, [$current_start, $current_end]);
          last;
        }
      }
    }
    
    for(my $i=2; $i <= $#record; $i ++)
    {
      push(@{$intergenic{$key}}, [$record[$i-1][1], $record[$i][0]]);
    }
  }
  
  my $curr_dir = getcwd();
  if(! -d "$curr_dir/table")
  { 
    mkdir "$curr_dir/table";
  }
  
  open(Promoter, ">$curr_dir/table/promoter.txt") or die "cannot create$curr_dir/table/promoter.txt";
  open(DistalPromoter, ">$curr_dir/table/distal.promoter.txt") or die "cannot create$curr_dir/table/distal.promoter.txt";
  open(UTR5, ">$curr_dir/table/5pUTR.txt") or die "cannot create $curr_dir/table/5pUTR.txt";
  open(UTR3, ">$curr_dir/table/3pUTR.txt") or die "cannot create $curr_dir/table/3pUTR.txt";
  open(Exon, ">$curr_dir/table/exon.txt") or die "cannot create $curr_dir/table/exon.txt";
  open(Intron, ">$curr_dir/table/intron.txt") or die "cannot create $curr_dir/table/intron.txt";
  open(intergenic, ">$curr_dir/table/intergenic.txt") or die "cannot create $curr_dir/table/intergenic.txt";
  
  print Promoter "chr\tstart\tend\n" if !$annotation;
  print DistalPromoter "chr\tstart\tend\n" if !$annotation;
  print UTR5 "chr\tstart\tend\n" if !$annotation;
  print UTR3 "chr\tstart\tend\n" if !$annotation;
  print Exon "chr\tstart\tend\n" if !$annotation;
  print Intron "chr\tstart\tend\n" if !$annotation;
  print intergenic "chr\tstart\tend\n" if !$annotation;
  
  print Promoter "chr\tstart\tend\tgene_id\n" if $annotation;
  print DistalPromoter "chr\tstart\tend\tgene_id\n" if $annotation;
  print UTR5 "chr\tstart\tend\tgene_id\n" if $annotation;
  print UTR3 "chr\tstart\tend\tgene_id\n" if $annotation;
  print Exon "chr\tstart\tend\tgene_id\n" if $annotation;
  print Intron "chr\tstart\tend\tgene_id\n" if $annotation;
  print intergenic "chr\tstart\tend\tgene_id\n" if $annotation;

	foreach my $key (keys %pomoter)
	{
    foreach my $v (@{$pomoter{$key}}) 
    {
      print Promoter "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print Promoter "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$distal_promoter{$key}}) 
    {
      print DistalPromoter "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print DistalPromoter "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$_5UTR{$key}}) 
    {
      print UTR5 "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print UTR5 "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$_3UTR{$key}}) 
    {
      print UTR3 "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print UTR3 "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$exon{$key}}) 
    {
      print Exon "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print Exon "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$intron{$key}}) 
    {
      print Intron "$key\t$$v[0]\t$$v[1]\n" if !$annotation;
      print Intron "$key\t$$v[0]\t$$v[1]\t$$v[2]\n" if $annotation;
    }
    
    foreach my $v (@{$intergenic{$key}}) 
    {
      print intergenic "$key\t$$v[0]\t$$v[1]\n";
    }
	}
  
  close(Promoter); close(DistalPromoter); close(UTR5); close(UTR3); close(Exon); close(Intron); close(intergenic); 
}

sub GetTSSBins
{
	my $genesFile = shift;
  my $length = shift;
  my $num_bins = shift;
	my %bins = ();
  
  my $step = $length / $num_bins;
	
	open(IN, "<$genesFile") or die "$!\n";
	<IN>; # omit header. 
	while(<IN>)
	{
		chomp;
		my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds) = (split("\t"))[1..10];
		
		next if($chrom =~ /_random/i || $chrom =~ /chrX/i || $chrom =~ /chrY/i || $chrom =~ /chrM/i || $chrom =~ /chrUn/i);
		
		if($strand eq "+")
		{
			for(my $i = 1; $i <= 2 * $num_bins; $i ++)
      {
        my $start = $txStart - $length + 1 + $i * $step; 
        my $end = $txStart - $length + 1 + ($i + 1) * $step;
        push(@{${$bins{join("_", "bin", $i)}}{$chrom}}, [$start, $end]);
      }
  	}
  	elsif($strand eq "-")
  	{
  		for(my $i = 1; $i <= 2 * $num_bins; $i ++)
      {
        my $start = $txEnd - $length + 1 + $i * $step; 
        my $end = $txEnd - $length + 1 + ($i + 1) * $step;
        push(@{${$bins{join("_", "bin", $i)}}{$chrom}}, [$start, $end]);
      }
  	}
	}
	close(IN);
	
  my $curr_dir = getcwd();
  
  if(! -d "$curr_dir/table")
  { 
    mkdir "$curr_dir/table";
  }
  
  foreach my $bin (keys %bins)
  {
    open(BIN, ">$curr_dir/table/$bin\.txt") or die "$!\n";
    print BIN "chr\tstart\tend\n";
    foreach my $chrom (keys %{$bins{$bin}})
    {
      foreach my $v (@{${$bins{$bin}}{$chrom}})
      {
        print BIN "$chrom\t$$v[0]\t$$v[1]\n";
      }
    }
    close(BIN);
  }
}

sub GetRegionBins
{
	my $peakFile = shift;
  my $length = shift;
  my $num_bins = shift;
  my $mode = shift;
	my %bins = ();
  
  my $step;
  my $start;
  my $st;
  my $en;
	
	open(IN, "<$peakFile") or die "$peakFile: $!\n";
	<IN>; # omit header. 
	while(<IN>)
	{
		chomp;
    
		my ($chrom, $chStart, $chEnd) = (split("\t"))[0..2];
    my $logic = $chEnd > $chStart;
		
    if($mode =~ /center/i) {
      $step = 2*$length / (2*$num_bins+1);
      
      $chStart = int(($chStart + $chEnd) / 2);
      for(my $i = 1; $i <= 2 * $num_bins + 1; $i ++)
      {
        $st = int($chStart - $length + 1 + ($i - 1) * $step); 
        $en = int($chStart - $length + $i * $step);
        
        if($logic) {push(@{$bins{join("_", "bin", $i)}}, [$chrom, $st, $en]);}
        else {push(@{$bins{join("_", "bin", 2 * $num_bins - $i + 2)}}, [$chrom, $st, $en]);}
      }
    }
    elsif($mode =~ /region/i){
      if($logic) {
        $start = $chStart;
        $step = ($chEnd - $chStart + 1) / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          push(@{$bins{join("_", "bin_mid", $i)}}, [$chrom, $st, $en]);
        }
        
        $start = $chStart - $length + 1;
        $step = $length / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          push(@{$bins{join("_", "bin_left", $i)}}, [$chrom, $st, $en]);
        }
        
        $start = $chEnd;
        $step = $length / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          push(@{$bins{join("_", "bin_right", $i)}}, [$chrom, $st, $en]);
        }
      }
      else {
        $start = $chEnd;
        $step = ($chStart - $chEnd + 1) / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          if($st <= $en) {push(@{$bins{join("_", "bin_mid", $num_bins - $i + 1)}}, [$chrom, $st, $en]);}
          else {push(@{$bins{join("_", "bin_mid", $num_bins - $i + 1)}}, [$chrom, $en, $st]);}
        }
        
        $start = $chStart;
        $step = $length / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          if($st <= $en) {push(@{$bins{join("_", "bin_left", $num_bins - $i + 1)}}, [$chrom, $st, $en]);}
          else {push(@{$bins{join("_", "bin_left", $num_bins - $i + 1)}}, [$chrom, $en, $st]);}
        }
        
        $start = $chEnd - $length + 1;
        $step = $length / $num_bins;
        for(my $i = 1; $i <= $num_bins; $i ++)
        {
          $st = int($start + ($i - 1) * $step); 
          $en = int($start + $i * $step - 1);
          if($st <= $en) {push(@{$bins{join("_", "bin_right", $num_bins - $i + 1)}}, [$chrom, $st, $en]);}
          else {push(@{$bins{join("_", "bin_right", $num_bins - $i + 1)}}, [$chrom, $en, $st]);}
        }
      }
    }
    
	}
	close(IN);
	
  $peakFile =~ s/.txt$//;
  $peakFile =~ s/.bed$//;
  
  open(BIN, ">${peakFile}.bins.txt") or die "$!\n";
  print BIN "chr\tstart\tend\tbin_info\n";
  foreach my $bin (sort{$a cmp $b} keys %bins)
  {
    foreach my $v (@{$bins{$bin}})
    {
      print BIN "$$v[0]\t$$v[1]\t$$v[2]\t$bin\n";
    }
  }
  close(BIN);
}

sub print_help() {
  warn "perl bins_generator.pl [-m/--mode center(default)|region] [-l/--length <int>] [-n/--nbins <int>] [-r/--reference <string>] <-t/--type GB|GI|TB|RB> [Input]";
  warn "where LENGTH is optional for some of functions: LENGTH>0\n";
  warn "where NBINS is also optional for some of functions: NBINS>0\n";
  warn "where REFERENCE denotes gene reference (such as: refseq.genes) which inlucde all annotations (more than 12 columns) for genes\n";
  warn "where TYPE is an argument required, which denotes which function will be executed; optional types inlcude: to divide gene-structured (Promoter, Exon, Intron, ...) into bins (GB); to divide the given regions around TSS into bins (TB); to get gene features (GI); to divide region and its upstream/downstream into bins (RB).\n"
}
__END__