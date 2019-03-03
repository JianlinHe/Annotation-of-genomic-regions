#!/usr/perl/bin

use strict;
use Getopt::Long;

my $help;
my $input;
my $mode;
my $pattern;
my $logic_header;
my $append = -1;
my $depth;
my $abrr;
my $is_avg;

GetOptions(
  'h|help' => \$help,
  'p|pattern=s' => \$pattern,
  'm|mode=s' => \$mode,
  't|title' => \$logic_header,
  'a|append=i' => \$append,
  'd|depth=i' => \$depth,
  'n|name=s' => \$abrr,
  'c|caculate-average' => \$is_avg,
);

if($help || @ARGV < 1){
  warn "no input arguments\n" if !$help && @ARGV < 1;
  warn "Usage: perl to_integrate_methylomes.pl [-h/--help] [-m/--mode] [-p/--pattern] [-t/--title] [-a/--append <int>] [-c/--caculate-average] [-n/--name] [-d/--depth] <inputs...>\n";
  warn "where MODE has two options: merge (default), annotate and format\n";
  warn "where PATTERN is a given pattern which is used in pattern matching. (default: null)\n";
  warn "where NAME is a given name which is used to distinguish different output. (default: none)\n";
  warn "where TITLE is of logical variable: yes (inlcuding a title in input files) or no (default)\n";
  warn "where APPEND is of numeric variable: 1 (multi-ML appending), 2 (multi-counts appending), 0 (default in MERGE, merging) or -1 (default in non-MERGE, not appending). When using MERGE mode, the option will work\n";
  warn "where DEPTH requires the minimal read coverage. (default: 5)\n";
  warn "where CACULATE-AVERAGE is of logical variable:. yes (to calculate average value) or no (default)\n";
  exit(0);
}

unless($mode) {
  $mode = "merge"; 
  warn "default mode: merge\n";
}

unless($pattern) {
  $pattern = "null";
}

if($mode eq "merge" && !$logic_header) {
  warn "in default that there is no title in input files\n";
}

if ($append>=0 && !$depth) {
  $depth = 5;
  warn "in default the minimal read coverage is set as: $depth\n";
}

if($mode eq "merge") {

  unless($append) {
    $append = 0;
    warn "in default: not appending data \n"; 
  }
  &merge(@ARGV);
} elsif($mode eq "annotate") {
  $abrr = "none" unless $abrr;
  
  if(@ARGV < 2) {
    die "MODE in annotation must include two input arguments: object file and annotated file\n";
  } elsif (@ARGV >= 3) {
    warn "the first two arguments will be used\n";
  }
  &annotate($ARGV[0], $ARGV[1]);
} elsif ($mode eq "format") {
  if(@ARGV < 1) {
    die "MODE in format must include two input arguments: an output file from Bismark\n";
  } elsif (@ARGV >= 3){
    warn "the first two arguments will be used\n";
  }
  
  &format_trans($ARGV[0], $ARGV[1]);
}


sub merge {
  my $dir = shift;
  my @files;
  my %common_sites;
  
  opendir(DIR, "$dir") or die "cannot open $dir...\n";
  while(my $file = readdir(DIR)) {
    if($pattern ne "null") {
      if(-f "$dir/$file" && $file =~ m/$pattern/) {
        push(@files, $file);
      }
    }
    else {
      if(-f "$dir/$file") {
        push(@files, $file);
      }
    }
  }
  closedir(DIR);
  
  warn "Dealing with: \n", join("\n", @files), "\n";
  my $head = "";
  foreach my $file (sort{$a cmp $b} @files) {
    open(IN, $file =~ /\.gz$/ ? "gunzip -c $dir/$file |" : "<$dir/$file") or die "cannot open $file ...\n"; ## format: chr position methylation_level #total_CG #methylated_CG
    
    $file =~ s/.txt.gz//;
    $file =~ s/.txt//;
    $file =~ s/_common.sites//;
    
    $head .= "\t$file";
    if($logic_header) {<IN>;}
    
    while(<IN>) {
      chomp;
      my @array = split("\t");
      
      if($append == 1) {
 
        if(!exists $common_sites{"$array[0]\t$array[1]"}) {
          if($array[3]>=$depth) {$common_sites{"$array[0]\t$array[1]"} = "$file:$array[2]";}
          else {$common_sites{"$array[0]\t$array[1]"} = "$file:NA";}
        }
        else {
          if($array[3]>=$depth) {$common_sites{"$array[0]\t$array[1]"} .= "\t$file:$array[2]";}
          else {$common_sites{"$array[0]\t$array[1]"} .= "\t$file:NA";}
        }
      }
      elsif($append == 2) {
 
        if(!exists $common_sites{"$array[0]\t$array[1]"}) {
          if($array[3]>=$depth) {$common_sites{"$array[0]\t$array[1]"} = "$file:$array[3],$array[4]";}
          else {$common_sites{"$array[0]\t$array[1]"} = "$file:NA,NA";}
        }
        else {
          if($array[3]>=$depth) {$common_sites{"$array[0]\t$array[1]"} .= "\t$file:$array[3],$array[4]";}
          else {$common_sites{"$array[0]\t$array[1]"} .= "\t$file:NA,NA";}
        }
      }
      elsif($append == 0) {
      
        if(!exists $common_sites{"$array[0]\t$array[1]"}) {
          $common_sites{"$array[0]\t$array[1]"} = "$array[2]\t$array[3]\t$array[4]";
        }
        else {
          my @tmpArray = split("\t", $common_sites{"$array[0]\t$array[1]"});
          $tmpArray[1] += $array[3]; $tmpArray[2] += $array[4]; 
          if($tmpArray[1] > 0) {
            $tmpArray[0] = $tmpArray[2] / $tmpArray[1];
          }
          else {
            $tmpArray[0] = "NA";
          }
          $common_sites{"$array[0]\t$array[1]"} = join("\t", @tmpArray);
        }
      }      
    }
    close(IN);
  }
  
  if($append == 1) {
  
    open(OUT, "| gzip -c - >$dir/common.CpGs.ge$depth"."X.txt.gz") or die "cannot creat $dir/common.CpGs.ge$depth"."X.txt.gz...\n";
    print OUT "chrom\tposition$head\n";
    
    foreach my $key (sort {$a cmp $b} keys %common_sites) {
      my @array = split("\t", $common_sites{$key});
      my %fh = ();
      foreach my $var (@array) {
        my @tempArray = split(":", $var);
        $fh{$tempArray[0]} = $tempArray[1] unless exists $fh{$tempArray[0]};
      }
      
      my $ml_val = "";
      foreach my $f (sort{$a cmp $b} @files) {
        if(exists $fh{$f}) {
          if($ml_val eq "") {$ml_val = $fh{$f};}
          else {$ml_val .= "\t$fh{$f}";}          
        }
        else {
          if($ml_val eq "") {$ml_val = "NA";}
          else {$ml_val .= "\tNA";}      
        }
      }
      
      print OUT join("\t", $key, $ml_val), "\n";
    }
    close(OUT);
  }
  elsif($append == 2) {
  
    open(OUT, "| gzip -c - >$dir/common.CpGs.ge$depth"."X.txt.gz") or die "cannot creat $dir/common.CpGs.ge$depth"."X.txt.gz...\n";
    my @array = split("\t", $head); shift @array; 
    $head = join("\t", map{$_."_totalCG\t".$_."_mCG"} @array);
    print OUT "chrom\tposition\t$head\n";
    
    foreach my $key (sort {$a cmp $b} keys %common_sites) {
      @array = split("\t", $common_sites{$key});
      my %fh = ();
      foreach my $var (@array) {
        my @tempArray = split(":", $var);
        $fh{$tempArray[0]} = $tempArray[1] unless exists $fh{$tempArray[0]};
      }
      
      my $ml_val = "";
      foreach my $f (sort{$a cmp $b} @files) {
        if(exists $fh{$f}) {
          $fh{$f}=~s/,/\t/g; 
          if($ml_val eq "") {$ml_val = $fh{$f};}
          else {$ml_val .= "\t$fh{$f}";}          
        }
        else {
          if($ml_val eq "") {$ml_val = "NA\tNA";}
          else {$ml_val .= "\tNA\tNA";}      
        }
      }
      
      print OUT join("\t", $key, $ml_val), "\n";
    }
    close(OUT);
  }
  elsif($append == 0) {
  
    open(OUT, "| gzip -c >$dir/${pattern}_merged.txt.gz") or die "cannot creat $dir/${pattern}_merged.txt.gz\n";
    print OUT "chrom\tposition\tavg_ML\ttotal_CG\tme_CG\n";
    foreach my $key (sort {$a cmp $b} keys %common_sites) {
      print OUT "$key\t$common_sites{$key}\n" if split("\t", $common_sites{$key}) == 3;
    }
    close(OUT);
  }
}

sub annotate {
  my ($obj_table, $anno_table) = @_;
  my %anno_list;
  
  open(IN, $anno_table=~/gz$/? "gunzip -c $anno_table | " : "<$anno_table") or die "Cannot open table $anno_table. \n";
  my $anno_header = <IN>;  
  chomp $anno_header;
  $anno_header =~ s/chrom\tposition\t//;
  
  my $elm_num = -1;  
  while(<IN>){
    chomp;
    my @tmpArray = split("\t");
    push(@{$anno_list{$tmpArray[0]}}, [(@tmpArray)[1..$#tmpArray]]);
    
    if($elm_num < 0) {
      $elm_num = @tmpArray - 2;
    }
  }
  close(IN);

  foreach my $chrom (keys %anno_list){
    @{$anno_list{$chrom}} = sort {$a -> [0] <=> $b -> [0]} @{$anno_list{$chrom}};
  }
  
  open(IN, $obj_table=~/gz$/? "gunzip -c $obj_table | " :  "<$obj_table") or die "Cannot open table $obj_table. \n";
  my $obj_header = <IN>; 
  chomp $obj_header;
  
  $obj_table =~ s/.txt//;
  $obj_table =~ s/.bed//;
  
  if($abrr ne "none") {
    open(OUT, "| gzip -c > ${obj_table}.${abrr}.anno.txt.gz") or die "connot create ${obj_table}.${abrr}.anno.txt.gz\n";
    #open(OUT, "> ${obj_table}.${abrr}.anno.txt") or die "connot create ${obj_table}.${abrr}.anno.txt\n";
    print OUT join("\t", $obj_header, $anno_header), "\n";
  }
  else {
    open(OUT, "| gzip -c > ${obj_table}.anno.txt.gz") or die "connot create ${obj_table}.anno.txt.gz\n";
    #open(OUT, "> ${obj_table}.anno.txt") or die "connot create ${obj_table}.anno.txt\n";
    print OUT join("\t", $obj_header, $anno_header), "\n";
  }
  
  while(<IN>) {
    chomp;
    my @tmpArray = split("\t");
    my @array = ();
    
    &binary_search(\@{$anno_list{$tmpArray[0]}}, $tmpArray[1], $tmpArray[2], \@array);
    if(@array > 0 && $elm_num == @array) {
      print OUT $_;
      for(my $i=0; $i <= $#array; $i ++) {
        print OUT "\t$array[$i]";
      }
      print OUT "\n";
    }
    else {
      print OUT $_;
      for(my $i=0; $i <= $elm_num - 1; $i ++) {
        print OUT "\t-1";
      }
      print OUT "\n";
    }
  }
  close(IN);
  close(OUT);
}

sub format_trans {
  my $bismark_output = shift;
  my $cpg_sites = shift;
  my %total_counts;
  my %methy_counts;
  my $chrom_idx = -1;
  my $position_idx = -1;
  my %cpg;
  
  open(IN, $cpg_sites =~ /.gz$/ ? "gzip -d -c $cpg_sites |" : "<$cpg_sites") or die "$cpg_sites: $!\n";
  while(<IN>) {
    chomp;
    next if(/chrom/gi && /start/gi && /end/gi);
    
    my @array = split("\t");
    if(!exists $cpg{join("\t", $array[0], $array[1])}) {
      $cpg{join("\t", $array[0], $array[1])} = join("\t", $array[0], $array[1]);
      $cpg{join("\t", $array[0], $array[1]+1)} = join("\t", $array[0], $array[1]);
    }
  }
  close(IN);
  
  open(IN, $bismark_output =~ /.gz$/ ? "gzip -d -c $bismark_output |" : "<$bismark_output") or die "$bismark_output: $!\n";
  while(<IN>) {
    chomp;
    next if(/methylation/gi && /extractor/gi);
    
    if($pattern eq "null" || /$pattern/gi) {
      my @array = split("\t");
      
      if($chrom_idx < 0 || $position_idx < 0) {
        for(my $i = 0; $i < @array; $i ++) {
          if($array[$i] =~ /^chr/i) {$chrom_idx = $i;}
          if($array[$i] =~ /^\d+$/ && $array[$i] > 0) {$position_idx = $i;}
        }
      }
      
      next if !exists $cpg{"$array[$chrom_idx]\t$array[$position_idx]"};
      my $ID = $cpg{"$array[$chrom_idx]\t$array[$position_idx]"};
      
      
      if(!exists $total_counts{$ID}) {
        $total_counts{$ID} = 1;
      } else {
        $total_counts{$ID} ++;
      }
      
      if($array[4] eq "X" || $array[4] eq "H" || $array[4] eq "Z") {
        if(!exists $methy_counts{$ID}) {
          $methy_counts{$ID} = 1;
        } else {
          $methy_counts{$ID} ++;
        }
      }
    }
    
  }
  close(IN);

  $bismark_output =~ s/.gz$//;
  $bismark_output =~ s/.txt$//;
  
  if($abrr ne "none") {
    open(OUT, "| gzip -c > ${bismark_output}.${abrr}.txt.gz") or die "${bismark_output}.${abrr}.txt.gz: $!\n"; 
  } else {
    open(OUT, "| gzip -c > ${bismark_output}.meInfo.txt.gz") or die "${bismark_output}.meInfo.txt.gz: $!\n"; 
  }
  
  print OUT "chrom\tposition\tmethylation_level\ttotal_counts\tmethylated_counts\n";
  foreach my $key (sort{$a cmp $b} keys %total_counts) {
    if(exists $methy_counts{$key}) {
      my $ml = $methy_counts{$key} / $total_counts{$key};
      print OUT "$key\t$ml\t$total_counts{$key}\t$methy_counts{$key}\n";
    } else {
      print OUT "$key\t0\t$total_counts{$key}\t0\n";
    }
  }
  close(OUT);
}

sub binary_search {
  my $array = shift;
  my $start = shift;
  my $end = shift;
  my $tmpArray = shift;
  
  my $left = 0;
  my $right = $#$array;
  my $mid = int(($left + $right) / 2);
  
  if($end < $$array[$left][0] || $start > $$array[$right][0]) {
    @$tmpArray = ();
  }
  else {
    while($right - $left > 1) {
      if($end < $$array[$mid][0]) {
        $right = $mid;
        $mid = int(($left + $right) / 2);
      }
      elsif($start > $$array[$mid][0]) {
        $left = $mid;
        $mid = int(($left + $right) / 2);
      }
      else {
        last;
      }
    }
    
    if($start <= $$array[$mid][0] && $end >= $$array[$mid][0]) {
      
      my @count;
      if($is_avg) {@count = map {0} @{$$array[$mid]}; shift @count;}
      for(my $i = 1; $i <= $#{$$array[$mid]}; $i ++) {
        push(@$tmpArray, $$array[$mid][$i]);
        if($is_avg && $$array[$mid][$i] ne "NA") {$count[$i - 1] ++;}
      }
      
      $left = $mid - 1;
      $right = $mid + 1;
      while($left >= 0 && $start <= $$array[$left][0] && $end >= $$array[$left][0]) {
        for(my $i = 1; $i <= $#{$$array[$left]}; $i ++) {
          if($$array[$left][$i] ne "NA") {
            if($$tmpArray[$i-1] eq "NA") {$$tmpArray[$i-1] = $$array[$left][$i];}
            else {$$tmpArray[$i-1] += $$array[$left][$i];}
            
            if($is_avg) {$count[$i-1] ++;}
          }
        }
        $left --;
      }
      
      while($right <= $#$array && $start <= $$array[$right][0] && $end >= $$array[$right][0]) {
        for(my $i = 1; $i <= $#{$$array[$right]}; $i ++) {
          if($$array[$right][$i] ne "NA") {
            if($$tmpArray[$i-1] eq "NA") {$$tmpArray[$i-1] = $$array[$right][$i];}
            else {$$tmpArray[$i-1] += $$array[$right][$i];}
            
            if($is_avg) {$count[$i - 1] ++;}
          }
        }
        $right ++;
      }
      
      if($is_avg) {
        for(my $i = 0; $i <= $#$tmpArray; $i ++) {
          if($$tmpArray[$i] ne "NA") {$$tmpArray[$i] /= $count[$i];}
        }
      }
    }
    else {
      @$tmpArray = ();
    }
  }
}