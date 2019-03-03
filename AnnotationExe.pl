#!/usr/perl/bin

use strict;
use Getopt::Long;

my $help;
my $table_dir;
my $obj_table;

my @tables;
my %obj_list;
my %obj_list_ref;
my %com_list;
my %list;
my %tList;
my @tmpList;
my $chrom;
my $tab;
my $current_tab;
my ($i, $j, $pre_j);
my @obj_value;
my $tmp;
my @tmpArray;
my $flag;
my $tmp_ind;
my $mode;
my $com_ncol;
my $avg = 0;
my $pattern = "\\S+";
my $notBedCoordi;

GetOptions(
  'h|help' => \$help,
  'd|dir=s' => \$table_dir,
  't|table=s' => \$obj_table,
  'm|mode=s' => \$mode,
  'f|flag' => \$flag,
  'a|cal-average' => \$avg,
  'p|pattern=s' => \$pattern,
  'n|notBedCoordi' => \$notBedCoordi, 
);

sub print_help {
  print "perl AnnotationExe.pl [-h/--help] [-m/--mode] [-f/--flag] [-a/--cal-average] <-d/--dir> <dir> <-t/--table> <table>\n";
  print "where the argument [-m/--mode] is with optional value: 'sensitive' or 'stable' (default)\n";
  print "where the argument [-f/--flag] is a Boolean variable (default: False)\n";
  print "where the argument [-a/--cal-average] is a Boolean variable (default: False), which depends on [-f/--flag]\n";
  exit;
}

warn "Must specify valid directory which stores tables annoated CSM loci. \n" unless -d $table_dir;
warn "Must specify valid file path which represents the table to cross. \n" unless -f $obj_table;

if($help || !-d $table_dir || !-f $obj_table) {
  print_help();
}

unless($mode){
  $mode = "stable"; # optional value: 'sensitive' or 'stable' (default)
}

unless($flag){
  $flag = "0"; # optional value: '0' (default) or '1'
}

opendir(DIR, $table_dir) or die "Cannot open directory $table_dir \n";
while(my $file = readdir(DIR)){
  push @tables, $file if -f "$table_dir/$file" && $file =~ /$pattern/;
}
closedir(DIR);

open(IN, $obj_table =~ /gz$/ ? "gunzip -c $obj_table |" : "<$obj_table") or die "Cannot open table $obj_table. \n";
my $header = <IN>; # skim title 
chomp $header;

%obj_list_ref = ();
while(<IN>){
  chomp;
  s/[\r|\n]//g;
  
  @tmpArray = split("\t");
  if(!exists $obj_list_ref{$_}) {
    $obj_list_ref{$_} = 1;
    push(@{$obj_list{$tmpArray[0]}}, [(@tmpArray)[1..$#tmpArray]]);
  }
}
close(IN);
undef %obj_list_ref;

foreach $chrom (keys %obj_list){
  @{$obj_list{$chrom}} = sort {$a -> [0] <=> $b -> [0]} @{$obj_list{$chrom}};
}

@tables = sort{$a cmp $b} @tables;
for(my $ind = 0; $ind < @tables; $ind ++){
  $tab = $tables[$ind];
  $tab =~ /^(\S+)\.\S+$/; 
  $tables[$ind] = $1;
  $com_ncol = 4;
  
  %tList = (); %com_list = ();
  $current_tab = "$table_dir/$tab";
  open(IN, $current_tab =~ /gz$/ ? "gunzip -c $current_tab |" : "<$current_tab") or die "Cannot open table $current_tab. \n";
  <IN>; # skim title
  while(<IN>){
    chomp;
    s/[\r|\n]//g;
    
    @tmpArray = split("\t");
    if(!exists $tList{$_}) {
      $tList{$_} = $_;
      if(!$notBedCoordi) {push(@{$com_list{$tmpArray[0]}}, [(@tmpArray)[1..$#tmpArray]]);}
      else {push(@{$com_list{$tmpArray[0]}}, [$tmpArray[1], $tmpArray[1]+1, (@tmpArray)[2..$#tmpArray]]);}
    }
    
    if($com_ncol > @tmpArray) {$com_ncol = scalar @tmpArray;}
  }
  close(IN);
  
  foreach $chrom (keys %com_list){
    @{$com_list{$chrom}} = sort {$a -> [0] <=> $b -> [0]} @{$com_list{$chrom}};
  }
    
  foreach $chrom (keys %obj_list){
  
    $pre_j = 0; $j = 0; 
    for($i = 0; $i < @{$obj_list{$chrom}}; $i++){
    
      if($mode =~ /sensitive/i){
        @obj_value = @{${$obj_list{$chrom}}[$i]};
      }
      else{
        $tmp = (${$obj_list{$chrom}}[$i][0] + ${$obj_list{$chrom}}[$i][1]) / 2;
        @obj_value = ($tmp, $tmp);
      }
      
      @tmpList = ();
      $tmp = join("\t", $chrom, @{${$obj_list{$chrom}}[$i]});
      
      while($obj_value[0] > ${$com_list{$chrom}}[$j][1]){
        $j ++;
        
        if($j >= @{$com_list{$chrom}}) {

          if($flag){
          
            if(exists $list{$tmp}){
              $list{$tmp} .= "\tnull";
            }
            else{
              $list{$tmp} = "null";
            }
          }
          else{
          
            if(exists $list{$tmp}){
              $list{$tmp} .= "\t0";
            }
            else{
              $list{$tmp} = "0";
            }
          }
          
          last;
        }
      }
      
      if(! exists $com_list{$chrom} || $obj_value[1] < ${$com_list{$chrom}}[$j][0] || $j >= @{$com_list{$chrom}}){

        if($flag){
        
          if(exists $list{$tmp}){
            $list{$tmp} .= "\tnull";
          }
          else{
            $list{$tmp} = "null";
          }
        }
        else{
        
          if(exists $list{$tmp}){
            $list{$tmp} .= "\t0";
          }
          else{
            $list{$tmp} = "0";
          }
        }
        next;
      }
      
      if(!($obj_value[1] < ${$com_list{$chrom}}[$j][0] || $obj_value[0] > ${$com_list{$chrom}}[$j][1])) {
        $pre_j = $j;
      }
      
      while(!($obj_value[1] < ${$com_list{$chrom}}[$j][0] || $obj_value[0] > ${$com_list{$chrom}}[$j][1])){
      
        if($flag){
          if($com_ncol<4) {
            push(@tmpList, "$chrom:${$com_list{$chrom}}[$j][0]-${$com_list{$chrom}}[$j][1]");
          } else {
            push(@tmpList, join(",", (@{${$com_list{$chrom}}[$j]})[2..$#{${$com_list{$chrom}}[$j]}])); 
            #push(@tmpList, "$chrom:${$com_list{$chrom}}[$j][0]-${$com_list{$chrom}}[$j][1]"); 
          }
        }
        else{
          if(!$avg) {
            if($com_ncol<4) {
              if(@tmpList <= 0) {push(@tmpList, 1);} else {$tmpList[0] ++;}
            } else {
              if(@tmpList <= 0) {
                push(@tmpList, ${$com_list{$chrom}}[$j][2]);
              }
              else {
                $tmpList[0] += ${$com_list{$chrom}}[$j][2];
              }
            }
          } else {
            if(@tmpList <= 0) {
              push(@tmpList, 1); 
              push(@tmpList, ${$com_list{$chrom}}[$j][2]);
            }
            else {
              $tmpList[0] ++;
              $tmpList[1] += ${$com_list{$chrom}}[$j][2];
            }
          }
        }
      
        $j ++;
        last if($j >= @{$com_list{$chrom}});
      }
      
      if(@tmpList > 0){
        
        if($flag){
          %tList = ();
          foreach my $el (@tmpList) {
            $tList{$el} = $el if !exists $tList{$el};
          }
          
          if(exists $list{$tmp}){
          
            $list{$tmp} .= "\t".join(";", keys %tList);
          }
          else{
            $list{$tmp} = join(";", keys %tList);
          }
        }
        else {
          if(exists $list{$tmp}){
          
            if(!$avg) {$list{$tmp} .= "\t".$tmpList[0];}
            else {$list{$tmp} .= "\t".($tmpList[1]/$tmpList[0]);}
          }
          else{
            if(!$avg) {$list{$tmp} = $tmpList[0];}
            else {$list{$tmp} = ($tmpList[1]/$tmpList[0]);}
          }
        }
        
        $tmp_ind = $pre_j;
	      $pre_j = $j;
	      $j = $tmp_ind;
	    }
	  }
  }
  
  undef %tList;
  undef %com_list;
  undef @obj_value;
}

$obj_table =~ s/\.txt$//;
if($pattern ne "\\S+") {
  open(OUT, ">${obj_table}.${pattern}.annotated.txt") or die "Cannot create file ${obj_table}.${pattern}.annotated.txt. \n";
} else {
  open(OUT, ">${obj_table}.annotated.txt") or die "Cannot create file ${obj_table}.annotated.txt. \n";
}

print OUT "$header\t", join("\t", @tables), "\n";
foreach my $key (sort{$a cmp $b} keys %list){
  my @array = split("\t", $list{$key});
  if(@array != @tables) {
    
    if($flag) {
      my $tmp = join"", map {"\tnull"} @tables;
      print OUT "$key$tmp\n";
    }
    else {
      my $tmp = join"", map {"\t0"} @tables;
      print OUT "$key$tmp\n";
    }
  }
  else {
    print OUT join("\t", $key, (@array)[0..$#tables]), "\n";
  }
}
close(OUT);
exit(0);