# Annotation-of-genomic-regions
The Perl script may annotate any BED-like data to a given BED-like data with a faster way. 

Customized Perl script: AnnotationExe.pl. 

Usage:

perl AnnotationExe.pl [-h/--help] [-m/--mode <mode>] [-f/--flag] [-a/--cal-average] <-d/--dir dir> <-t/--table table>
where the argument [-m/--mode] is with optional value: 'sensitive' or 'stable' (default)
where the argument [-f/--flag] is a Boolean variable (default: False)
where the argument [-a/--cal-average] is a Boolean variable (default: False), which depends on [-f/--flag]



Note: the argument [-m/--mode] is with optional value: 'sensitive' or 'stable' (default). Here ‘stable’ means that the overlapping section between two loci has at least a half of the target locus; ‘sensitive’ means that the overlapping section between two loci is at least 1bp. The tool supports to input several separate tables (each table includes 3 (such as ‘chrom start end’) or 4 columns and the last column represents annotation information) which are included in a folder specified by <-d/--dir>. The table to be annotated is required to include at least 3 columns (the first 3 columns: ‘chrom start end’), which is specified by <-t/--table>. If one uses [-f/--flag] argument, it means that the tool will print details of annotation(s); Or it will print the count of overlapping. If the fourth column of annotation(s) is numerical, one may use [-a/--cal-average] argument to calculate the average of value across all loci overlapped for a given target locus.

For example : perl AnnotationExe.pl –d /myData/tables/ -t myTable.txt 
