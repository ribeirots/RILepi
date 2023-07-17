#! /usr/bin/env perl
use strict;


#my $file = shift;

#my $split_num = "4000000";
my $split_num = $ARGV[0];
my $jobfolder = $ARGV[1];

my $reads = "";
chdir($jobfolder);
#opendir DS, "$jobfolder" or die "Can not open dataset<.>\n";
opendir DS, "./" or die "Can not open dataset<.>\n";
my $split = "no";
foreach my $subfile (readdir DS)
{
	next if $subfile =~ /^\.\.?$/;
	if($subfile =~ /(.*?)\.fastq/) {
		$reads = $1;
 		print "reading $jobfolder/$subfile\n";
		#my $total_lines_cmd = "cat $jobfolder/$subfile | wc -l";
		my $total_lines_cmd = "cat $subfile | wc -l";
		print "cmd=$total_lines_cmd\n";
		my $total_lines = `$total_lines_cmd`;
		print "total number of lines $total_lines\n";
		$total_lines = $total_lines / 4;
		print "total reads in file $subfile = $total_lines\n";
		if ($total_lines == 0) {
			print "this file $subfile is empty";
			exit;
		}
		#emtpy file exit now
		if ($total_lines == 0) {
			exit;
		}
		elsif ($total_lines > $split_num) {
			$split = "yes";
			my $splitname = "$reads" . "_Block";
			#my $splitcmd = "cat $jobfolder/$subfile |split -a 5 -d -l $split_num - $splitname";
			my $splitcmd = "cat $subfile |split -a 5 -d -l $split_num - $splitname";       
			print "running $splitcmd\n";
			system($splitcmd);
			print "+++++++++++++++++++++++++++++++++++++++\n";
		}
		else {
			print "remaning file because we dont need to split!\n";
			my $file = "$subfile " . "$reads" . "_Block00000.fastq";
			#my $mv_cmd = "cp $jobfolder/$file";
			my $mv_cmd = "cp $file";
			system($mv_cmd);
			#chdir("$jobfolder");
			#system("gzip $subfile/$file");
			system("gzip $file");
		}
	}
}
closedir(DS);


if ($split eq "yes") {
my @chunks = `ls *Block*`;
foreach my $name (@chunks) {
	chomp($name);
	print "<$name>\n";
	system("mv $name $name.fastq");
	system("gzip $name.fastq");
}
}
exit(0);
