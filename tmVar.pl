#!/usr/bin/env perl

#===================================================
# Author: Chih-Hsuan Wei
# Software: tmVar
# Last update date: 2013/04/25
#===================================================

#=====
#Mutation Extraction
sub ME
{
	my $input=$_[0];
	my $filename=$_[1];
	my $output_route=$_[2];
	my $RegEx=$_[3];
	my $Mutation_Extraction=$_[4];
	my $setup=$_[5];
	my $format=$_[6];
	if($Mutation_Extraction eq "True")
	{
		if(-e $input)
		{

			if($format eq "BioC")
			{
				require './lib/ExtractFeature_BioC.pm';
				ME::ExtractFeature_BioC($input,"tmp/".$filename.".data","tmp/location_".$filename.".data");
			}
			else
			{
				require './lib/ExtractFeature.pm';
				ME::ExtractFeature_PubTator($input,"tmp/".$filename.".data","tmp/location_".$filename.".data");
			}
			close input;

			if($^O=~/Win/)
			{
				open(CRF,"|crf_test -m CRF/MentionExtractionUB.Model -o tmp/".$filename.".output tmp/".$filename.".data")|| die ("\nAn Error in CRF module. \nPlease re-install CRF module:\n\n\tsh Installation.sh\n\n");	close(CRF);
				close CRF;
			}
			else
			{
				open(CRF,"|crf_test -m CRF/MentionExtractionUB.Model -o tmp/".$filename.".output tmp/".$filename.".data")|| die ("\nAn Error in CRF module. \nPlease re-install CRF module:\n\n\tsh Installation.sh\n\n");	close(CRF);
				close CRF;
			}

			if($format eq "BioC")
			{
				require './lib/Translation_BioC.pm';
				ME::Translation_BioC($input,"tmp/".$filename.".output","tmp/location_".$filename.".data","tmp/".$filename.".ME");
			}
			else
			{
				require './lib/Translation.pm';
				ME::Translation_PubTator($input,"tmp/".$filename.".output","tmp/location_".$filename.".data","tmp/".$filename.".ME");
			}

			require './lib/PostProcessing.pm';
			ME::PostProcessing("tmp/".$filename.".ME","tmp/".$filename.".PostME",$RegEx);
			ME::ExtractStateFeature("tmp/".$filename.".PostME","tmp/".$filename.".PostME.data",$filename);

			if($format eq "BioC")
			{
				require './lib/OutputFormat_BioC.pm';
				ME::OutputFormat_BioC($input,"tmp/".$filename.".PostME","tmp/".$filename.".PostME.output",$output_route."/".$filename.".BioC.xml",$setup);
			}
			else
			{
				require './lib/OutputFormat.pm';
				ME::OutputFormat_PubTator("tmp/".$filename.".PostME","tmp/".$filename.".PostME.output",$output_route."/".$filename.".PubTator",$setup);
			}

			open clear,"|rm tmp/*".$filename."*";	close clear;
			#open clear,"|rm input/".$filename;	close clear;

		}
		else
		{
			die "\nCannot find input data.\n";
		}
	}
}

sub main
{
	my $folder_route;
	my $output_route;
	my $setup;
	for(my $i=0;$i<@ARGV;$i++)
	{
		if($ARGV[$i] eq "-i")
		{
			$i++;
			$folder_route=$ARGV[$i];
		}
		elsif($ARGV[$i] eq "-o")
		{
			$i++;
			$output_route=$ARGV[$i];
		}
		elsif($ARGV[$i] eq "-s")
		{
			$i++;
			$setup=$ARGV[$i];
		}
		elsif($ARGV[$i]=~/^-i(.+)$/)
		{
			$folder_route=$1;
		}
		elsif($ARGV[$i]=~/^-o(.+)$/)
		{
			$output_route=$1;
		}
		elsif($ARGV[$i]=~/^-s(.+)$/)
		{
			$setup=$1;
		}
	}
	my %setup_hash=();

	if($folder_route eq "")
	{
		print "Instruction Format:\n\n\tperl tmVar.pl -i [input dir] -o [output dir] -s [setup]\n";
		print "\te.g. perl tmVar.pl -i input -o output\n";
	}
	elsif($output_route eq "")
	{
		print "Instruction Format:\n\n\tperl tmVar.pl -i [input dir] -o [output dir] -s [setup]\n";
		print "\te.g. perl tmVar.pl -i input -o output\n";
	}
	else
	{
		if($setup eq "")
		{
			$setup="setup.txt";
			print "-s [setup]: setup.txt\n";
		}
		$setup_hash{"Mutation_Extraction"} = "True";
		my $RegEx="RegEx/RegEx.NL";

		#=====
		#tmVar processing
		opendir(DIR, $folder_route);
		@class = grep(/[a-z0-9]/,readdir(DIR));
		closedir(DIR);
		foreach $filename(@class)
		{
			my $countX=0;
			my $format="";
			my $checkformat=0;
			open input,"<".$folder_route."/".$filename;
			while(<input>)
			{
				my $tmp=$_;
				if($tmp=~/<infon key=\"type\">/){$format="BioC"; $checkformat=1; last;}
				elsif($tmp=~/^(.+)\|(.+)\|(.*)$/){$format="PubTator"; $checkformat=1; last;}
			}
			close input;
			if($checkformat==0){print "Format Error. It should follow PubTator/BioC format.\n";exit;}
			elsif($format eq "BioC") {print "Input format: BioC\n";}
			elsif($format eq "PubTator")
			{
				my $pmid="";
				my $cut=0;
				open input,"<".$folder_route."/".$filename;
				while(<input>)
				{
					my $tmp=$_;
					if($tmp=~/^([^\t]+)\|([^\t]+)\|(.*)$/)
					{
						if($pmid eq $1 && $cut==1)
						{
							print "\nPubTator Format Error: Blank lines are necessary for different pmids but no extra lines between the same pmids.\n\n";
							print "An example is shown in below:\n\n";
							print "21738389|t|A novel DFNB31 mutation associated with Usher ...\n21738389|a|PURPOSE: To identify the genetic defect of a consanguineous ...\n\n22051099|t|Variation in the CXCR1 gene (IL8RA) is not associated with ...\n22051099|a|BACKGROUND: The chemokine receptor 1 CXCR-1 (or IL8R-alpha) ...\n";
							exit;
						}
						else
						{
							$pmid=$1;
							$cut=0;
						}
					}
					else
					{
						$cut=1;
					}
				}
				close input;

				print "Input format: PubTator\n";
			}

			open input,"<".$folder_route."/".$filename;
			while(<input>)
			{
				my $tmp=$_;
				if($tmp=~/<infon key=\"type\">/){$format="BioC";}
				if($format eq "BioC")
				{
					my @cont=($tmp=~/<collection>/g);
					$cont=@cont;
					$countX=$countX+$cont;
				}
				else
				{
					if($tmp=~/^[\r\n]+$/){$countX++;}
				}
			}
			close input;

			my ($sec1,$min1,$hour1,$day1,$mon,$year)=localtime(time);

			#=====
			#Mutation Extraction
			if($^O=~/Win/ && $format eq "BioC")
			{
				print "Running tmVar on $countX docs in $filename ... Failed: BioC format is not workable for ".$^O.".                    \n";
			}
			else
			{
				ME($folder_route."/".$filename,$filename,$output_route,$RegEx,$setup_hash{"Mutation_Extraction"},$setup,$format);

				my ($sec2,$min2,$hour2,$day2,$mon,$year)=localtime(time);
				$hour1=$hour1+$day1*24;
				$hour2=$hour2+$day2*24;
				my $timecost=((($hour2-$hour1)*60)+($min2-$min1)*60)+($sec2-$sec1);
				print "Running tmVar on $countX docs in $filename ... Finished in $timecost seconds. \n";
			}
		}
	}
}

main();
