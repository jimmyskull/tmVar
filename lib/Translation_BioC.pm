#===================================================
# Author:Chih-Hsuan Wei
# Issue: ME
# Description: Translate CRF result to PubTator Format
#===================================================

package ME;

BEGIN {
        push(@INC, 'lib/BioC');
}
use BioC_full;

use encoding 'big5', STDIN => 'big5', STDOUT => 'big5';

sub Translation_BioC
{
	my $PubTator_intput=$_[0];
	my $CRF_output=$_[1];
	my $location=$_[2];
	my $PubTator_output=$_[3];
	
	my %sentence_hash=();
	my %article_hash=();
	my %TypeOrder_hash=();
	my @pmidARR;
	my $BioC_text="";
	
	my $input_collection = new BioC_full::Collection();
	my $input_xml = new BioC_full::Connector_libxml();
	$input_xml->start_read($PubTator_intput, $input_collection);
	
	my $document = new BioC_full::Document();
	while ( $input_xml->read_next($document) ) # documents
	{
		my $pmid = $document->{id};
		for(my $i=0;$i<$document->{passages}->size();$i++) # passages
		{
			my $type = $document->{passages}->get($i)->{infons}->get("type"); #infon
			my $offset = $document->{passages}->get($i)->{offset}; #offset
			my $sentence = $document->{passages}->get($i)->{text}; #text
			$sentence=~s/[\n\r]//g;
			$sentence_hash{$pmid."\t".$type."_".$offset}=$sentence;
			$article_hash{$pmid}=$article_hash{$pmid}.$sentence." ";
			$TypeOrder_hash{$pmid}=$TypeOrder_hash{$pmid}.$type."_".$offset."\t";
			
			push(@pmidARR,$pmid);
		}
	}
	
	my %output_hash=();
	my $count=0;
	open output,"<".$CRF_output;
	while(<output>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$output_hash{$count}=$tmp;
		$count++;
	}
	
	my %location_hash=();
	$count=0;
	open location,"<".$location;
	while(<location>)
	{
		my $location=$_;
		$location=~s/[\n\r]//g;
		$location_hash{$count}=$location;
		$count++;
	}
	
	my %printSTR_hash=();
	for(my $i=0;$i<$count;$i++)
	{
		my $output=$output_hash{$i};
		my $location=$location_hash{$i};
		
		my $start=100000;
		my $last=0;
		my $pmid="";
		my $mention_tmp="";
		#reference seq	A
		#mutation type	T
		#Position	P
		#WildType	W
		#Mutant	M
		#FrameShift	F
		#FrameShiftPosition	S
		#DuplicationTime	D
		#SNP	R
		#Otherpartofmutation	I
		my %identifier_hash=();
		if($output=~/	([ATPWMFSDIR])$/ && $output!~/^;/)
		{
			my $prestate="";
			while($output=~/	([ATPWMFSDIR])$/ && $output!~/^;/)
			{
				$state=$1;
				if($location=~/^([^\t]+)	([^\t]+)	([0-9]+)	([0-9]+)/)
				{
					$pmid=$1;
					$mention=$2;
					$mention_tmp=$mention_tmp.$mention;
					if($state ne $prestate && $identifier_hash{$state} ne "")
					{
						$identifier_hash{$state}=$identifier_hash{$state}.",".$mention;
					}
					else
					{
						$identifier_hash{$state}=$identifier_hash{$state}.$mention;
					}
					if($3<$start){$start=$3;}
					if($4>$last){$last=$4;}
				}
				$i++;
				$output=$output_hash{$i};
				$location=$location_hash{$i};
				$prestate=$state;
			}		

			#identifier
			my $identifier="";
			my $Type="";
			if($identifier_hash{"D"} ne "") #dup
			{
				$identifier=$identifier_hash{"A"}."|".$identifier_hash{"T"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"}."|".$identifier_hash{"D"};
			}
			elsif($identifier_hash{"T"} ne "")
			{
				$identifier=$identifier_hash{"A"}."|".$identifier_hash{"T"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"};
			}
			elsif($identifier_hash{"S"} ne "")
			{
				$identifier=$identifier_hash{"A"}."|".$identifier_hash{"W"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"}."|".$identifier_hash{"F"}."|".$identifier_hash{"S"};
				$Type="ProteinMutation";
			}
			elsif($identifier_hash{"F"} ne "")
			{
				$identifier=$identifier_hash{"A"}."|".$identifier_hash{"W"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"}."|".$identifier_hash{"F"};
				$Type="ProteinMutation";
			}
			elsif($identifier_hash{"R"} ne "")
			{
				$identifier=$mention_tmp;
				$Type="SNP";
			}
			else
			{
				$identifier=$identifier_hash{"A"}."|".$identifier_hash{"W"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"};
			}
			
			#Type (ProteinMutation|DNAMutation|SNP)
			if($Type eq "")
			{
				if($identifier_hash{"T"}=~/[Dd]elta/)
				{
					$Type="DNAMutation";
				}
				elsif($identifier_hash{"P"}=~/^([Ee]x|EX|[In]ntron|IVS|[Ii]vs)/)
				{
					$Type="DNAMutation";
				}
				elsif($identifier_hash{"M"}!~/^[ATCGatcgu]+$/ || ($identifier_hash{"W"}!~/^[ATCGatcgu]+$/))
				{
					$Type="ProteinMutation";
				}
				else #default
				{
					$Type="DNAMutation";
				}
				if($identifier_hash{"A"}=~/(c|r|m|g|C)/)
				{
					$Type="DNAMutation";
				}
				elsif($identifier_hash{"A"}=~/p/)
				{
					$Type="ProteinMutation";
				}
			}
			
			if( ( length($identifier_hash{"W"})==3 || length($identifier_hash{"M"})==3 ) 
				&& length($identifier_hash{"W"}) != length($identifier_hash{"M"})
				&& $identifier_hash{"W"} ne "" && $identifier_hash{"M"} ne "" && $identifier_hash{"W"}!~/,/ && $identifier_hash{"M"}!~/,/
				&& ($identifier_hash{"W"}!~/^[ATCG]+$/ || $identifier_hash{"M"}!~/^[ATCG]+$/
				&& $identifier_hash{"T"} eq "")
				){} #remove Wildtype & Mutant are not the same length
			elsif(  $identifier_hash{"M"} eq "" && $identifier_hash{"T"} eq "" && $identifier_hash{"F"} eq "" && $identifier_hash{"P"} ne ""){} #Arg235
			#elsif(	$identifier_hash{"P"}=~/,/ && ( $identifier_hash{"W"}!~/,/ || $identifier_hash{"M"}!~/,/ )	){} #remove mutation with more than one positions
			elsif(	$identifier_hash{"T"}=~/Delta/ && $identifier_hash{"M"}=~/[A-Z]/ ){} #DeltaG
			#elsif(	$identifier_hash{"P"}=~/,/){} #Multiple positions
			elsif(	$identifier_hash{"P"}=~/^-/ && $Type eq "ProteinMutation"){} #negative protein mutation
			elsif(	$identifier_hash{"W"}=~/^[BJOUZ]/ || $identifier_hash{"M"}=~/^[BJOUZ]/){} #not a mutation
			elsif($identifier_hash{"W"} eq "" || $identifier_hash{"W"} ne $identifier_hash{"M"})
			{
				$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".($start-1)."	".$last."	".substr($article_hash{$pmid},$start-1,$last-($start-1))."	$Type	".$identifier."\n";
			}
			if($output!~/	O$/)
			{
				$i--;
			}
			$mention_tmp="";
		}		
	}
	close output;
	close location;
	
	open output,">".$PubTator_output;
	my %printed_hash=();
	foreach my $pmid (@pmidARR)
	{
		if(not exists $printed_hash{$pmid})
		{
			my @type=split("\t",$TypeOrder_hash{$pmid});
			foreach my $type(@type)
			{
				print output $pmid."|".$type."|".$sentence_hash{$pmid."\t".$type}."\n";
			}
			print output $printSTR_hash{$pmid}."\n";
			$printed_hash{$pmid}=1;
		}
	}
	close output;
	return 1;
}

return 1;
