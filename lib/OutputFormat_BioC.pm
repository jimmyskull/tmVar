#===================================================
# Author:Chih-Hsuan Wei
# Issue: ME
# Description: Post-Processing of Mutation Extraction
#===================================================

package ME;

BEGIN {
        push(@INC, 'lib/BioC');
}
use BioC_full;

use encoding 'big5', STDIN => 'big5', STDOUT => 'big5';

sub OutputFormat_BioC
{
	my $BioC_input=$_[0];
	my $PubTator_input=$_[1];
	my $CRF_output=$_[2];
	my $BioC_output=$_[3];
	my $setup=$_[4];
	
	my %setup_hash=();
	if($setup eq ""){$setup="setup.txt";}
	open input,"<".$setup;
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\W]+)[ ]*=[ ]*([^\W]+)/)
		{
			$setup_hash{$1}=$2;
		}
	}
	close input;

	
	my %output_hash=();
	my $count=0;
	open input,"<".$CRF_output;
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$output_hash{$count}=$tmp;
		$count++;
	}
	close input;
	
	my %STR_hash=();
	my $str_count=1;
	for(my $i=0;$i<$count;$i++)
	{
		my $output=$output_hash{$i};
		my $mention_tmp="";
		my %identifier_hash=();
		my $prestate="";
		
		my $keepi=$i;
		my $MW_right=0;
		my $MW_left=10000;
		my $OtherType=0;
		while($output=~/^([^\t]+).*([ATPWMFSDIR])$/)
		{
			my $mention=uc($1);
			if($output=~/T$/)
			{
				$OtherType=1;
			}
			if($output=~/^([^\t]+).*([MW])$/)
			{
				if($MW_left>length($mention_tmp)){$MW_left=length($mention_tmp);}
				if($MW_right<(length($mention_tmp)+length($mention))){$MW_right=length($mention_tmp)+length($mention);}
			}
			$mention_tmp=$mention_tmp.$mention;
			$keepi++;
			$output=$output_hash{$keepi};
		}
		
		my $boundary=0;
		my $WMstate="F";	#Forward
		if($mention_tmp=~/^(.+?)(FOR|INPLACEOF|INSTEADOF|MUTANTOF|MUTANTSOF|RATHERTHAN|RATHERTHAN|REPLACEMENTOF|REPLACEMENTSOF|REPLACES|REPLACING|RESIDUEATPOSITION|RESIDUEFOR|RESIDUEINPLACEOF|RESIDUEINSTEADOF|SUBSTITUTIONAT|SUBSTITUTIONFOR|SUBSTITUTIONOF|SUBSTITUTIONSAT|SUBSTITUTIONSFOR|SUBSTITUTIONSOF|SUBSTITUTEDFOR|TOREPLACE)(.+)$/
			&& $MW_left<length($1) && length($1) < $MW_right )
		{
			$boundary=length($1);
			$WMstate="B";	#Backward
		}
		elsif($mention_tmp=~/^(.+?)(>|TO|INTO|OF|BY|WITH|AT)(.+)$/)
		{
			$boundary=length($1);
		}
		$mention_tmp="";
		my $start=0;
		
		$output=$output_hash{$i};
		while($output=~/^([^\t]+).*([ATPWMFSDIR])$/)
		{
			$mention=uc($1);
			$state=$2;
			$mention_tmp=$mention_tmp.$mention;
			if($state=~/[WM]/ && $boundary != 0)
			{
				if($WMstate eq "F" && $OtherType==0)
				{
					if($start<=$boundary){$state="W";}
					else{$state="M";}
				}
				elsif($WMstate eq "B" && $OtherType==0)
				{
					if($start<=$boundary){$state="M";}
					else{$state="W";}
				}				
				if($state ne $prestate && $identifier_hash{$state} ne "")
				{
					$identifier_hash{$state}=$identifier_hash{$state}.",".$mention;
				}
				else
				{
					$identifier_hash{$state}=$identifier_hash{$state}.$mention;
				}
			}
			else
			{
				if($state ne $prestate && $identifier_hash{$state} ne "")
				{
					$identifier_hash{$state}=$identifier_hash{$state}.",".$mention;
				}
				else
				{
					$identifier_hash{$state}=$identifier_hash{$state}.$mention;
				}
			}
			$i++;
			$output=$output_hash{$i};
			$prestate=$state;
			$start=$start+length($mention);
		}		
		
		my $translate="";
		my $NotAaminoacid="";
		
		my %nametothree = (GUANINE => "G", ADENINE => "A", CYTOSINE => "C", THYMINE => "T", ALANINE => "ALA", ALANINE => "ALA", ALANINE => "ALA", ALANINE => "ALA", ALANINE => "ALA", ARGININE => "ARG", ASPARAGINE => "ASN", "ASPARTICACID" => "ASP", "ASPARTATE" => "ASP", CYSTEINE => "CYS", GLUTAMINE => "GLN", "GLUTAMICACID" => "GLU", GLYCINE => "GLY", HISTIDINE => "HIS", ISOLEUCINE => "ILE", LEUCINE => "LEU", LYSINE => "LYS", METHIONINE => "MET", PHENYLALANINE => "PHE", PROLINE => "PRO", SERINE => "SER", THREONINE => "THR", TRYPTOPHAN => "TRP", TYROSINE => "TYR", VALINE => "VAL", STOP => "XAA");
		my %threetone = (ALA => "A", ARG => "R", ASN => "N", ASP => "D", CYS => "C", GLN => "Q", GLU => "E", GLY => "G", HIS => "H", ILE => "I", LEU => "L", LYS => "K", MET => "M", PHE => "F", PRO => "P", SER => "S", THR => "T", TRP => "W", TYR => "Y", VAL => "V", ASX => "B", GLX => "Z", XAA => "X",TER => "X");
		
		my @tmp=split(",",$identifier_hash{"M"});
		$identifier_hash{"M"}="";
		foreach my $tmp (@tmp)
		{
			if(exists $nametothree{$tmp})	
			{
				$tmp=$nametothree{$tmp};
				$translate="Y";	
			}
			if($identifier_hash{"M"} eq ""){$identifier_hash{"M"}=$tmp;}
			else{$identifier_hash{"M"}=$identifier_hash{"M"}.",".$tmp;}
		}
		my @tmp=split(",",$identifier_hash{"M"});
		$identifier_hash{"M"}="";
		foreach my $tmp (@tmp)
		{
			if(exists $threetone{$tmp})	
			{	
				$tmp=$threetone{$tmp};
				$translate="Y";	
			}
			elsif(length($tmp)>1)
			{
				$NotAaminoacid="N";
			}
			if($identifier_hash{"M"} eq ""){$identifier_hash{"M"}=$tmp;}
			else{$identifier_hash{"M"}=$identifier_hash{"M"}.",".$tmp;}
		}
			
		my @tmp=split(",",$identifier_hash{"W"});
		$identifier_hash{"W"}="";
		foreach my $tmp (@tmp)
		{
			if(exists $nametothree{$tmp})	
			{
				$tmp=$nametothree{$tmp};
				$translate="Y";	
			}
			if($identifier_hash{"W"} eq ""){$identifier_hash{"W"}=$tmp;}
			else{$identifier_hash{"W"}=$identifier_hash{"W"}.",".$tmp;}
		}
		my @tmp=split(",",$identifier_hash{"W"});
		$identifier_hash{"W"}="";
		foreach my $tmp (@tmp)
		{
			if(exists $threetone{$tmp})	
			{	
				$tmp=$threetone{$tmp};
				$translate="Y";	
			}
			elsif(length($tmp)>1)
			{
				$NotAaminoacid="N";
			}
			if($identifier_hash{"W"} eq ""){$identifier_hash{"W"}=$tmp;}
			else{$identifier_hash{"W"}=$identifier_hash{"W"}.",".$tmp;}
		}
		
		
		if($identifier_hash{"A"} eq "CDNA"){	$identifier_hash{"A"}="c";	}
		elsif($identifier_hash{"P"}=~/^[\+\-]/ && $identifier_hash{"A"} eq ""){	$identifier_hash{"A"}="c";	}
		$identifier_hash{"A"}=lc($identifier_hash{"A"});
		$identifier_hash{"F"}=~s/\*/X/g;
		$identifier_hash{"R"}=lc($identifier_hash{"R"});
		$identifier_hash{"R"}=~s/[\[\]]//g;
		
		if($identifier_hash{"P"}=~/^([0-9]+)\-([0-9]+)$/)
		{
			if($2>$1)
			{
				$identifier_hash{"P"}=$1."_".$2;
			}
		}
		if($identifier_hash{"T"} eq "DELTA") {$identifier_hash{"T"}="DEL";	$translate="N";}
		elsif($identifier_hash{"T"} eq "DELINS") {$identifier_hash{"T"}="INDEL";	$translate="N";}
		elsif($identifier_hash{"T"}=~/INS.*DEL/) {$identifier_hash{"T"}="INDEL";	$translate="N";}
		elsif($identifier_hash{"T"}=~/DEL.*INS/) {$identifier_hash{"T"}="INDEL";	$translate="N";}
		
		if($identifier_hash{"W"}=~/(.+),(.+)/ && ($identifier_hash{"M"} eq "" && $identifier_hash{"T"} eq ""))
		{
			$identifier_hash{"W"}=$1;
			$identifier_hash{"M"}=$2;
		}
		elsif($identifier_hash{"M"}=~/(.+),(.+)/ && ($identifier_hash{"W"} eq "" && $identifier_hash{"T"} eq ""))
		{
			$identifier_hash{"W"}=$1;
			$identifier_hash{"M"}=$2;
		}
		
		$identifier_hash{"P"}=~s/[\-\[]+$//g;
		$identifier_hash{"P"}=~s/^(POSITION|NUCLEOTIDE|\:)//g;
		if($Type eq "ProteinMutation"){$identifier_hash{"P"}=~s/^CODON//g;}
		
		if($identifier_hash{"M"}=~/(DEL|INS|DUP)/)
		{
			$identifier_hash{"T"}=$identifier_hash{"M"};
			$identifier_hash{"M"}="";
		}
		elsif($identifier_hash{"W"}=~/(DEL|INS|DUP)/)
		{
			$identifier_hash{"T"}=$identifier_hash{"W"};
			$identifier_hash{"W"}="";
		}
		
		#identifier
		my $Type="";
		my $identifier="";
		if($identifier_hash{""} eq "DUP") #dup
		{
			$identifier=$identifier_hash{"A"}."|".$identifier_hash{"T"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"}."|".$identifier_hash{"D"};
		}
		elsif($identifier_hash{"T"} ne "") #DEL/INS/INDEL
		{
			$identifier=$identifier_hash{"A"}."|".$identifier_hash{"T"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"};
		}
		elsif($identifier_hash{"F"} ne "") #Frameshit
		{
			$identifier=$identifier_hash{"A"}."|FS|".$identifier_hash{"W"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"}."|".$identifier_hash{"S"};
			$Type="ProteinMutation";
		}
		elsif($identifier_hash{"R"} ne "")
		{
			$identifier=$identifier_hash{"R"};
			$Type="SNP";
		}
		elsif($mention_tmp=~/^I([RSrs][Ss][0-9].+)$/)
		{
			$identifier=$1;
			$identifier=~s/[\W\-\_]//g;
			$identifier=lc($identifier);
			$Type="SNP";
		}
		else
		{
			$identifier=$identifier_hash{"A"}."|SUB|".$identifier_hash{"W"}."|".$identifier_hash{"P"}."|".$identifier_hash{"M"};
		}
		
		#Type (ProteinMutation|DNAMutation|SNP)
		if($Type eq "")
		{
			if($translate eq "Y")
			{
				$Type="ProteinMutation";
			}
			elsif($translate eq "N")
			{
				$Type="DNAMutation";
			}
			if($identifier_hash{"P"}=~/^([Ee]x|EX|[In]ntron|IVS|[Ii]vs)/)
			{
				$Type="DNAMutation";
				$identifier=~s/^[^\|]*\|/c\|/g;
			}
			elsif((($identifier_hash{"M"}=~/[ISQMNPKDFHLRWVEYX]/) || ($identifier_hash{"W"}=~/[ISQMNPKDFHLRWVEYX]/)))
			{
				$Type="ProteinMutation";
				$identifier=~s/^[^\|]*\|/p\|/g;
			}
			elsif($Type eq "")
			{
				$Type="DNAMutation";
			}
			if($identifier_hash{"A"}=~/(c|r|m|g|C)/)
			{
				$Type="DNAMutation";
				$identifier=~s/^[^\|]*\|/$identifier_hash{"A"}\|/g;
			}
			elsif($identifier_hash{"A"}=~/p/)
			{
				$Type="ProteinMutation";
			}	
		}
		if($Type eq "ProteinMutation")
		{
			$identifier=~s/^[^\|]*\|/p\|/g;
		}
		
		if(($identifier_hash{"M"} eq $identifier_hash{"W"}) && ($identifier_hash{"W"} ne "") && $Type eq "DNAMutation") # remove genotype
		{
		}
		elsif($NotAaminoacid eq "N" && $Type eq "ProteinMutation") #E 243 ASPARTATE
		{
		}
		elsif($identifier_hash{"W"}=~/,/ && $identifier_hash{"M"}=~/,/ && $identifier_hash{"P"} eq "") #T,C/T,C
		{
		}
		elsif($identifier_hash{"W"} eq "" && $identifier_hash{"M"} eq "" && $identifier_hash{"T"} eq "" && $Type ne "SNP") #exons 5
		{
		}
		elsif($mention_tmp!~/^I[RSrs][Ss]/ && $Type eq "SNP") #exons 5
		{
		}
		else
		{
			$STR_hash{$str_count}=$Type."	".$identifier;
		}
		$str_count++;
	}
		
	my %sentence_hash=();
	my %article_hash=();
	my %printSTR_hash=();
	my %TypeOrder_hash=();
	my $Infer_Refseq="";
	my $pmid="";
	my @pmidARR;
	open input,"<".$PubTator_input;
	$count=1;
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^([^\t\|]+)\|([^\t\|]+)\|(.*)$/)
		{
			$pmid=$1;
			$type=$2;
			$sentence=$3;
			$sentence=~s/[\n\r]//g;
			$sentence_hash{$pmid."\t".$type}=$sentence;
			$article_hash{$pmid}=$article_hash{$pmid}.$sentence." ";
			$TypeOrder_hash{$pmid}=$TypeOrder_hash{$pmid}.$type."\t";
			push(@pmidARR,$pmid);
		}
		elsif($tmp=~/^([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)/)
		{
			my $pmid=$1;
			my $start=$2;
			my $last=$3;
			my $mention=$4;
			my $type=$5;
			if(exists $STR_hash{$count})
			{
				if($STR_hash{$count}=~/^(SNP|DNAMutation|ProteinMutation)	([cgmr])\|/){$Infer_Refseq=$2;}
				$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$STR_hash{$count}."\n";
			}
			$count++;
		}
		else
		{
			my @split_STR=split("\n",$printSTR_hash{$pmid});
			$printSTR_hash{$pmid}="";
			foreach my $STR(@split_STR)
			{
				$STR=~s/[\n\r]//g;
				if($STR=~/^(.+)	(.+)	(.+)	(.+)	(.+)	(.+)/)
				{
					my $start=$2;
					my $last=$3;
					my $mention=$4;
					my $type=$5;
					my $id=$6;
					if($Infer_Refseq ne "" && $id=~/^\|[^\|]*\|[^\|]+\|/)
					{
						$id=$Infer_Refseq.$id;
					}
					
					if($type eq "ProteinMutation" && $setup_hash{ProteinMutation} ne "False")
					{
						if($id=~/^[^\|]*\|([^\|]*)\|([^\|]*)\|([^\|]*)/)#c|T|799+2|A
						{
							my $tmp1=$1;my $tmp2=$2;my $tmp3=$3;
							if($setup_hash{Substitution} ne "False" && $tmp1 eq "SUB")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Deletion} ne "False" && $tmp1 eq "DEL")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Insertion} ne "False" && $tmp1 eq "INS")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Duplication} ne "False" && $tmp1 eq "DUP")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{INDEL} ne "False" && $tmp1 eq "INDEL")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Frameshift} ne "False" && $tmp1 eq "FS")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
						}
					}
					elsif($type eq "DNAMutation" && $setup_hash{DNAMutation} ne "False")
					{
						if($id=~/^[^\|]*\|([^\|]*)\|([^\|]*)\|([^\|]*)/)#c|T|799+2|A
						{
							my $tmp1=$1;my $tmp2=$2;my $tmp3=$3;
							if($setup_hash{Substitution} ne "False" && $tmp1 eq "SUB")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Deletion} ne "False" && $tmp1 eq "DEL")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Insertion} ne "False" && $tmp1 eq "INS")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Duplication} ne "False" && $tmp1 eq "DUP")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{INDEL} ne "False" && $tmp1 eq "INDEL")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
							elsif($setup_hash{Frameshift} ne "False" && $tmp1 eq "FS")
							{
								$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
							}
						}
					}
					elsif($type eq "SNP" && $setup_hash{SNP} ne "False")
					{
						$printSTR_hash{$pmid}=$printSTR_hash{$pmid}.$pmid."	".$start."	".$last."	".$mention."	".$type."	".$id."\n";
					}
				}
			}
			$Infer_Refseq="";
		}
	}
	close input;
	
	my %mention_hash=();
	my %mention2id_hash=();
	my %mention2type_hash=();
	foreach my $pmid(keys %printSTR_hash)
	{
		my @split_STR=split("\n",$printSTR_hash{$pmid});
		foreach my $STR(@split_STR)
		{
			$STR=~s/[\n\r]//g;
			if($STR=~/^(.+)	(.+)	(.+)	(.+)	(.+)	(.+)/)
			{
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $type=$5;
				my $id=$6;
				$mention_hash{$pmid."\t".$mention}=length($mention);
				$mention2id_hash{$pmid."\t".$mention}=$id;
				$mention2type_hash{$pmid."\t".$mention}=$type;
			}
		}
	}
	
	#Read BioC
	my $input_collection = new BioC_full::Collection();
	my $input_xml = new BioC_full::Connector_libxml();
	$input_xml->start_read($BioC_input, $input_collection);
	
	#Write BioC
	my $annotation_collection = new BioC_full::Collection();	# create annotation collection object
	my $annotate_xml = new BioC_full::Connector_libxml();	# create libxml Connector for output
	$annotation_collection ->{source}=$input_collection->{source};
	$annotation_collection ->{date}=$input_collection->{date};
	$annotation_collection ->{key}=$input_collection->{key};
	$annotate_xml->start_write( $BioC_output, $annotation_collection );	# start to write contents of annotation collection to stdout
	
	my $document = new BioC_full::Document();
	while ( $input_xml->read_next($document) ) # documents
	{
		my $document_output = new BioC_full::Document();
		my $pmid = $document->{id};
		$document_output ->{id}=$pmid;
		my $count=0;
		for(my $i=0;$i<$document->{passages}->size();$i++) # passages
		{
			my $psg = new BioC_full::Passage();	# create Passage object
			$psg->{infons}->set("type", $document->{passages}->get($i)->{infons}->get("type"));	# copy type information
			$psg->{offset} = $document->{passages}->get($i)->{offset};	# copy offset element contents
			$psg->{text} = $document->{passages}->get($i)->{text};	# copy text
			
			##====
			##Added existed annotation from input to io
			#for(my $j=0;$j<$document->{passages}->get($i)->{annotations}->size();$j++)	
			#{
			#	$psg->{annotations}->push( $document->{passages}->get($i)->{annotations}->get($j) );
			#	$count++;
			#}
			#for(my $j=0;$j<$document->{passages}->get($i)->{relations}->size();$j++)
			#{
			#	$psg->{relations}->push( $document->{passages}->get($i)->{relations}->get($j) );
			#}
			
			#====
			#insert entities into BioC pasaages
			my @tmpMen = reverse sort {$mention_hash{$a} <=> $mention_hash{$b}} keys %mention_hash;
			foreach my $tmpMen (@tmpMen)
			{
				if($tmpMen=~/^($pmid)	(.+)/)
				{
					my $mention=$2;
					my $mention_tmp=$2;
					$mention_tmp=~s/[\+\-\(\)\[\]]/\[\\W\\\-\\\_\]/g;
					my $sentence=$psg->{text};
					while($sentence=~/^(.*)($mention_tmp)(.*)/)
					{
						my $pre=$1;
						my $mention_tmp=$2;
						my $post=$3;
						
						my $annotation = new BioC_full::Annotation(); #define annotation
						$annotation->{id} = "".$count;
						$annotation->{infons}->set("type", $mention2type_hash{$tmpMen});
						$annotation->{infons}->set("tmVar", $mention2id_hash{$tmpMen});
						$annotation->add_location( (length($pre)+$psg->{offset}), length($mention) );
						$annotation->{text} = $mention;
						$psg->{annotations}->push( $annotation );	# can I use shift?
						
						$count++;
						$mention_tmp=~s/./@/g;
						$sentence=$pre.$mention_tmp.$post;
					}
				}
			}
			$document_output->{passages}->push( $psg );
		}
		
		##====
		##Added existed annotation from input to io
		#for(my $j=0;$j<$document->{relations}->size();$j++)
		#{
		#	$document_output->{relations}->push( $document->{relations}->get($j) );
		#}
		
		$annotate_xml->write_next( $document_output );
	}
	
	$annotate_xml->end_write();
	
	return 1;
}
return 1;
