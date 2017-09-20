#===================================================
# Author:Chih-Hsuan Wei
# Issue: ME
# Description: Post-Processing of Mutation Extraction
#===================================================

package ME;

use utf8;

sub PostProcessing
{
	my $input=$_[0];
	my $output=$_[1];
	my $RegEx=$_[2];

	my %sentence_hash=();
	my %TypeOrder_hash=();
	my %annotation_hash=();
	my @pmidlist;
	my $article="";
	my %mention_pattern_hash=();
	my %type_pattern_hash=();
	my %mention_hash=();
	my %type_hash=();
	my %identifier_hash=();

	my $pmid="";
	open input,"<".$input;
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\t\|]+)\|([^\t\|]+)\|(.*)$/)
		{
			$pmid=$1;
			$type=$2;
			$sentence=$3;
			$sentence=~s/[\n\r]//g;
			$sentence_hash{$pmid."\t".$type}=$sentence;
			$article=$article.$sentence." ";
			$TypeOrder_hash{$pmid}=$TypeOrder_hash{$pmid}.$type."\t";
			push(@pmidlist,$pmid);
		}
		elsif($tmp=~/^([^\t]+)	([0-9]+)	([0-9]+)	(.+)	(SNP|DNAMutation|ProteinMutation)	(.*)/ || $tmp=~/^([^\t]+)	([0-9]+)	([0-9]+)	(.+)	(SNP|DNAMutation|ProteinMutation)/)
		{
			$pmid=$1;
			my $start=$2;
			my $last=$3;
			my $mention=$4;
			my $type=$5;
			my $identifier=$6;
			my @id_split=split(/\|/,$identifier);
			$mention_hash{$start."\t".$last}=$mention;
			$type_hash{$start."\t".$last}=$type;
			$identifier_hash{$start."\t".$last}=$identifier;
		}
		else
		{
			$article=" ".$article;

			#=====
			#Split RSnumber
			foreach my $posit (keys %mention_hash)
			{
				#print $mention_hash{$posit}."\n";
				if($mention_hash{$posit}=~/^([RrSs][Ss][ ]*[0-9]+)[ ]*(and|\/|,|or)[ ]*([RrSs][Ss][ ]*[0-9]+)$/)
				{
					my $tmp1=$1;
					my $tmp2=$3;
					my ($start,$last)=($posit=~/(.+)\t(.+)/);
					my $start1=$start;
					my $last1=$start1+length($tmp1);
					my $start2=$last-length($tmp2);
					my $last2=$last;
					$mention_hash{$start1."\t".$last1}=$tmp1;
					$type_hash{$start1."\t".$last1}="SNP";
					$mention_hash{$start2."\t".$last2}=$tmp2;
					$type_hash{$start2."\t".$last2}="SNP";
					delete $mention_hash{$posit};
					delete $type_hash{$posit};
				}
				elsif($mention_hash{$posit}=~/^(.*[^0-9])[ ]*([RrSs][Ss][ ]*[0-9]+)$/)
				{
					my $tmp1=$1;
					my $tmp2=$2;
					my ($start,$last)=($posit=~/(.+)\t(.+)/);
					my $start1=$start;
					my $last1=$start1+length($tmp1);
					my $start2=$last-length($tmp2);
					my $last2=$last;
					$mention_hash{$start1."\t".$last1}=$tmp1;
					$type_hash{$start1."\t".$last1}="DNAMMutation";
					$mention_hash{$start2."\t".$last2}=$tmp2;
					$type_hash{$start2."\t".$last2}="SNP";
					delete $mention_hash{$posit};
					delete $type_hash{$posit};
				}
				elsif($mention_hash{$posit}=~/^([RrSs][Ss][ ]*[0-9]+)[ ]*([^0-9].*)$/)
				{
					my $tmp1=$1;
					my $tmp2=$2;
					my ($start,$last)=($posit=~/(.+)\t(.+)/);
					my $start1=$start;
					my $last1=$start1+length($tmp1);
					my $start2=$last-length($tmp2);
					my $last2=$last;
					$mention_hash{$start1."\t".$last1}=$tmp1;
					$type_hash{$start1."\t".$last1}="SNP";
					$mention_hash{$start2."\t".$last2}=$tmp2;
					$type_hash{$start2."\t".$last2}="DNAMMutation";
					delete $mention_hash{$posit};
					delete $type_hash{$posit};
				}
			}

			#=====
			#boundary
			foreach my $posit(keys %mention_hash)
			{
				my ($start,$last)=($posit=~/(.+)\t(.+)/);
				my $mention=$mention_hash{$posit};
				if($mention=~/^[0-9]/ && substr($article,$start,1)=~/([\-\+])/)	#17000021	251	258	1858C>T --> +1858C>T
				{
					delete $mention_hash{$posit};
					$mention=$1.$mention;
					$mention_hash{($start-1)."\t".$last}=$mention;
					delete $mention_hash{$posit};
					$type_hash{($start-1)."\t".$last}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{($start-1)."\t".$last}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				if($mention=~/^[\W\-\_][RrSs][Ss][0-9]+/)	#_rs7207916
				{
					delete $mention_hash{$posit};
					$mention=~s/^\(//g;
					$mention_hash{($start+1)."\t".$last}=$mention;
					delete $mention_hash{$posit};
					$type_hash{($start+1)."\t".$last}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{($start+1)."\t".$last}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				if($mention=~/^\(/ && $mention!~/\)/) # delete (
				{
					delete $mention_hash{$posit};
					$mention=~s/^\(//g;
					$mention_hash{($start+1)."\t".$last}=$mention;
					delete $mention_hash{$posit};
					$type_hash{($start+1)."\t".$last}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{($start+1)."\t".$last}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				elsif($mention=~/\)$/ && $mention!~/\(/) # delete )
				{
					delete $mention_hash{$posit};
					$mention=~s/\)$//g;
					$mention_hash{$start."\t".($last-1)}=$mention;
					delete $mention_hash{$posit};
					$type_hash{$start."\t".($last-1)}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{$start."\t".($last-1)}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				elsif(substr($article,$start,1) eq "(" && $mention=~/\)/) # add (
				{
					delete $mention_hash{$posit};
					$mention="(".$mention;
					$mention_hash{($start-1)."\t".$last}=$mention;
					delete $mention_hash{$posit};
					$type_hash{($start-1)."\t".$last}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{($start-1)."\t".$last}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				elsif(substr($article,$last+1,1) eq ")" && $mention=~/\(/) # add )
				{
					delete $mention_hash{$posit};
					$mention=$mention.")";
					$mention_hash{$start."\t".($last+1)}=$mention;
					delete $mention_hash{$posit};
					$type_hash{$start."\t".($last+1)}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{$start."\t".($last+1)}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
				elsif($mention=~/^\(/ && $mention=~/\)$/) # delete (  )
				{
					delete $mention_hash{$posit};
					$mention=~s/^\(//g;
					$mention=~s/\)$//g;
					$mention_hash{($start+1)."\t".($last-1)}=$mention;
					delete $mention_hash{$posit};
					$type_hash{($start+1)."\t".($last-1)}=$type_hash{$posit};
					delete $type_hash{$posit};
					$identifier_hash{($start+1)."\t".($last-1)}=$identifier_hash{$posit};
					delete $identifier_hash{$posit};
				}
			}

			#=====
			#Self-pattern
			foreach my $posit(keys %mention_hash) #build pattern
			{
				if($type_hash{$posit} ne "SNP")
				{
					my $mention_pattern=$mention_hash{$posit};
					my $mention=$mention_hash{$posit};
					$mention_pattern=~s/([\W\-\_])/\\$1/g;
					$mention_pattern=~s/[0-9]+/\[0\-9\]\+/g;
					$mention_pattern=~s/(IVS|EX)/@@@@/g;
					$mention_pattern=~s/(rs|ss)/@@@/g;
					$mention_pattern=~s/[a-z]/\[a\-z\]/g;
					$mention_pattern=~s/[A-Z]/\[A\-Z\]/g;
					$mention_pattern=~s/@@@@/\(IVS\|EX\)/g;
					$mention_pattern=~s/@@@/rs/g;
					$mention_pattern_hash{$mention_pattern}=length($mention);
					$type_pattern_hash{$mention_pattern}=$type_hash{$posit};
				}
			}

			#=====
			#Append-pattern
			my $countRegEx=1;
			open RegEx,"<$RegEx";
			while(<RegEx>)
			{
				my $tmp=$_;
				$tmp=~s/[\n\r]//g;
				if($tmp=~/(.+)	(.+)/)
				{
					my $RegEx=$1;
					my $type=$2;
					$mention_pattern_hash{$RegEx}=$countRegEx;
					$type_pattern_hash{$RegEx}=$type;
					$countRegEx++;
				}
			}
			close RegEx;

			#=====
			#Extract mutation by pattern match
			my @rank = sort {$mention_pattern_hash{$a} <=> $mention_pattern_hash{$b}} keys %mention_pattern_hash;
			foreach my $mention_pattern(@rank)
			{
				my $article_tmp=$article;
				while($article_tmp=~/^(.+[\W\-\_])($mention_pattern)([\W\-\_].+)$/ || $article_tmp=~/^(.+[\W\-\_a-z0-9])($mention_pattern)([\W\-\_].+)$/)
				{
					my $pre=$1;
					my $mention=$2;
					my $post=$3;
					my $start=length($pre)-1;
					my $last=length($pre)+length($mention)-1;
					if($article_tmp=~/^(.+[\W\-\_])($mention_pattern)([\W\-\_].+)$/ || ($article_tmp=~/^(.+[\W\-\_a-z0-9])($mention_pattern)([\W\-\_].+)$/ && $mention=~/^[A-Z]/))
					{
						if(not exists $mention_hash{$start."\t".$last})
						{
							$mention_hash{$start."\t".$last}=$mention;
							$type_hash{$start."\t".$last}=$type_pattern_hash{$mention_pattern};
						}
						$mention=~s/./@/g;
						$article_tmp=$pre.$mention.$post;
					}
					else
					{
						$mention=~s/./@/g;
						$article_tmp=$pre.$mention.$post;
					}
				}
			}

			#=====
			#filter cover part
			foreach my $positA (keys %mention_hash)
			{
				my ($startA,$lastA)=($positA=~/(.+)\t(.+)/);
				foreach my $positB (keys %mention_hash)
				{
					my ($startB,$lastB)=($positB=~/(.+)\t(.+)/);
					if(($startA>$startB && $lastA<=$lastB) || ($startA>=$startB && $lastA<$lastB))
					{
						delete $mention_hash{$positA};
						delete $type_hash{$positA};
					}
				}
			}

			#=====
			#filter short mention
			foreach my $posit (keys %mention_hash)
			{
				#print $mention_hash{$posit}."\n";
				if($mention_hash{$posit}=~/^[A-Z][0-9][A-Z]$/)
				{
					delete $mention_hash{$posit};
					delete $type_hash{$posit};
				}
			}

			my %start_hash=();
			foreach my $posit(keys %mention_hash)
			{
				my ($start,$last)=($posit=~/(.+)\t(.+)/);
				$start_hash{$start*10000+$last}=$posit;
			}

			my @rank = sort {$start_hash{$a} <=> $start_hash{$b}} keys %start_hash;
			foreach my $start(@rank)
			{
				my $posit=$start_hash{$start};
				$annotation_hash{$pmid}=$annotation_hash{$pmid}.$pmid."\t".$posit."\t".$mention_hash{$posit}."\t".$type_hash{$posit}."\t".$identifier_hash{$posit}."\n";
			}

			%mention_pattern_hash=();
			%type_pattern_hash=();
			%mention_hash=();
			%type_hash=();
			%identifier_hash=();
			$article="";
		}
	}
	close input;

	open output,">".$output;
	my %printed_hash=();
	foreach my $pmid (@pmidlist)
	{
		if(not exists $printed_hash{$pmid})
		{
			my @type=split("\t",$TypeOrder_hash{$pmid});
			foreach my $type(@type)
			{
				print output $pmid."|".$type."|".$sentence_hash{$pmid."\t".$type}."\n";
			}
			print output $annotation_hash{$pmid}."\n";
			$printed_hash{$pmid}=1;
		}
	}
	close output;
	return 1;
}
sub ExtractStateFeature
{
	my $input=$_[0];
	my $output=$_[1];
	my $filename=$_[2];

	open input,"<".$input;
	my $pmid="";
	my %RegEx_locationmap_hash=();
	my %RegEx_lastmap_hash=();
	my %RegEx_mention_hash=();
	my %RegEx_type_hash=();
	my %RegEx_result_hash=();

	my %Result_hash=();
	my %locationmap_hash=();
	my %lastmap_hash=();
	my %mention_hash=();
	my %type_hash=();
	my %character_hash=();

	open output,">".$output;
	my $mode_RegEx="";
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\t\|]+)\|([^\t\|]+)\|(.*)$/)
		{
		}
		elsif($tmp=~/^([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)/)
		{
			%Result_hash=();
			%locationmap_hash=();
			%lastmap_hash=();
			%mention_hash=();
			%type_hash=();
			%character_hash=();

			$pmid=$1;
			$start=1;
			$last=$3;
			$mention=$4;
			$type=$5;
			$identifier=$6;
			if($identifier=~/^([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|(fs[^\|]*)\|([^\|]*)$/i)
			{
				my @type=split(",",$1);
				my @W=split(",",$2);
				my @P=split(",",$3);
				my @M=split(",",$4);
				my @F=split(",",$5);
				my @S=split(",",$6);
				my $tmpmention=$mention;
				foreach my $type(@type){	if($mention=~/^(.*?)($type)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="A";	}	}	}
				foreach my $W(@W){	if($tmpmention=~/^(.*)($W)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="W"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $P(@P){	$P=~s/([\W\_\-])/\\$1/g;	if($tmpmention=~/^(.*)($P)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="P"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $M(@M){	if($tmpmention=~/^(.*)($M)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="M"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $F(@F){	$F=~s/([\W\_\-])/\\$1/g;	if($mention=~/^(.*)($F)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="F";	}	}	}
				foreach my $S(@S){	if($mention=~/^(.*)($S)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="S";	}	}	}
			}
			elsif($identifier=~/^([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|(fs[^\|]*)$/i)
			{
				my @type=split(",",$1);
				my @W=split(",",$2);
				my @P=split(",",$3);
				my @M=split(",",$4);
				my @F=split(",",$5);
				my $tmpmention=$mention;
				foreach my $type(@type){	if($mention=~/^(.*?)($type)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="A";	}	}	}
				foreach my $W(@W){	if($tmpmention=~/^(.*)($W)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="W"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $P(@P){	$P=~s/([\W\_\-])/\\$1/g;	if($tmpmention=~/^(.*)($P)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="P"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $M(@M){	if($tmpmention=~/^(.*)($M)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="M"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $F(@F){	$F=~s/([\W\_\-])/\\$1/g;	if($mention=~/^(.*)($F)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="F";	}	}	}
			}
			elsif($identifier=~/^([^\|]*)\|([^\|]*(ins|del|Del|dup|\-)[^\|]*)\|([^\|]*)\|([^\|]*)$/i) #substitution,insertion, deletion
			{
				my @type=split(",",$1);
				my @T=split(",",$2);
				my @P=split(",",$4);
				my @M=split(",",$5);
				my $tmpmention=$mention;
				foreach my $type(@type){	if($mention=~/^(.*?)($type)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="A";	}	}	}
				foreach my $P(@P){	$P=~s/([\W\_\-])/\\$1/g;	if($tmpmention=~/^(.*)($P)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="P"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $T(@T){	if($tmpmention=~/^(.*?)($T)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="T"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $M(@M){	if($tmpmention=~/^(.*)($M)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="M"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
			}
			elsif($identifier=~/^([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)$/i)
			{
				my @type=split(",",$1);
				my @W=split(",",$2);
				my @P=split(",",$3);
				my @M=split(",",$4);
				my $tmpmention=$mention;
				foreach my $type(@type){	if($mention=~/^(.*?)($type)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="A";	}	}	}
				foreach my $W(@W){	if($tmpmention=~/^(.*)($W)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="W"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $P(@P){	$P=~s/([\W\_\-])/\\$1/g;	if($tmpmention=~/^(.*)($P)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="P"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $M(@M){	if($tmpmention=~/^(.*)($M)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="M"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
			}
			elsif($identifier=~/^([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)$/i) # dup
			{
				my @type=split(",",$1);
				my @T=split(",",$2);
				my @P=split(",",$3);
				my @M=split(",",$4);
				my @D=split(",",$5);
				my $tmpmention=$mention;
				foreach my $type(@type){	if($mention=~/^(.*?)($type)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="A";	}	}	}
				foreach my $P(@P){	$P=~s/([\W\_\-])/\\$1/g;	if($tmpmention=~/^(.*)($P)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="P"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $T(@T){	if($tmpmention=~/^(.*)($T)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="T"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $M(@M){	if($tmpmention=~/^(.*)($M)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="M"; my $str1=$1;my $str2=$2;my $str3=$3;$str2=~s/./@/g;$tmpmention=$str1.$str2.$str3;	}	}	}
				foreach my $D(@D){
					$D=~s/([\W\_\-])/\\$1/g;
					if($mention=~/^(.*)($D)(.*)$/)
					{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="D";	}	}
				}
			}
			elsif($identifier=~/^((rs|RS|Rs|reference SNP no\. )[0-9]+)$/i)
			{
				my @RS=split(",",$1);
				foreach my $RS(@RS){	if($mention=~/^(.*)($RS)(.*)$/)	{	for(my $i=length($1);$i<(length($1)+length($2));$i++){$character_hash{$i}="R";	}	}	}
			}

			#=====
			#Start state:I
			print output "I I I I I I I I I I I I I I I I I I\n";

			#=====
			#Translate all tokens to states
			my $tmp=$mention;
			my $start=0;
			while($tmp=~/^([A-Z]+)(.*)$/ || $tmp=~/^([a-z]+)(fs.*)$/ || $tmp=~/^([a-z]+)(.*)$/ || $tmp=~/^([0-9]+)(.*)$/ || $tmp=~/^([\W\-\_])(.*)$/ )
			{
				my $pre=$1;
				my $post=$2;
				my $end=length($pre)+$start;

				if ($pre ne " ")
				{
					my $state="I";
					if (exists $character_hash{$end-1}){$state=$character_hash{$end-1};}

					my $token=$pre;

					#Number of Numbers [0-9]
					my $tmp=$token;
					$tmp=~s/[^0-9]//g;
					my $Num_num="";
					if(length($tmp)>2){$Num_num="N:3+";}elsif(length($tmp)<=2 && length($tmp)>=1){$Num_num="N:1-2";}else{$Num_num="N:".length($tmp);}

					#Number of Uppercase [A-Z]
					$tmp=$token;
					$tmp=~s/[^A-Z]//g;
					my $Num_Uc="";
					if(length($tmp)>3){$Num_Uc="U:3+";}elsif(length($tmp)<=2 && length($tmp)>=1){$Num_Uc="N:1-2";}else{$Num_Uc="U:".length($tmp);}

					#Number of Lowercase [a-z]
					$tmp=$token;
					$tmp=~s/[^a-z]//g;
					my $Num_Lc="";
					if(length($tmp)>3){$Num_Lc="L:3+";}elsif(length($tmp)<=2 && length($tmp)>=1){$Num_Lc="N:1-2";}else{$Num_Lc="L:".length($tmp);}

					#Number of ALL char
					$tmp=$token;
					$tmp=~s/[^a-z]//g;
					my $Num_All="";
					if(length($tmp)>3){$Num_All="A:3+";}elsif(length($tmp)<=2 && length($tmp)>=1){$Num_All="N:1-2";}else{$Num_All="A:".length($tmp);}

					#specific character (;:,.->+_)
					$tmp=$token;
					my $SpecificC="";
					if($tmp=~/[\;\:\,\.\-\>\+\_]/){$SpecificC="-SpecificC1-";}
					elsif($tmp=~/[\(\)]/){$SpecificC="-SpecificC2-";}
					elsif($tmp=~/[\{\}]/){$SpecificC="-SpecificC3-";}
					elsif($tmp=~/[\[\]]/){$SpecificC="-SpecificC4-";}
					elsif($tmp=~/[\\\/]/){$SpecificC="-SpecificC5-";}
					else{$SpecificC="__nil__";}

					#mutation level
					my $Mlevel="";
					if($token eq "p"){$Mlevel="-ProteinLevel-";}
					elsif($token=~/^[cgmr]$/){$Mlevel="-DNALevel-";}
					else{$Mlevel="__nil__";}

					#mutation type
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $Mtype="";
					if($lc_tmp=~/^(deletion|insertion|duplication|repeat|inversion|delta)$/){$Mtype="-Mtype- -MtypeFull-";}
					elsif($lc_tmp=~/^(eletion|elta|nsertion|uplication|epeat|nversion)$/){$Mtype="-Mtype- -MtypeFull_suffix-";}
					elsif($lc_tmp=~/^(del|ins|delins|indel|dup|inv)$/){$Mtype="-Mtype- -MtypeTri-";}
					elsif($lc_tmp eq "\/" && substr($mention,$i-3,3)=~/(ins|del)/i && substr($mention,$i+1,3)=~/(ins|del)/i){$Mtype="-Mtype- -MtypeTri-";}
					elsif(substr($mention,$i-1,2) eq "\/\\" || substr($mention,$i,2) eq "\/\\"){$Mtype="-Mtype- -MtypeTri-";}
					else{$Mtype="__nil__ __nil__";}

					#DNA symbols
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $DNASym="";
					if($lc_tmp=~/^(adenine|guanine|thymine|cytosine)$/){$DNASym="-DNASym- -DNASymFull-";}
					elsif($lc_tmp=~/^(denine|uanine|hymine|ytosine)$/){$DNASym="-DNASym- -DNASymFull_suffix-";}
					elsif($lc_tmp=~/^[atcgu]+$/){$DNASym="-DNASym- -DNASymChar-";}
					else{$DNASym="__nil__ __nil__";}

					#Protein symbols
					my $uc_tmp=$token;
					my $lc_tmp=lc($token);
					my $ProteinSym="";
					if($lc_tmp=~/^(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift)$/){$ProteinSym="-ProteinSym- -ProteinSymFull-";}
					elsif($lc_tmp=~/^(lutamine|lutamic|eucine|aline|soleucine|ysine|lanine|lycine|spartate|ethionine|hreonine|istidine|spartic|sparticacid|rginine|sparagine|ryptophan|roline|henylalanine|ysteine|erine|lutamate|yrosine|top|rameshift)$/){$ProteinSym="-ProteinSym- -ProteinSymFull_suffix-";}
					elsif($lc_tmp=~/^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr)$/){$ProteinSym="-ProteinSym- -ProteinSymTri-";}
					elsif($lc_tmp=~/^(ys|le|er|ln|et|sn|ro|ys|sp|hr|phe|la|ly|is|eu|rg|rp|al|lu|yr)$/){$ProteinSym="-ProteinSym- -ProteinSymTri_suffix-";}
					elsif($token=~/^[CISQMNPKDTFAGHLRWVEYX]$/ && $next!~/^[ylsrhpieraYLSRHPIERA]$/){$ProteinSym="-ProteinSym- -ProteinSymChar-";}
					elsif($token=~/^[CISGMPLTHAVF]$/ && $next=~/^[ylsrhpieraYLSRHPIERA]$/){$ProteinSym="-ProteinSym- -ProteinSymChar-";}
					else{$ProteinSym="__nil__ __nil__";}

					#IVS/EX
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $IVSEX="";
					if($lc_tmp=~/^(ivs|ex)$/){$IVSEX="-IVSEX-";}
					elsif($tmp eq "E" && substr($mention,$i+1,1) eq "x"){$IVSEX="-IVSEX-";}
					elsif(substr($mention,$i-1,1) eq "E" && $tmp eq "x"){$IVSEX="-IVSEX-";}
					else{$IVSEX="__nil__";}

					#FSX feature
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $FSXfeature="";
					if($lc_tmp=~/^(fs|fsx|x|\*)$/){$FSXfeature="-FSX-";}
					elsif(substr($mention,$i-1,1)=~/[sS]/ && $tmp eq "X"){$FSXfeature="-FSX-";}
					else{$FSXfeature="__nil__";}

					#position type
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $PositionType="";
					if($lc_tmp=~/^(nucleotide|codon|amino|acid|position|bp|b)$/){$PositionType="-PositionType-";}
					else{$PositionType="__nil__";}

					#sequence location
					$tmp=$token;
					my $lc_tmp=lc($tmp);
					my $SeqLocat="";
					if($lc_tmp=~/^(intron|exon|promoter|utr)$/){$SeqLocat="-SeqLocat-";}
					else{$SeqLocat="__nil__";}

					#RS
					$tmp=lc($token);
					my $RScode="";
					if($tmp eq "rs"){$RScode="-RScode-";}
					else{$RScode="__nil__";}

					my $result="";
					if($wildtype ne "" && $i>=$wildtype_start && $i<$wildtype_last)
					{
						$result=$Inside_hash{"WildType"};
					}
					elsif($mutant ne "" && $i>=$mutant_start && $i<$mutant_last)
					{
						$result=$Inside_hash{"Mutant"};
					}
					elsif($position ne "" && $i>=$position_start && $i<$position_last)
					{
						$result=$Inside_hash{"position"};
					}
					elsif($FSX ne "" && $i>=$FSX_start && $i<$FSX_last)
					{
						$result=$Inside_hash{"FSX"};
					}
					else
					{
						$result="O";
					}

					print output "$token $Num_num $Num_Uc $Num_Lc $Num_All $SpecificC $Mlevel $Mtype $DNASym $ProteinSym $IVSEX $FSXfeature $PositionType $SeqLocat $RScode\n";

				}
				$tmp=$post;
				if ($char ne "")
				{
					$start=$end+1;
				}
				else
				{
					$start=$end;
				}
			}
			print output "                  \n";
		}
	}
	close input;
	close output;

	open(CRF,"|crf_test -m CRF/ComponentExtraction.Model -o tmp/".$filename.".PostME.output tmp/".$filename.".PostME.data")|| die ("can't open\n");	close(CRF);

	return 1;
}
return 1;
