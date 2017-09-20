#===================================================
# Author:Chih-Hsuan Wei
# Issue: ME
# Description: Extracting Features
#===================================================

package ME;

require 'lib/Text/English.pm';
#require 'lib/Lingua/EN/Tagger.pm';
#my $p = new Lingua::EN::Tagger;
use encoding 'big5', STDIN => 'big5', STDOUT => 'big5';

sub ExtractFeature_PubTator
{
	my $input=$_[0];
	my $output=$_[1];
	my $output_loc=$_[2];

	my @stemmed_tokens;
	my @POS_tokens;
	my %Begin_hash=();	
	my %Inside_hash=();
	my %End_hash=();
	$Begin_hash{"DNAMutation"}="D";			$Inside_hash{"DNAMutation"}="E";			$End_hash{"DNAMutation"}="F";
	$Begin_hash{"ProteinMutation"}="P";		$Inside_hash{"ProteinMutation"}="Q";		$End_hash{"ProteinMutation"}="R";
	$Begin_hash{"ChromosomalMutation"}="C";	$Inside_hash{"ChromosomalMutation"}="B";	$End_hash{"ChromosomalMutation"}="A";
	$Begin_hash{"SNP"}="S";					$Inside_hash{"SNP"}="T";					$End_hash{"SNP"}="U";
	
	my %article_hash=();
	my $title="";
	my $abstract="";
	my $article="";
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
	close output;
	open locationA,">".$output_loc;
	close locationA;
	my $mode_RegEx="";
	my $mode=0;
	open input,"<".$input;
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$tmp=~s/[^0-9a-zA-Z\~\!\@\#\$\%\^\&\*\(\)\_\+\{\}\|\:\"\<\>\?\`\-\=\[\]\\\;\'\,\.\/\t ]/ /g;
		if($tmp=~/^([^\t\|]+)\|([^\t\|]+)\|(.*)$/)
		{
			$pmid=$1;
			$type=$2;
			$sentence=$3;
			$sentence=~s/[\n\r]//g;
			$article_hash{$pmid."\t".$type}=$sentence;
			$article=$article.$sentence." ";
		}
		elsif($tmp=~/^([0-9]+)\t([^\t]+)\t([^\t]+)$/)
		{
			$pmid=$1;
			$title=$2;
			$abstract=$3;
			$title=~s/[\n\r]//g;$abstract=~s/[\n\r]//g;
			$article_hash{$pmid."\tt"}=$title;
			$article_hash{$pmid."\ta"}=$abstract;
			$article=$title." ".$abstract." ";
			$mode=1;
		}
		elsif($tmp=~/^([0-9]+)\t([^\t]+)$/)
		{
			$pmid=$1;
			$type="t";
			$sentence=$2;
			$sentence=~s/[\n\r]//g;
			$article_hash{$pmid."\t".$type}=$sentence;
			$article=$article.$sentence." ";
			$mode=1;
		}
		if(length($tmp)<=1 || $mode==1)
		{
			my $count_post=0;
			#RegEx
			my $article_tmp=" ".$article;
			while ($article_tmp=~/^(.*\W)([cgrm]\.[ATCGatcgu \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(inv|del|ins|dup|tri|qua|con|delins|indel)[ATCGatcgu0-9\_\.\:]*)(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)(IVS[ATCGatcgu \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(del|ins|dup|tri|qua|con|delins|indel)[ATCGatcgu0-9\_\.\:]*)(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([cgrm]\.[ATCGatcgu \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)(IVS[ATCGatcgu \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([cgrm]\.[ATCG][0-9]+[ATCG])()(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([ATCGU][0-9]+[ATCG])()(\W.*)$/
			|| $article_tmp=~/^(.*\W)([0-9]+(del|ins|dup|tri|qua|con|delins|indel)[ATCG]*)()(\W.*)$/
			|| $article_tmp=~/^(.*\W)(Delta[0-9]+)()(\W.*)$/
			)
			{
				my $pre=$1;
				$pre=substr($pre,1,);
				my $mention=$2;
				my $tmp=$2;
				my $post=$4;
				$tmp=~s/./@/g;
				$article_tmp=" ".$pre.$tmp.$post;
				if($mention !~/^[cgr][\W\-\_]*\.[\W\-\_]*[0-9]+[\W\-\_]*$/)
				{
					$mention=~s/[\s\;\.\,\:]+$//g;
					if($mention!~/\[/){$mention=~s/\]$//g;}
					if($mention!~/\(/){$mention=~s/\)$//g;}
					if($mention!~/\{/){$mention=~s/\}$//g;}
					if($mention!~/\]/){$mention=~s/^\[//g;}
					if($mention!~/\)/){$mention=~s/^\(//g;}
					if($mention!~/\}/){$mention=~s/^\{//g;}
					$RegEx_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$RegEx_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
					$RegEx_mention_hash{(length($pre)+1)." ".(length($pre)+length($mention))}=$mention;
					$RegEx_type_hash{(length($pre)+1)." ".(length($pre)+length($mention))}="DNAMutation";
				}
				if($count_post==0)
				{
					my $post_tmp=$post;
					while ($post_tmp=~/^(.*\W)([cgrm]\.[ATCGatcgu \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(inv|del|ins|dup|tri|qua|con|delins|indel)[ATCGatcgu0-9\_\.\:]*)(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)(IVS[ATCGatcgu \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(inv|del|ins|dup|tri|qua|con|delins|indel)[ATCGatcgu0-9\_\.\:]*)(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([cgrm]\.[ATCGatcgu \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)(IVS[ATCGatcgu \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([cgrm]\.[ATCG][0-9]+[ATCG])()(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([ATCG][0-9]+[ATCG])()(\W.*)$/
					|| $post_tmp=~/^(.*\W)([0-9]+(del|ins|dup|tri|qua|con|delins|indel)[ATCGU]*)()(\W.*)$/
					)
					{
						my $pre_post=$1;
						$pre_post=substr($pre_post,1,);
						my $mention_post=$2;
						my $tmp=$2;
						my $post_post=$4;
						$tmp=~s/./@/g;
						$post_tmp=" ".$pre_post.$tmp.$post_post;
						if($mention_post !~/^[cgr][\W\-\_]*\.[\W\-\_]*[0-9]+[\W\-\_]*$/)
						{
							$mention_post=~s/[\s\;\.\,\:]+$//g;
							if($mention_post!~/\[/){$mention_post=~s/\]$//g;}
							if($mention_post!~/\(/){$mention_post=~s/\)$//g;}
							if($mention_post!~/\{/){$mention_post=~s/\}$//g;}
							if($mention_post!~/\]/){$mention_post=~s/^\[//g;}
							if($mention_post!~/\)/){$mention_post=~s/^\(//g;}
							if($mention_post!~/\}/){$mention_post=~s/^\{//g;}
							$RegEx_locationmap_hash{length($pre_post)+1}=(length($pre_post)+1)." ".(length($pre_post)+length($mention_post));
							$RegEx_lastmap_hash{length($pre_post)+1}=(length($pre_post)+length($mention_post));
							$RegEx_mention_hash{(length($pre_post)+1)." ".(length($pre_post)+length($mention_post))}=$mention_post;
							$RegEx_type_hash{(length($pre_post)+1)." ".(length($pre_post)+length($mention_post))}="DNAMutation";
						}
					}
				}
				$count_post++;
			}
			
			$count_post=0;
			$article_tmp=" ".$article;
			while ($article_tmp=~/^(.*\W)([p]\.[CISQMNPKDTFAGHLRWVEYX \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(inv|del|ins|dup|tri|qua|con|delins|indel|fsX|fs X|fsx|fs x|fs)[CISQMNPKDTFAGHLRWVEYX \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]*)(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([p]\.[CISQMNPKDTFAGHLRWVEYX \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([p]\.[A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}[A-Z][a-z]{0,2})()(\W.*)$/ 
			|| $article_tmp=~/^(.*\W)([p]\.[A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}(fs|fsx|fsX))(\W.*)$/
			|| $article_tmp=~/^(.*\W)([A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}[A-Z][a-z]{0,2})()(\W.*)$/
			|| $article_tmp=~/^(.*\W)([A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}(fs|fsx|fsX))(\W.*)$/
			)
			{
				my $pre=$1;
				$pre=substr($pre,1,);
				my $mention=$2;
				my $tmp=$2;
				my $post=$4;
				$tmp=~s/./@/g;
				$article_tmp=" ".$pre.$tmp.$post;
				if($mention !~/^[p][\W\-\_]*\.[\W\-\_]*[0-9]+[\W\-\_]*$/)
				{
					$mention=~s/[\s\;\.\,\:]+$//g;
					if($mention!~/\[/){$mention=~s/\]$//g;}
					if($mention!~/\(/){$mention=~s/\)$//g;}
					if($mention!~/\{/){$mention=~s/\}$//g;}
					if($mention!~/\]/){$mention=~s/^\[//g;}
					if($mention!~/\)/){$mention=~s/^\(//g;}
					if($mention!~/\}/){$mention=~s/^\{//g;}
					$RegEx_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$RegEx_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
					$RegEx_mention_hash{(length($pre)+1)." ".(length($pre)+length($mention))}=$mention;
					$RegEx_type_hash{(length($pre)+1)." ".(length($pre)+length($mention))}="ProteinMutation";
				}
				if($count_post==0)
				{
					my $post_tmp=$post;
					while ($post_tmp=~/^(.*\W)([p]\.[CISQMNPKDTFAGHLRWVEYX \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]+(inv|del|ins|dup|tri|qua|con|delins|indel|fsX|fs X|fsx|fs x|fs)[CISQMNPKDTFAGHLRWVEYX \/\>\<\?\(\)\[\]\;\:\*\_\-\+0-9]*)(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([p]\.[CISQMNPKDTFAGHLRWVEYX \/\>\?\(\)\[\]\;\:\*\_\-\+0-9]+)()(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([p]\.[A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}[A-Z][a-z]{0,2})()(\W.*)$/ 
					|| $post_tmp=~/^(.*\W)([p]\.[A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}(fs|fsx|fsX))(\W.*)$/
					|| $post_tmp=~/^(.*\W)([A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}[A-Z][a-z]{0,2})()(\W.*)$/
					|| $post_tmp=~/^(.*\W)([A-Z][a-z]{0,2}[\W\-]{0,1}[0-9]+[\W\-]{0,1}(fs|fsx|fsX))(\W.*)$/
					)
					{
						my $pre_post=$1;
						$pre_post=substr($pre_post,1,);
						my $mention_post=$2;
						my $tmp=$2;
						my $post_post=$4;
						$tmp=~s/./@/g;
						$post_tmp=" ".$pre_post.$tmp.$post_post;
						if($mention_post !~/^[cgr][\W\-\_]*\.[\W\-\_]*[0-9]+[\W\-\_]*$/)
						{
							$mention_post=~s/[\s\;\.\,\:]+$//g;
							if($mention_post!~/\[/){$mention_post=~s/\]$//g;}
							if($mention_post!~/\(/){$mention_post=~s/\)$//g;}
							if($mention_post!~/\{/){$mention_post=~s/\}$//g;}
							if($mention_post!~/\]/){$mention_post=~s/^\[//g;}
							if($mention_post!~/\)/){$mention_post=~s/^\(//g;}
							if($mention_post!~/\}/){$mention_post=~s/^\{//g;}
							$RegEx_locationmap_hash{length($pre_post)+1}=(length($pre_post)+1)." ".(length($pre_post)+length($mention_post));
							$RegEx_lastmap_hash{length($pre_post)+1}=(length($pre_post)+length($mention_post));
							$RegEx_mention_hash{(length($pre_post)+1)." ".(length($pre_post)+length($mention_post))}=$mention_post;
							$RegEx_type_hash{(length($pre_post)+1)." ".(length($pre_post)+length($mention_post))}="ProteinMutation";
						}
					}
				}
				$count_post++;
			}			
			
			$article_tmp=" ".$article;
			while ($article_tmp=~/^(.*\W)([Rr][Ss][0-9]+)(\W.*)$/)
			{
				my $pre=$1;
				$pre=substr($pre,1,);
				my $mention=$2;
				my $tmp=$2;
				my $post=$4;
				$tmp=~s/./@/g;
				$article_tmp=" ".$pre.$tmp.$post;
				$RegEx_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
				$RegEx_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
				$RegEx_mention_hash{(length($pre)+1)." ".(length($pre)+length($mention))}=$mention;
				$RegEx_type_hash{(length($pre)+1)." ".(length($pre)+length($mention))}="SNP";
			}			
			
			#Extract the location for each token
			open locationA,">>".$output_loc;
			$article_tmp = $article;
			my $count_token=0;
			my $start=0;
			my $mention_region=0;
			my $mention_Regex_region=0;
			#while($article_tmp=~/^(rs[0-9]+)()([CISQMNPKDTFAGHLRWVEYX][\W\-\_][CISQMNPKDTFAGHLRWVEYX]\s.*)$/ || $article_tmp=~/^(.*?)([\W\-\_])(.*)$/)
			while($article_tmp=~/^([A-Z]+)(.*)$/ || $article_tmp=~/^([a-z]+)(fs.*)$/ || $article_tmp=~/^([a-z]+)(.*)$/ || $article_tmp=~/^([0-9]+)(.*)$/ || $article_tmp=~/^([\W\-\_])(.*)$/ )
			{
				my $pre=$1;
				my $post=$2;
				my $end=length($pre)+$start;
				
				if ($pre ne " ")
				{
					#for location
					if(exists $locationmap_hash{$start+1})
					{
						$mention_region=$start+1;
						my $state="I";
						if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
						print locationA $pmid."	".$pre."	".($start+1)."	".$end."	".$state."\n";
						#$Result_hash{$count_token}=$Begin_hash{$type_hash{$locationmap_hash{$mention_region}}};
						$Result_hash{$count_token}=$state;
					}
					elsif($mention_region!=0)
					{
						if($lastmap_hash{$mention_region}>$end)
						{
							my $state="I";
							if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
							print locationA $pmid."	".$pre."	".($start+1)."	".$end."	".$state."\n";
							#$Result_hash{$count_token}=$Inside_hash{$type_hash{$locationmap_hash{$mention_region}}};
							$Result_hash{$count_token}=$state;
						}
						elsif($lastmap_hash{$mention_region}==$end)
						{
							my $state="I";
							if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
							print locationA $pmid."	".$pre."	".($start+1)."	".$end."	".$state."\n";
							#$Result_hash{$count_token}=$End_hash{$type_hash{$locationmap_hash{$mention_region}}};
							$Result_hash{$count_token}=$state;
						}
						else
						{
							print locationA $pmid."	".$pre."	".($start+1)."	".$end."\n";
							#$Result_hash{$count_token-1}=$End_hash{$type_hash{$locationmap_hash{$mention_region}}};
							$Result_hash{$count_token}=$state;
							$mention_region=0;
						}
					}
					else
					{
						print locationA $pmid."	".$pre."	".($start+1)."	".$end."\n";
					}
					
					#for RegEx
					if(exists $RegEx_locationmap_hash{$start+1})
					{
						$mention_Regex_region=$start+1;
						$RegEx_result_hash{$count_token}=$Begin_hash{$RegEx_type_hash{$RegEx_locationmap_hash{$mention_Regex_region}}};
					}
					elsif($mention_Regex_region!=0)
					{
						if($RegEx_lastmap_hash{$mention_Regex_region}>=$end)
						{
							$RegEx_result_hash{$count_token}=$Inside_hash{$RegEx_type_hash{$RegEx_locationmap_hash{$mention_Regex_region}}};
						}
						else
						{
							$RegEx_result_hash{$count_token-1}=$End_hash{$RegEx_type_hash{$RegEx_locationmap_hash{$mention_Regex_region}}};
							$mention_Regex_region=0;
						}
					}
					$count_token++;
				}
				
				$article_tmp=$post;
				if ($char ne "")
				{
					$start=$end+1;
				}
				else
				{
					$start=$end;
				}
			}
			if($article ne "")
			{
				print locationA "\n";
			}
			close locationA;
			

			my $sentence=$article;
			$sentence=~s/([0-9])([A-Za-z])/$1 $2/g;
			$sentence=~s/([A-Z])([a-z])/$1 $2/g;
			$sentence=~s/([a-z])([A-Z])/$1 $2/g;
			$sentence=~s/([A-Za-z])([0-9])/$1 $2/g;
			$sentence=~s/([a-z])(fs)/$1 $2/g;
			$sentence=~s/([\W\-\_])/ $1 /g;
			$sentence =~ s/[ ]+/ /g;
			
			#tokens
			my @tokens = split(" ",$sentence);
			
			#Stemming
			my @stemmed_tokens = Text::English::stem( @tokens ); 
			
			#POS	
			my @POS_tokens=();
			#my $tagged_text = $p->add_tags( $sentence ); 
			#my $count_pos=0;
			#while($tagged_text =~/^<([a-z]+?)>(.+?)<\/(.+?)>(.+)$/)
			#{
			#	$postag=$1;
			#	$POS_tokens[$count_pos]=$postag;
			#	$count_pos++;
			#	$tagged_text=$4;
			#	$tagged_text =~ s/^ //g;
			#}
			
			open outputA,">>".$output;
			$count_token=0;
			my $last_token="";
			foreach $token(@tokens)
			{
				$insert_query="";
				if($stemmed_tokens[$count_token] eq "")
				{
					$stemmed_tokens[$count_token]=$token;
				}
				if($POS_tokens[$count_token] eq "")
				{
					$POS_tokens[$count_token]=$token;
				}
				$temmed=$stemmed_tokens[$count_token];
				$POStag=$POS_tokens[$count_token];
				
				#Number of Numbers [0-9] 
				my $tmp=$token;
				$tmp=~s/[^0-9]//g;
				my $Num_num="";
				if(length($tmp)>3){$Num_num="N:4+";}else{$Num_num="N:".length($tmp);}
				
				#Number of Uppercase [A-Z]
				$tmp=$token;
				$tmp=~s/[^A-Z]//g;
				my $Num_Uc="";
				if(length($tmp)>3){$Num_Uc="U:4+";}else{$Num_Uc="U:".length($tmp);}
				
				#Number of Lowercase [a-z]
				$tmp=$token;
				$tmp=~s/[^a-z]//g;
				my $Num_Lc="";
				if(length($tmp)>3){$Num_Lc="L:4+";}else{$Num_Lc="L:".length($tmp);}
				
				#Number of ALL char
				$tmp=$token;
				$tmp=~s/[^a-z]//g;
				my $Num_All="";
				if(length($tmp)>3){$Num_All="A:4+";}else{$Num_All="A:".length($tmp);}
				
				#specific character (;:,.->+_)
				$tmp=$token;
				my $SpecificC="";
				if($tmp=~/[\;\:\,\.\-\>\+\_]/){$SpecificC="-SpecificC1-";}
				elsif($tmp=~/[\(\)]/){$SpecificC="-SpecificC2-";}
				elsif($tmp=~/[\{\}]/){$SpecificC="-SpecificC3-";}
				elsif($tmp=~/[\[\]]/){$SpecificC="-SpecificC4-";}
				elsif($tmp=~/[\\\/]/){$SpecificC="-SpecificC5-";}
				else{$SpecificC="__nil__";}
				
				#chromosomal keytokens
				$tmp=$token;
				my $ChroKey="";
				if($tmp=~/^(q|p|q[0-9]+|p[0-9]+|qter|pter|XY|t)$/){$ChroKey="-ChroKey-";}else{$ChroKey="__nil__";}
				
				#Mutation type
				$tmp=lc($token);
				my $MutatType="";
				if($tmp=~/(del|ins|dup|tri|qua|con|delins|indel)/){$MutatType="-MutatType-";}else{$MutatType="__nil__";}
				if($tmp=~/(fs|fsX|fsx)/){$MutatType="-FrameShiftType-";}else{$MutatType="__nil__";}
				
				#Mutation word
				$tmp=lc($token);
				my $MutatWord="";
				if($tmp=~/^(deletion|delta|elta|insertion|repeat|inversion|deletions|insertions|repeats|inversions)$/){$MutatWord="-MutatWord-";}
				else{$MutatWord="__nil__";}
				
				#Mutation article & basepair
				$tmp=lc($token);
				my $MutatArticle="";
				if($tmp=~/^(single|a|one|two|three|four|five|six|seven|eight|nine|ten|[0-9]+|[0-9]+\.[0-9]+)$/){$MutatArticle="-Base-";}
				if($tmp=~/(kb|mb)/){$MutatArticle="-Byte-";}
				elsif($tmp=~/(base|bases|pair|amino|acid|acids|codon|postion|postions|bp|nucleotide|nucleotides)/){$MutatArticle="-bp-";}
				else{$MutatArticle="__nil__";}
				
				#Type1
				$tmp=lc($token);
				my $Type1="";
				if($tmp=~/^[cgrm]$/){$Type1="-Type1-";}
				elsif($tmp=~/^(ivs|ex|orf)$/){$Type1="-Type1_2-";}
				else{$Type1="__nil__";}
				
				#Type2
				$tmp=$token;
				my $Type2="";
				if($tmp eq "p"){$Type2="-Type2-";}else{$Type2="__nil__";}
				
				#DNA symbols
				$tmp=$token;
				my $DNASym="";
				if($tmp=~/^[ATCGUatcgu]$/){$DNASym="-DNASym-";}else{$DNASym="__nil__";}
				
				#Protein symbols
				my $uc_tmp=$token;
				my $lc_tmp=lc($token);
				my $ProteinSym="";
				if($lc_tmp=~/(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift)/){$ProteinSym="-ProteinSymFull-";}
				elsif($lc_tmp=~/^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)$/){$ProteinSym="-ProteinSymTri-";}
				elsif($lc_tmp=~/^(ys|le|er|ln|et|sn|ro|ys|sp|hr|he|la|ly|is|eu|rg|rp|al|lu|yr)$/ && $last_token=~/^[CISQMNPKDTFAGHLRWVEYX]$/){$ProteinSym="-ProteinSymTriSub-";}
				elsif($uc_tmp=~/^[CISQMNPKDTFAGHLRWVEYX]$/){$ProteinSym="-ProteinSymChar-";}
				else{$ProteinSym="__nil__";}
				
				#RS
				$tmp=$token;
				my $RScode="";
				if($tmp=~/^(rs|RS|Rs)[0-9]/){$RScode="-RScode-";}
				elsif($tmp=~/^(rs|RS|Rs)$/){$RScode="-RScode-";}
				else{$RScode="__nil__";}
				
				#Patterns
				my $Pattern1=$token;
				if($Pattern1=~/[\W\-\_]/){$Pattern1="__nil__";}
				else{
					$Pattern1=~s/[A-Z]/A/g;
					$Pattern1=~s/[a-z]/a/g;
					$Pattern1=~s/[0-9]/0/g;
					$Pattern1="P1:".$Pattern1;
				}
				my $Pattern2=$token;
				if($Pattern2=~/[\W\-\_]/){$Pattern2="__nil__";}
				else{
					$Pattern2=~s/[A-Za-z]/a/g;
					$Pattern2=~s/[0-9]/0/g;
					$Pattern2="P2:".$Pattern2;
				}
				my $Pattern3=$token;
				if($Pattern3=~/[\W\-\_]/){$Pattern3="__nil__";}
				else{
					$Pattern3=~s/[A-Z]+/A/g;
					$Pattern3=~s/[a-z]+/a/g;
					$Pattern3=~s/[0-9]+/0/g;
					$Pattern3="P3:".$Pattern3;
				}
				my $Pattern4=$token;
				if($Pattern4=~/[\W\-\_]/){$Pattern4="__nil__";}
				else{
					$Pattern4=~s/[A-Za-z]+/a/g;
					$Pattern4=~s/[0-9]+/0/g;
					$Pattern4="P4:".$Pattern4;
				}
				
				#prefix-pattern
				my $prefix="";
				$temp=$token;
				if(length($temp)>=1){ $prefix=$prefix.substr($temp,0,1);}else { $prefix=$prefix."__nil__";}
				if(length($temp)>=2){ $prefix=$prefix." ".substr($temp,0,2);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=3){ $prefix=$prefix." ".substr($temp,0,3);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=4){ $prefix=$prefix." ".substr($temp,0,4);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=5){ $prefix=$prefix." ".substr($temp,0,5);}else { $prefix=$prefix." __nil__";}
				
				#suffix-pattern
				my $suffix="";
				$temp=$token;
				if(length($temp)>=1){ $suffix=$suffix.substr($temp,-1,1);}else { $suffix=$suffix."__nil__";}
				if(length($temp)>=2){ $suffix=$suffix." ".substr($temp,-2,2);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=3){ $suffix=$suffix." ".substr($temp,-3,3);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=4){ $suffix=$suffix." ".substr($temp,-4,4);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=5){ $suffix=$suffix." ".substr($temp,-5,5);}else { $suffix=$suffix." __nil__";}		
				
				if((not exists $Result_hash{$count_token}) || $Result_hash{$count_token} eq ""){$Result_hash{$count_token}="O";}
				print outputA $token." ".$temmed." ".$POStag." ".$Num_num." ".$Num_Uc." ".$Num_Lc." ".$Num_All." ".$SpecificC." ".$ChroKey." ".$MutatType." ".$MutatWord." ".$MutatArticle." ".$Type1." ".$Type2." ".$DNASym." ".$ProteinSym." ".$RScode." ".$Pattern1." ".$Pattern2." ".$Pattern3." ".$Pattern4." ".$prefix." ".$suffix;
				if(not exists $RegEx_result_hash{$count_token}){$RegEx_result_hash{$count_token}="O";}
				print outputA " ".$RegEx_result_hash{$count_token}."\n"; #RegEx result
				$count_token++;
				$last_token=$token;
			}
			if($article ne "")
			{
				print outputA "                                  \n";
			}
			close outputA;
			
			%RegEx_locationmap_hash=();
			%RegEx_lastmap_hash=();
			%RegEx_mention_hash=();
			%RegEx_type_hash=();
			%RegEx_result_hash=();
			
			%Result_hash=();
			%locationmap_hash=();
			%lastmap_hash=();
			%mention_hash=();
			%type_hash=();
			%character_hash=();
			$article="";
			$mode=0;
		}
	}
	close input;
	return 1;
}

return 1;
