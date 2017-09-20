[Directory]
A. Introduction of folders
B. Installation
C. Instruction
D. FULL Usage
E. Retrieve articles via E-utilities (For PubMed abstracts & PMC fulltext articles)

#======================================================================================#

A. [Introduction of folders]
	
	Library & Module:
		[lib]
		[CRF++]
	Model & Regular Expression tables:
		[model]
		[RegEx]
	Input data folder:
		[input]
	Output data folder:
		[output]
	TMP folder:
		[tmp]
	
B. [Installation] 

	*Users don't need to install CRF++ module if their operation system is windows.
	
	*Installation
		*Windows Environment
			Users need to install Perl. User can download ActivePerl via this link: http://www.activestate.com/activeperl/downloads
		
		*Linux Environment
			Users need to execute "Installation.sh" to install crf_test first.
	
			$ sh Installation.sh
					
C. [Instruction]

	$ perl tmVar.pl -i [Input] -o [Output] -s [Setup]
	
	Example: perl tmVar.pl -i input -o output -s setup.txt
	
	Input: User can provide the input data folder route ("input" or "C:\input"). 
	Output: User can provide the output data folder route ("output" or "C:\output"). 
	Setup: User can choose the mutation types for display. Default is "setup.txt". ("setup.txt" or "C:\setup.txt"). 


D. [FULL Usage] 

	tmVar is developed to find mutation mentions and individual components. 

	INPUT:
	
		Input file folder. Each input file should follow the PubTator format below or BioC format(http://bioc.sourceforge.net/). 

	RUN: 

		perl tmVar.pl -i input -o output

	OUTPUT:

		Output file folder. Each output file should follow the PubTator format below or BioC format(http://bioc.sourceforge.net/). 
				
	Input File Format:
	
		<PMID>|t|<TITLE>
		<PMID>|a|<ABSTRACT>
		...
				
		Example: 
		22016685|t|A novel missense mutation Asp506Gly in Exon 13 of the F11 gene in an asymptomatic Korean woman with mild factor XI deficiency.
		22016685|a|Factor XI (FXI) deficiency is a rare autosomal recessive coagulation disorder most commonly found in Ashkenazi and Iraqi Jews, but it is also found in other ethnic groups. It is a trauma or surgery-related bleeding disorder, but spontaneous bleeding is rarely seen. The clinical manifestation of bleeding in FXI deficiency cases is variable and seems to poorly correlate with plasma FXI levels. The molecular pathology of FXI deficiency is mutation in the F11 gene on the chromosome band 4q35. We report a novel mutation of the F11 gene in an 18-year-old asymptomatic Korean woman with mild FXI deficiency. Pre-operative laboratory screen tests for lipoma on her back revealed slightly prolonged activated partial thromboplastin time (45.2 sec; reference range, 23.2-39.4 sec). Her FXI activity (35%) was slightly lower than the normal FXI activity (reference range, 50-150%). Direct sequence analysis of the F11 gene revealed a heterozygous A to G substitution in nucleotide 1517 (c.1517A>G) of exon 13, resulting in the substitution of aspartic acid with glycine in codon 506 (p.Asp506Gly). To the best of our knowledge, the Asp506Gly is a novel missense mutation, and this is the first genetically confirmed case of mild FXI deficiency in Korea.
		
	Output File Format:

		<PMID>|t|<TITLE>
		<PMID>|a|<ABSTRACT>
		<PMID><tab><OFFSET_START><tab><OFFSET_END><tab><MENTION><tab><TYPE><tab><Normalized Form>
		...
				
		Example: 
		22016685|t|A novel missense mutation Asp506Gly in Exon 13 of the F11 gene in an asymptomatic Korean woman with mild factor XI deficiency.
		22016685|a|Factor XI (FXI) deficiency is a rare autosomal recessive coagulation disorder most commonly found in Ashkenazi and Iraqi Jews, but it is also found in other ethnic groups. It is a trauma or surgery-related bleeding disorder, but spontaneous bleeding is rarely seen. The clinical manifestation of bleeding in FXI deficiency cases is variable and seems to poorly correlate with plasma FXI levels. The molecular pathology of FXI deficiency is mutation in the F11 gene on the chromosome band 4q35. We report a novel mutation of the F11 gene in an 18-year-old asymptomatic Korean woman with mild FXI deficiency. Pre-operative laboratory screen tests for lipoma on her back revealed slightly prolonged activated partial thromboplastin time (45.2 sec; reference range, 23.2-39.4 sec). Her FXI activity (35%) was slightly lower than the normal FXI activity (reference range, 50-150%). Direct sequence analysis of the F11 gene revealed a heterozygous A to G substitution in nucleotide 1517 (c.1517A>G) of exon 13, resulting in the substitution of aspartic acid with glycine in codon 506 (p.Asp506Gly). To the best of our knowledge, the Asp506Gly is a novel missense mutation, and this is the first genetically confirmed case of mild FXI deficiency in Korea.
		22016685	26	35	Asp506Gly	ProteinMutation	p|D|506|G
		22016685	1109	1118	c.1517A>G	DNAMutation	c|A|1517|G
		22016685	1206	1217	p.Asp506Gly	ProteinMutation	p|D|506|G
		22016685	1254	1263	Asp506Gly	ProteinMutation	p|D|506|G

E. [Retrieve articles]
	
	$ perl PreProcessing.pl -t [type] -i [input] -o [output]
	
	Example: perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt
	
	type: User can choose to one of the input type : PMID | PMCID | PMIDlist | PMCIDlist (PMIDlist and PMCIDlist are files which store the pmid or pmcid list.)
	Input: User can provide the input parameter, such as 22016685(for PMID type), PMIDlist.txt(for PMIDlist type). 
	Input: User can provide the output file name (input.txt). 
	
