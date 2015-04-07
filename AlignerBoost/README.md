AlignerBoost is a generalized software toolkit for boosting Next-Gen sequencing mapping precision
using a Bayesian based mapping quality framework.
================================================================================================== 

AlignerBoost works with any NGS aligners that can produce standard SAM/BAM alignment outputs,
i.e. Bowtie, Bowtie2, BWA, SOAP2, YOABS, YA, and also supports RNA-seq aligners (i.e. Tophat/Tophat2, STAR). 

AlignerBoost works by tuning NGS aligners to report all potential alignments,
then utilizes a Bayesian-based framework to accurately estimate the mapping quality
of ambiguously mapped reads.

AlignerBoost can dramatically increase mapping precision without a significant loss of
sensitivity under various experimental strategies.

AlignerBoost is SNP-aware, and higher quality alignments can be achieved if provided with known SNPs.
-----------------------------------------------------------------------------------------------------

-	Download and installation
	You can download the latest executable release from GitHub at: https://github.com/Grice-Lab/AlignerBoost/releases
	You can also download or fork and pull the source codes from GitHub at: https://github.com/Grice-Lab/AlignerBoost
	AlignerBoost is pure Java based, and you can run it without the need for installation on Unix/Linux, Mac OS X,
	and Windows by simply type "java -jar AlignerBoost.jar" in the shell/terminal.

-	Dependencies
