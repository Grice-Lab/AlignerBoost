AlignerBoost manual
===================

AlignerBoost is a generalized software toolkit for boosting Next-Gen sequencing
mapping precision using a Bayesian based mapping quality framework.

AlignerBoost works with any NGS aligners that can produce standard SAM/BAM alignment outputs,
i.e. Bowtie, Bowtie2, BWA, SOAP2, YOABS, YA, and also supports RNA-seq aligners (i.e. Tophat/Tophat2, STAR). 

AlignerBoost works by tuning NGS aligners to report all potential alignments,
then utilizes a Bayesian-based framework to accurately estimate the mapping quality
of ambiguously mapped reads.

AlignerBoost can dramatically increase mapping precision without a significant loss of
sensitivity under various experimental strategies.

AlignerBoost is SNP-aware, and higher quality alignments can be achieved if provided with known SNPs.

Download and installation
-------------------------
You can download the latest executable release from GitHub at: https://github.com/Grice-Lab/AlignerBoost/releases.
You can also download or fork and pull the source codes from GitHub at: https://github.com/Grice-Lab/AlignerBoost.
AlignerBoost is pure Java based, and you can run it without the need for installation on Unix/Linux, Mac OS X, and Windows by simply type "java -jar AlignerBoost.jar" in the shell/terminal.

Dependencies
------------
AlignerBoost does not dependent on any 3rd party library directly. However, if you are using AlignerBoost's
best practice to generate executable shell scripts, you do need to have your NGS aligner
of choice available in the PATH to be able to run these scripts. You might also need other programs in PATH
for some other AlignerBoost pre-processing functionality. See "examples/README.example" for best practice.

Customized SAM format tags
--------------------------
AlignerBoost uses a set of customized tags in generated SAM/BAM files to store auxiliary alignment information
calculated during its filter process. These tags are listed below.
Note: X?: global tags, Y? seed region reated tags, Z? entire alignment related tags
  Tag  Type  Description
  XA   i     alignment length, including M,=,X,I,D,S but not H,P,N
  XL   i     insert length, including M,=,X,I,D but not S,H,P,N, determined by Cigar or 1DP
  XF   i     actual insert from (start) relative to reference
  XI   f     alignment identity as 1 - (YX + YG) / XL
  XH   Z     alignment likelihood given this mapping locus and base quality, in string format to preserve double precision
  XV   i     known SNVs (if any) used in calculating XH
  XP   Z     alignment posterior probability in string format to preserve double precision
  YL   i     seed length
  YX   i     No. of seed mismatches
  YG   i     No. of seed indels
  ZX   i     No. of all mismatches
  ZG   i     No. of all indels