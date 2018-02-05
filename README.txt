
This program takes a directory full of nascent fastq files as inputs. It outputs bams, bedgraphs, and tdfs. 

to run:
python create_scipts_to_map_on_fiji_argparse.py <indir> <outdir> <genome> <email>

for example:
python create_scipts_to_map_on_fiji_argparse.py /scratch/Users/allenma/pro/fastqfiles/ /scratch/Users/allenma/pro/ hg19 Mary.A.Allen@colorado.edu

important: 
You must git clone this somewhere on /scratch/. Additionaly, Both your indir and outdir must be on scratch.


More information: 
This script will make a slurm script for processing GRO-seq data and submit it to the queue on
fiji. Buy defalut the program will map the reads using bowtie2, then convert
the sam to a sorted bam file, then create bedgraphs, then make tdf files which
can be used in the program IGV.

positional arguments:
  indir                 Directory with fastq or bam files. Must be located on
                        /scratch/.
  outdir                Directory files will be output into. The program will
                        create this direcoty if it does not exist. This
                        direcory must be located on /scratch/.
  genome                Which genome do you wish to map too? Allowed values
                        are hg19, dm3, mm10,
  email                 Email address is required. If you mistype your email
                        then IT gets your emails. Please don't do that.

optional arguments:
  -h, --help            show this help message and exit
  -f {True,False}, --flipreads {True,False}
  -q {True,False}, --checkquality {True,False}
  -t {True,False}, --trimgalore {True,False}
  -m {True,False}, --fastqtosam {True,False}
  -s2b {True,False}, --samtobam {True,False}
  -b2sb {True,False}, --bamtosortedbam {True,False}
  -ib {True,False}, --sortedbamtobai {True,False}
  -bg {True,False}, --sortedbamtobedgraph {True,False}
  -c {True,False}, --readcountcorrectbedgraph {True,False}
  -igv {True,False}, --igvcreate {True,False}
