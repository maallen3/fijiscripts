import argparse
import sys
import os
import glob 
#this is built for mapping single end nascent data!!! This program expects a indir and a outdir. They must both end in a slash.



#First choose which steps do you want to run

inputfiletype=".fastq"
#Flip reads is required for some protocols such as PRO-seq and the NEB random prime kit
flipreads = True
#Check read quality?
checkquality=True
#Trim adaptors?
trimgalore = True
fastqtosam = True
samtobam = True
bamtosortedbam = True
sortedbamtobai = True
sortedbamtobedgraph = True
readcountcorrectbedgraph = True
igvcreate = True
millionsmappedcount = True

submit=True
#What bowtie options do you want to use
bowtieoptions = "--very-sensitive"
processers="1"
if fastqtosam==True:
	processers="32"
queue="long"
memamount="200"
memtype="gb" #gb, mb, kb
time="48:00:00" #hours:min:seconds
igvdir="/opt/igvtools/2.3.75/"

#-------- these are variables you probably don't need to check
bowtieindexdic = {"hg19":"/scratch/Shares/dowell/pubgro/genomefiles/bowtiebwaindexs/hg19_Bowtie2_index", "mm10":"/scratch/Shares/dowell/pubgro/genomefiles/bowtiebwaindexs/mm10_Bowtie2_index", "dm3":"/scratch/Shares/dowell/pubgro/genomefiles/bowtiebwaindexs/dm3.fa.Bowtie2"}
memdic = {"gb":"G", "kb":"K", "mb":"M"}
samtoolsmem = memamount+memdic[memtype]
sbatchmem = memamount+memtype

#--------------all scripts below this line

def flipfile(passfile, outdir, wfname):
	wf = open(wfname, "a")
	rootname = passfile.split("/")[-1]
        rootname = rootname.strip(".fastq")
	outfastqdir = outdir+"flipped/"
        ensure_dir(outfastqdir)
	outfastq = outfastqdir+rootname+".flip.fastq"
	wf.write("date\n")
	wf.write("\n")
	wf.write("module load fastx-toolkit/0.0.13\n")
	wf.write("/opt/fastx-toolkit/0.0.13/bin/fastx_reverse_complement -Q33 -i "+passfile+" -o "+outfastq+" \n")	 
	wf.write("echo flipped\n")
	wf.write("date\n")
	wf.close()
	return outfastq

def checkqualityfile(passfile, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("module load fastqc/0.11.5\n")
	qualoutdir = outdir+"qual/"
	ensure_dir(qualoutdir)
	wf.write("fastqc "+passfile+" -o "+qualoutdir+"\n")
	wf.write("echo qual\n")
	wf.write("date\n")
        wf.close()

def trimgalorefile(passfile, outdir, wfname):
	rootname = passfile.split("/")[-1]
	rootname = rootname.strip(".fastq")
	cutoutdir = outdir+"cutadapt/"
	ensure_dir(cutoutdir)
	#MA_DMSO_Groseq_S7_R1_001.flip_trimmed.fq
	#MA_DMSO_Groseq_S7_R1_001.flip.fastq  
	trimmedfastq = cutoutdir+rootname+"_trimmed.fq"
	trimmedfastq2 = cutoutdir+rootname+".trimmed.fastq"
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("\n")
        wf.write("module load python/2.7.14/cutadapt/1.12\n")
	wf.write("module load trim_galore/0.4.3\n")
        wf.write("echo $PATH\n")
        wf.write("echo $PYTHONPATH\n")
        wf.write("/opt/trim_galore/0.4.3/trim_galore --path_to_cutadapt /opt/cutadapt/python-2.7.14/1.12/bin/cutadapt -o "+cutoutdir+" "+passfile+"\n")
	wf.write("mv "+trimmedfastq+" "+trimmedfastq2+"\n")
	wf.write("echo trimmed\n")
        wf.write("date\n\n\n")
	wf.close()
        return trimmedfastq2


def bowtie2file(passfile, outdir, wfname):
        wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("module load samtools/1.3.1\n")
	wf.write("module load bowtie/2.2.9\n\n")
	samdir = outdir+"sams/"
	ensure_dir(samdir)
	rootname = passfile.split("/")[-1]
	rootname = rootname.strip(".fastq")
	samfile = samdir+rootname+".sam"
	stderrfile = samdir+rootname+".stderr"
	bowtieindex = bowtieindexdic[genome]
	wf.write("bowtie2 -p "+processers+" "+bowtieoptions+" -x "+bowtieindex+" -U "+passfile+" >"+samfile+" 2>"+stderrfile+" \n")
	wf.write("echo mapped\n")
	wf.write("date\n")
	wf.close()
        return samfile


def samtobamfile(passfile, outdir, wfname):
        wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("wc -l "+passfile+" > "+passfile+".wc\n")
	rootname=passfile.split("/")[-1]
	rootname=rootname.strip(".sam")
	bamoutdir = outdir+"bams/"
	ensure_dir(bamoutdir)
	bamfile = bamoutdir+rootname+".bam"
	wf.write("samtools view -S -b -o "+bamfile+" "+passfile+" 2>"+bamfile+".err\n")
	wf.write("samtools flagstat "+bamfile+" > "+bamfile+".flagstat 2>"+bamfile+".flagstat.err\n")
	wf.write("echo bam\n")
	wf.write("date\n")
        wf.close()
        return bamfile


def bamtosortedbamfile(bamfile, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	sortedbamfile = bamfile.strip(".bam")+".sorted.bam"
	sortedbamdir = outdir+"sortedbams/"
	ensure_dir(sortedbamdir)
	sortedbamfile = sortedbamdir+sortedbamfile.split("/")[-1]
	wf.write("samtools sort -m "+samtoolsmem+" "+bamfile+" >"+sortedbamfile+"\n")
	wf.write("samtools flagstat "+sortedbamfile+" >"+sortedbamfile+".flagstat 2>"+sortedbamfile+".flagstat.err\n")
	wf.write("echo sorted.bam\n")
        wf.write("date\n")
	return sortedbamfile
                        


def sortedbamtobai(sortedbam, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("samtools index "+sortedbam+"\n")
        wf.write("echo indexed.bam\n")
        wf.write("date\n")

def sortedbamtobedgraphfile(sortedbam, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("module load bedtools/2.25.0\n")
	wf.write("echo bedgraph\n")
	rootname = sortedbam.split("/")[-1]
	rootname = rootname.strip(".sorted.bam")
	Bedgraphoutdir = outdir+"/bedgraphs/"
	ensure_dir(Bedgraphoutdir)
	wf.write("genomeCoverageBed -bg -strand + -ibam "+sortedbam+" -g "+genome+" > "+Bedgraphoutdir+rootname+".pos.BedGraph\n")
	wf.write("genomeCoverageBed -bg -strand - -ibam "+sortedbam+" -g "+genome+" | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' > "+Bedgraphoutdir+rootname+".neg.BedGraph\n")
	wf.write("cat "+Bedgraphoutdir+rootname+".pos.BedGraph "+Bedgraphoutdir+rootname+".neg.BedGraph > "+Bedgraphoutdir+rootname+".unsorted.BedGraph\n")
	wf.write("sortBed -i "+Bedgraphoutdir+rootname+".unsorted.BedGraph >"+Bedgraphoutdir+rootname+".BedGraph\n")
        wf.write("date\n")
	Bedgraphfile = Bedgraphoutdir+rootname+".BedGraph"
	return Bedgraphfile

                        



def readcountcorrectbedgraphfile(Bedgraphfile, sortedbam, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("echo readcountcorrectedbedgraph\n")
	wf.write("module load python/2.7.14\n")
	#need to copy readcountcorrectBG.py to out dir
	rootname = Bedgraphfile.strip(".BedGraph")
	flagstatfile = sortedbam+".flagstat"
	readcountcorrBedgraphfile = rootname +".mp.BedGraph" 
	wf.write("python "+outdir+"readcountcorrectBG.py "+Bedgraphfile+" "+flagstatfile+" "+readcountcorrBedgraphfile+"\n")	
        wf.write("date\n")
	return readcountcorrBedgraphfile
                        

def igvcreatefile(readcountcorrBedgraphfile, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	tdffiledir=outdir+"tdfs/"
	ensure_dir(tdffiledir)
	rootname = readcountcorrBedgraphfile.split("/")[-1]
	rootname = rootname.strip(".mp.BedGraph")
	tdffile = tdffiledir+rootname+".tdf"
	igvgenomefile = igvdir+"/genomes/"+genome+".chrom.sizes" 
	wf.write(igvdir+"/igvtools toTDF "+readcountcorrBedgraphfile+" "+tdffile+" "+igvgenomefile+"\n")
        wf.write("echo tdf\n")
        wf.write("date\n")



def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
	

def main(indir, outdir):
	#make outdir
	ensure_dir(outdir)
	#make e_and_o
	os.system("cp --no-clobber readcountcorrectBG.py "+outdir+"\n")
	ensure_dir(outdir+"e_and_o/")
	indirfiles=[infile for infile in glob.glob(os.path.join(indir, '*'+inputfiletype))]
	wf2 = open(outdir+"rsync_out.sh", "w")
	for filename in indirfiles:
		passfile = filename
		rootname = filename.split("/")[-1]
		rootname = rootname.strip(inputfiletype)
		wfname = outdir+rootname+"_map.slurm"
		wf = open(wfname, "w")
		wf.write("#!/bin/bash\n")
		wf.write("#SBATCH --job-name="+rootname+"  # Job name\n")
		wf.write("#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)\n")
		wf.write("#SBATCH --mail-user="+email+"# Where to send mail\n")
		wf.write("#SBATCH --nodes=1\n")
		wf.write("#SBATCH --ntasks="+processers+"# Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.\n")
		wf.write("#SBATCH --time="+time+" # Time limit hrs:min:sec\n")
		wf.write("#SBATCH -p "+queue+"\n")
		wf.write("#SBATCH --mem="+sbatchmem+" # Memory limit\n")
		wf.write("#SBATCH --output="+outdir+"e_and_o/"+rootname+".%j.out\n")
		wf.write("#SBATCH --error="+outdir+"e_and_o/"+rootname+".%j.err\n") 
		wf.write("\n\n\npwd; hostname; date\n")
		wf.close()
		if flipreads==True:
			passfile = flipfile(passfile, outdir, wfname) 
		if checkquality==True:
			checkqualityfile(passfile, outdir, wfname)
		if trimgalore==True:
			passfile = trimgalorefile(passfile, outdir, wfname)
			checkqualityfile(passfile, outdir, wfname)
		if fastqtosam==True:
			passfile = bowtie2file(passfile,  outdir, wfname)
		if samtobam==True:
			passfile = samtobamfile(passfile, outdir, wfname)
		if bamtosortedbam==True:
			sortedbam = bamtosortedbamfile(passfile, outdir, wfname)
		if sortedbamtobai==True:
			sortedbamtobai(sortedbam, wfname)
		if sortedbamtobedgraph==True:
			Bedgraphfile = sortedbamtobedgraphfile(sortedbam, outdir, wfname)
		if readcountcorrectbedgraph==True:
			readcountcorrBedgraphfile = readcountcorrectbedgraphfile(Bedgraphfile, sortedbam, outdir, wfname)
		if igvcreate==True:
			igvcreatefile(readcountcorrBedgraphfile, outdir, wfname)
		wf = open(wfname, "a")
		wf.write("date\n")
		wf.close()
		if submit==True:	
			os.system('sbatch '+wfname)	
	
	if submit==True:
		print "To check if your scripts are still runing, type squeue or qstat on the command line."
		print "If your scripts error they are in "+outdir
		print "Once your jobs have ended check these files for quality reports in "+outdir+"qual/."
		print "Once jobs are complete you can use the files in rsync_out.sh to move the imporant files to projects."
	else:
		print "instructions"
	wf2.close()




if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Make a slurm script for processing GRO-seq data and submit it to the queue on fiji.\n Buy defalut the program will map the reads using bowtie2, then convert the sam to a sorted bam file, then create bedgraphs, then make tdf files which can be used in the program IGV.',usage='python create_scipts_to_map_on_fiji_argparse.py indir outdir genome email')
	parser.add_argument("-f","--flipreads", choices=[True, False], default=False)
	parser.add_argument("-q","--checkquality", choices=[True, False], default=True)
	parser.add_argument("-t","--trimgalore", choices=[True, False], default=False)
	parser.add_argument("-m","--fastqtosam", choices=[True, False], default=True)
	parser.add_argument("-s2b","--samtobam", choices=[True, False], default=True)
	parser.add_argument("-b2sb","--bamtosortedbam", choices=[True, False], default=True)
	parser.add_argument("-ib","--sortedbamtobai", choices=[True, False], default=True)
	parser.add_argument("-bg","--sortedbamtobedgraph", choices=[True, False], default=True)
	parser.add_argument("-c","--readcountcorrectbedgraph", choices=[True, False], default=True)
	parser.add_argument("-igv","--igvcreate", choices=[True, False], default=True)
	parser.add_argument("indir", help="Directory with fastq or bam files. Must be located on /scratch/.", type=str)
	parser.add_argument("outdir", help="Directory files will be output into. The program will create this direcoty if it does not exist. This direcory must be located on /scratch/.", type=str)
	genomeschoices = bowtieindexdic.keys()
	parser.add_argument("genome", metavar="genome",help="Which genome do you wish to map too? Allowed values are "+', '.join(genomeschoices)+',', type=str, choices=genomeschoices)
	parser.add_argument("email", help="Email address is required. If you mistype your email then IT gets your emails. Please don't do that. ", type=str)
	args = parser.parse_args()
	print args
##this is built for mapping single end nascent data!!! This program expects a indir and a outdir. They must both end in a slash.



#First choose which steps do you want to run

#inputfiletype=".fastq"
#Flip reads is required for some protocols such as PRO-seq and the NEB random prime kit
#flipreads = True
##Check read quality?
#checkquality=True
#Trim adaptors?
#trimgalore = True
#fastqtosam = True
#samtobam = True
#bamtosortedbam = True
#sortedbamtobai = True
#sortedbamtobedgraph = True
#readcountcorrectbedgraph = True
#igvcreate = True
#millionsmappedcount = True
#	if len(sys.argv)<2:
#		print "to run"
#                print "python create_scipts_to_map_on_fiji.py <indir> <outdir> <genome> <your_email>"
#                print "\n\n"
#		print "genome options are hg19, mm10, dm3"
#		print "Remember that the indir(in direcory) and outdir(out directory) must be on /scratch/!!!"
#		print "But you should ALWAYs keep a copy of raw data on /projects/ becuase /projects/ is backed up and scratch is not."
#	else:
#	        indir = sys.argv[1]
#	        outdir = sys.argv[2]
#        	genome = sys.argv[3]
#        	email=sys.argv[4]
#		if indir.startswith("/scratch/") and outdir.startswith("/scratch/"):
#
#			main(indir, outdir)
#		else:
#			print "Your in directory and out directory must be on /scratch/", indir, outdir


 
#test
#python create_scipts_to_map_on_fiji.py /scratch/Users/allenma/171025_NB501447_0179_fastq/Allen_Dowell-371/ /scratch/Users/allenma/pro/ hg19 allenma@colorado.edu 
