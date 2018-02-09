import argparse
import sys
import os
import glob 
#this is built for mapping single end nascent data!!! This program expects a indir and a outdir. They must both end in a slash.



#First choose which steps do you want to run

inputfiletype=".fastq"
#Flip reads is required for some protocols such as PRO-seq and the NEB random prime kit
#flipreads = True
#Check read quality?
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
#submit=True


#What bowtie options do you want to use
bowtieoptions = "--very-sensitive"
processers="1" #this will move to 32 if bowtie is being used
queue="long"
memamount="200"
memtype="gb" #gb, mb, kb
time="48:00:00" #hours:min:seconds
igvdir="/opt/igvtools/2.3.75/"

#-------- these are variables you probably don't need to check
genomesdir="/scratch/Shares/public/genomes/"
def createbowtieindexdic():
	b = {}
	bowtiefiles = [infile for infile in glob.glob(os.path.join(genomesdir, '*/*/*/Sequence/Bowtie2Index/genome*'))]
	for bowtiefile in bowtiefiles:
		dname = os.path.dirname(bowtiefile)
		rootname=bowtiefile.split("/")[7]
		b[rootname]=dname+"/genome"	
	return b
def creategenomefadic():
	g = {}
	genomefastafiles = [infile for infile in glob.glob(os.path.join(genomesdir, '*/*/*/Sequence/WholeGenomeFasta/genome.fa'))]
	for genomefa in genomefastafiles:
		rootname=genomefa.split("/")[7]
		g[rootname] = genomefa

	return g
genomefadic = creategenomefadic()
bowtieindexdic = createbowtieindexdic()

memdic = {"gb":"G", "kb":"K", "mb":"M"}
samtoolsmem = memamount+memdic[memtype]
sbatchmem = memamount+memtype
bedtoolsgenomedir="/scratch/Shares/public/genomes/bedtools_genome_files/"
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


def bowtie2file(passfile, outdir, wfname, genome):
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

def sortedbamtobedgraphfile(sortedbam, outdir, wfname, genome):
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("module load bedtools/2.25.0\n")
	wf.write("echo bedgraph\n")
	rootname = sortedbam.split("/")[-1]
	rootname = rootname.strip(".sorted.bam")
	Bedgraphoutdir = outdir+"/bedgraphs/"
	ensure_dir(Bedgraphoutdir)
	bedtoolsgenomefile=bedtoolsgenomedir+genome+".chrom.sizes.genome"
	wf.write("genomeCoverageBed -bg -strand + -ibam "+sortedbam+" -g "+bedtoolsgenomefile+" > "+Bedgraphoutdir+rootname+".pos.BedGraph\n")
	wf.write("genomeCoverageBed -bg -strand - -ibam "+sortedbam+" -g "+bedtoolsgenomefile+" | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' > "+Bedgraphoutdir+rootname+".neg.BedGraph\n")
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
                        

def igvcreatefile(readcountcorrBedgraphfile, outdir, wfname, genome):
	wf = open(wfname, "a")
        wf.write("date\n")
	tdffiledir=outdir+"tdfs/"
	ensure_dir(tdffiledir)
	rootname = readcountcorrBedgraphfile.split("/")[-1]
	rootname = rootname.strip(".mp.BedGraph")
	tdffile = tdffiledir+rootname+".tdf"
	#igvgenomefile = igvdir+"/genomes/"+genome+".chrom.sizes" 
	#igvgenomefile=bedtoolsgenomedir+genome+".chrom.sizes.genome"
	igvgenomefile = genomefadic[genome]
	wf.write(igvdir+"/igvtools toTDF "+readcountcorrBedgraphfile+" "+tdffile+" "+igvgenomefile+"\n")
        wf.write("echo tdf\n")
        wf.write("date\n")



def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
	

def main(indir, outdir, email, genome):
	if args.Turnoff_fastqtosam==False:
        	processers="32"
	#make outdir
	ensure_dir(outdir)
	#make e_and_o
	os.system("cp --no-clobber readcountcorrectBG.py "+outdir+"\n")
	ensure_dir(outdir+"e_and_o/")
	ensure_dir(outdir+"qsubscripts/")
	indirfiles=[infile for infile in glob.glob(os.path.join(indir, '*'+inputfiletype))]
	for filename in indirfiles:
		passfile = filename
		rootname = filename.split("/")[-1]
		rootname = rootname.strip(inputfiletype)
		wfname = outdir+"qsubscripts/"+rootname+"_map.slurm"
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
		if args.flipreads==True:
			passfile = flipfile(passfile, outdir, wfname) 
		if args.Turnoff_checkquality==False:
			checkqualityfile(passfile, outdir, wfname)
		if args.Turnoff_trimgalore==False:
			passfile = trimgalorefile(passfile, outdir, wfname)
			checkqualityfile(passfile, outdir, wfname)
		if args.Turnoff_fastqtosam==False:
			passfile = bowtie2file(passfile,  outdir, wfname, args.genome)
		if args.Turnoff_samtobam==False:
			passfile = samtobamfile(passfile, outdir, wfname)
		if args.Turnoff_bamtosortedbam==False:
			sortedbam = bamtosortedbamfile(passfile, outdir, wfname)
		if args.Turnoff_sortedbamtobai==False:
			sortedbamtobai(sortedbam, wfname)
		if args.Turnoff_sortedbamtobedgraph==False:
			Bedgraphfile = sortedbamtobedgraphfile(sortedbam, outdir, wfname, genome)
		if args.Turnoff_readcountcorrectbedgraph==False:
			readcountcorrBedgraphfile = readcountcorrectbedgraphfile(Bedgraphfile, sortedbam, outdir, wfname)
		if args.Turnoff_igvcreatetdf==False:
			igvcreatefile(readcountcorrBedgraphfile, outdir, wfname, genome)
		wf = open(wfname, "a")
		wf.write("date\n")
		wf.close()
		if args.Turnoff_submit==False:	
			os.system('sbatch '+wfname)	
	
	if args.Turnoff_submit==False:
		print "To check if your scripts are still runing, type squeue or qstat on the command line."
		print "If your scripts error the scripts are in "+outdir+"qsubscripts/ and the error files are in "+outdir+"e_and_o/"
		print "Once your jobs have ended check the files in "+outdir+"qual/ for quality reports."
		print "Make sure you rsync what you want backed up to /projects/."
	else:
		print "You now need to check your shell scripts in "+outdir+"qsubscripts/"
		print "Once you decide you like them submit using sbatch"
		




if __name__=="__main__":
	parser = argparse.ArgumentParser(description='This program will make a slurm script for processing GRO-seq data and submit it to the queue on fiji.\n Buy defalut the program will map the reads using bowtie2, then convert the sam to a sorted bam file, then create bedgraphs, then make tdf files which can be used in the program IGV.',usage='python create_scipts_to_map_on_fiji_argparse.py <indir> <outdir> <genome> <email>')
	parser.add_argument("-f","--flipreads", action='store_true', dest="flipreads", default=False)
	parser.add_argument("-q","--Turnoff_checkquality", action='store_true', default=False)
	parser.add_argument("-t","--Turnoff_trimgalore", action='store_true', default=False)
	parser.add_argument("-m","--Turnoff_fastqtosam", action='store_true', default=False)
	parser.add_argument("-s2b","--Turnoff_samtobam", action='store_true', default=False)
	parser.add_argument("-b2sb","--Turnoff_bamtosortedbam", action='store_true', default=False)
	parser.add_argument("-ib","--Turnoff_sortedbamtobai", action='store_true', default=False)
	parser.add_argument("-bg","--Turnoff_sortedbamtobedgraph", action='store_true', default=False)
	parser.add_argument("-c","--Turnoff_readcountcorrectbedgraph", action='store_true', default=False)
	parser.add_argument("-igv","--Turnoff_igvcreatetdf", action='store_true', default=False)
	parser.add_argument("-s","--Turnoff_submit", action='store_true', default=False)
	parser.add_argument("indir", help="Directory with fastq or bam files. Must be located on /scratch/.", type=str)
	parser.add_argument("outdir", help="Directory files will be output into. The program will create this direcoty if it does not exist. This direcory must be located on /scratch/.", type=str)
	genomeschoices = bowtieindexdic.keys()
	parser.add_argument("genome", metavar="genome",help="Which genome do you wish to map too? Allowed values are "+', '.join(genomeschoices)+',', type=str, choices=genomeschoices)
	parser.add_argument("email", help="Email address is required. If you mistype your email then IT gets your emails. Please don't do that. ", type=str)
	args = parser.parse_args()
	print args
##this is built for mapping single end nascent data!!! This program expects a indir and a outdir. They must both end in a slash.
	if args.indir.startswith("/scratch/") and args.outdir.startswith("/scratch/"):
			ensure_dir(args.outdir)
                        main(args.indir, args.outdir, args.email, args.genome)
	else:
		print "indir and outdir must be in /scratch/"


 
#test
#python create_scipts_to_map_on_fiji_argparse.py /scratch/Users/allenma/171025_NB501447_0179_fastq/Allen_Dowell-371/ /scratch/Users/allenma/pro/ hg19 allenma@colorado.edu 
