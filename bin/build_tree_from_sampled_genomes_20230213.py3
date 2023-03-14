#!/usr/bin/python3
#-*- coding:utf8 -*-


import csv
import re
import os
import logging
import time
import argparse
import multiprocessing as mp
from collections import *
import pandas as pd
from Bio import SeqIO
import subprocess
import yaml
from io import StringIO
import glob
import warnings
import random

warnings.filterwarnings("ignore")

clearline = '\x1b[1A\x1b[2K'

parser = argparse.ArgumentParser(description="Acinetobacter project",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c","--configfile",help="Configfile",default="build_tree_from_sampled_genomes_20230213.yaml")
parser.add_argument("-n",help="CPU number to use",type=int)
pargs = parser.parse_args()
lock = mp.Lock()

def msg(txt,clearline=False,error=False):
	if error:
		logging.error(txt)
	else:
		logging.info(txt)
	if error:
		txt = "\x1b[31m" + time.strftime("%H:%M:%S ") + " ERROR " + txt + "\x1b[0m"
	else:
		txt = time.strftime("%H:%M:%S ") + txt
	if clearline:
		print('\x1b[1A\x1b[2K' + txt)
	else:
		print(txt)

def csvrdr(f):
	return csv.reader(f,delimiter="\t")

def sysexec(cmd):
	logging.info(cmd)
	print(cmd)
	os.system(cmd)

def read_cmd(cmd,nostrip=False):
	if nostrip:
		return subprocess.check_output(cmd,shell=True,universal_newlines=True)
	else:
		return subprocess.check_output(cmd,shell=True,universal_newlines=True).strip()

def mkdir_force(d):
	if os.path.isdir(d):
		pass
	else:
		os.mkdir(d)

def concat_prodigal(genomes):
	msg("Concatting protein predictions")
	emptyfile(config['prodigal']['concat'])
	for genome in genomes:
		os.system('''cat {0}/{1}.proteins.filtered.faa >> {2}'''.format(config['prodigal']['folder'],genome,config['prodigal']['concat']))
	msg("Ready with concatenating prodigal results")



def mmseq_cmd(proteomefile,outprefix):
	cmd = "mmseqs easy-cluster {} {} {}/tmp -c {} --cov-mode {} --min-seq-id {} --threads {}".format(proteomefile,outprefix,config['resdir'],config['mmseqs']['mincov'],config['mmseqs']['covmode'],config['mmseqs']['minseqid'],threads)
	msg(cmd)
	os.system(cmd)




def emptyfile(out,header=None):
	if header:
		with open(out,"w") as g:
			g.write("\t".join(map(str,header))+"\n")
	else:
		if os.path.isfile(out):
			os.remove(out)
		os.system("touch " + out)

def alignment(cluster,path=""):
	if path == "":
		path = config['clusters']['path']
	msg("MAFFT started for {}".format(cluster))
	sysexec('mafft --quiet --thread 1 --threadtb 1 --threadit 1 --retree 1 --maxiterate 0 {1}/{0}.fna > {1}/{0}.aln'.format(cluster,path))

def concat_snpsites(clusters,genomes):
	msg("Concatenating clusters to {} and to {}".format(config['snpsites'],config['dropped_recombinant']))
	concatdic = defaultdict(str)
	concatdic_keep = defaultdict(str)
	for cl in clusters:
		cl_keep = not bool(int(read_cmd("wc -l {}/{}.recombinant_regions.strict.bed".format(config['clusters']['path'],cl)).split()[0]))
		##if cl.endswith("000073"):
		##	print(read_cmd("wc -l {}/{}.recombinant_regions.strict.bed".format(config['clusters']['path'],cl)).split())
		##	print(cl_keep)
		parsimony_aln = "{}/{}.masked.trimmed.parsimony_sites.aln".format(config['clusters']['path'],cl)
		logging.info("Concatenation: processing {}".format(cl))
		clusterdic = {}
		if not os.path.isfile(parsimony_aln):
			continue
		for rcd in SeqIO.parse(parsimony_aln,"fasta"):
			try:
				clusterdic[rcd.id.split("|")[0]] = str(rcd.seq).upper()
			except:
				print(rcd)
		alnhossz = len(rcd.seq)
		for key in genomes:
			if key in clusterdic:
				seq = clusterdic[key]
			else:
				seq = "-"*alnhossz
			if cl_keep:
				concatdic_keep[key] += seq
			concatdic[key] += seq
	hosszak = list({len(concatdic[key]) for key in genomes})
	hosszak2 = list({len(concatdic_keep[key]) for key in genomes})
	msg("Length of alignment: {}".format(" ".join(map(str,hosszak))))
	msg("Length of alignment (dropped recombinants): {}".format(" ".join(map(str,hosszak2))))
	with open(config['snpsites'],"w") as g:
		for key in genomes:
			g.write(">{}\n{}\n".format(key,concatdic[key]))
	with open(config['dropped_recombinant'],"w") as g:
		for key in genomes:
			g.write(">{}\n{}\n".format(key,concatdic_keep[key]))

def add_to_version(subcmd,stderr=False,newline=True):
	if stderr:
		os.system(subcmd + " 2>> " + config['versionlog'])
	else:
		os.system(subcmd + " >> " + config['versionlog'])
	if newline:
		os.system("echo >> " + config['versionlog'])

def get_versions():
	emptyfile(config['versionlog'])
	add_to_version("prodigal -v",True)
	add_to_version("echo MMseqs2",newline=False)
	add_to_version("mmseqs version")
	add_to_version("echo MAFFT",newline=False)
	add_to_version("mafft --version",True)
	add_to_version("trimal --version")
	add_to_version("run_treeshrink.py -v")

def replace_wildcard_in_dic(dic,old,new):
	outdict = {}
	for k,v in dic.items():
		if type(v) == dict:
			outdict[k] = replace_wildcard_in_dic(v,old,new) #Itt a rekurzio
		elif type(v) == str and old in v:
			outdict[k] = v.replace(old,new)
		else:
			outdict[k] = v
	return outdict




def parse_config():
	print("Parsing " + pargs.configfile)
	dic = yaml.safe_load(open(pargs.configfile))
	dic = replace_wildcard_in_dic(dic,"{projectname}",dic["projectname"])
	dic = replace_wildcard_in_dic(dic,"{resdir}",dic["resdir"])
	return dic

def prodigal(genome):
	if os.path.isfile("{}/{}.proteins.filtered.faa".format(config['prodigal']['folder'],genome)):
		return
	print("Prodigal for "+ genome)
	os.system("prodigal -i {1}/{0}.fna -a {2}/{0}.proteins.nonpartial.faa -c -m -d {2}/{0}.proteins.nonpartial.fna -f gff -o {2}/{0}.proteins.nonpartial.gff -q".format(genome,config['input']['genomedir'],config['prodigal']['folder']))
	logging.info("Ready with prodigal (nonpartial) for " + genome)
	msg("Copying proteins of " + genome)
	sysexec("seqtk comp {0}/{1}.proteins.nonpartial.faa | awk '{{if ($2>={2}) print $1}}' > {0}/{1}.orf.list".format(config['prodigal']['folder'],genome,config['prodigal']['minorf']))
	rcd_id_dict = {}
	with StringIO(read_cmd("seqtk subseq {0}/{1}.proteins.nonpartial.faa {0}/{1}.orf.list".format(config['prodigal']['folder'],genome))) as fasta_io, open("{}/{}.proteins.filtered.faa".format(config['prodigal']['folder'],genome),"w") as g:
		ctg_names = defaultdict(lambda:0)
		for rcd in SeqIO.parse(fasta_io, "fasta"):
			ctg = rcd.id.rsplit("_",1)[0]
			ctg_names[ctg] += 1
			rcd_id_dict[rcd.id] = "{}|ORF_{:06}".format(ctg,ctg_names[ctg])
			g.write(">{}\n{}\n".format(rcd_id_dict[rcd.id],rcd.seq))
	with StringIO(read_cmd("seqtk subseq {0}/{1}.proteins.nonpartial.fna {0}/{1}.orf.list".format(config['prodigal']['folder'],genome))) as fasta_io, open("{}/{}.proteins.filtered.fna".format(config['prodigal']['folder'],genome),"w") as g:
		g.writelines([">{}\n{}\n".format(rcd_id_dict[rcd.id],rcd.seq) for rcd in SeqIO.parse(fasta_io, "fasta")])

	


def get_my_genomes():
	msg("Reading {}".format(config['sample']['file']))
	return read_cmd("cut -f 1 {} | tail -n +2".format(config['sample']['file'])).split()

def collect_mmseqs_clusters(genomes): #Ez az eljaras nagyon memoriapazarlo, ha defaultdict-tel fut - azt el kell kerulni, mert megeszi az egesz memoriat!!
	print("Creating cluster prevalence tables")
	clusterfile = config['mmseqs']['prefix'] + "_cluster.tsv"
	#!cluster_aa_seqs = config['mmseqs']['prefix'] + "_all_seqs.fasta"
	protdir = config['prodigal']['folder']
	cluster_dic = {}
	prevalence_dic = {}
	msg("Reading " + clusterfile)
	with open(clusterfile,"r") as f:
		rdr = csv.reader(f,delimiter="\t")
		for row in rdr:
			if len(row)<2:
				continue
			if row[0] not in cluster_dic:
				cluster_dic[row[0]] = set()
			if row[0] not in prevalence_dic:
				prevalence_dic[row[0]] = {}
			genome = row[1].split("|")[0]
			if genome not in prevalence_dic[row[0]]:
				prevalence_dic[row[0]][genome] = 0
			cluster_dic[row[0]].add(row[1])
			prevalence_dic[row[0]][genome] += 1
	clusterreps = list(cluster_dic.keys())
	clusterreps.sort(key = lambda x: -len(cluster_dic[x]))
	msg("Writing cluster prevalence tables")
	with open(config['clusters']['clustermembers'],'w') as g, open(config['clusters']['prevalence'],'w') as h:
		g.write("Cluster ID\tNumber of genomes\tPercentage of genomes\tNumber of ORFs\tRepresentant\tCluster members\n")
		h.write("Cluster ID\tNumber of genomes\tPercentage of genomes\tNumber of ORFs\t" + "\t".join(genomes) + "\n")
		ix = 0
		for rep in clusterreps:
			ix += 1
			if ix%10000 == 0:
				msg("Now processing cluster nr. {}".format(ix))
			clname = "Cluster_{:07}".format(ix)
			prevrow = []
			genomenum = 0
			for genome in genomes:
				if genome not in prevalence_dic[rep]:
					prevrow.append(0)
				else:
					prevrow.append(prevalence_dic[rep][genome])
					genomenum += 1
			h.write("\t".join(map(str,[clname,genomenum,genomenum/len(genomes)*100,sum(prevrow)] + prevrow)) + "\n")
			g.write("\t".join(map(str,[clname,genomenum,genomenum/len(genomes)*100,len(cluster_dic[rep]),rep,";".join(sorted(cluster_dic[rep]))])) + "\n")

def core_clusters():
	sysexec(r"awk -F '\t' '{if ($3>=99) print $0}' " + config['clusters']['clustermembers'] + " > " + config["clusters"]["core_details"])
	sysexec(r"awk -F '\t' '{if ($3>=99) print $0}' " + config['clusters']['prevalence'] + " > " + config["clusters"]["core_prevalence"])
	sysexec("cut -f 1 {} | tail -n +2 > {}".format(config['clusters']['core_details'],config['clusters']['core_ids']))

def read_core_clusters():
	msg("Reading core clusters")
	cluster_orf_pairs = []
	with open(config['clusters']['core_details']) as f:
		f.readline()
		for line in f.readlines():
			row = line.strip().split()
			cluster_orf_pairs.append((row[0],row[5],))
	return cluster_orf_pairs

def collect_cluster_sequences(cluster,orfs):
	msg("Collecting sequences for " + cluster)
	clfile = "{}/{}.fna".format(config['clusters']['path'],cluster)
	emptyfile(clfile)
	for orf in orfs.split(";"):
		sysexec("grep -A 1 '{}' {}/{}.proteins.filtered.fna >> {}".format(orf,config['prodigal']['folder'],orf.split("|")[0],clfile))
	msg("Ready with " + cluster)



def build_tree(aln,tree,outprefix,outaln,threads="AUTO",mode="iqtree"):
	msg("Building tree for {}".format(aln))
	if mode == "iqtree":
		cmd = "iqtree2 -m GTR+F+R3 -T {0} -s {1} > {1}.iqtree.log".format(threads,aln)
		sysexec(cmd)
	elif mode == "rapidnj":
		sysexec("rapidnj -c {} {} -i fa > {}".format(threads,aln,tree))
	sysexec(r"sed -i s/\'//g {}".format(tree))
	sysexec("run_treeshrink.py -t {} -a {} -o {} -O {}".format(tree,aln,config['resdir'],outprefix))
	msg("Remove dropped sequences from the alignment")
	with open("{}/{}.txt".format(config['resdir'],outprefix)) as f:
		excluded_genomes = f.read().strip().split()
	with open(outaln,"w") as g:
		for rcd in SeqIO.parse(aln,"fasta"):
			if rcd.id not in excluded_genomes:
				g.write(">{}\n{}\n".format(rcd.id,rcd.seq))
	msg("Ready with {}/{}.nwk and {}".format(config['resdir'],outprefix,outaln))


def rename_genome(row,jobname_index,path_index):
	if os.path.isfile("{}/{}.fna".format(config['input']['genomedir'],row[0])):
		return
	msg("Copying " + row[0])
	with open("{}/{}.fna".format(config['input']['genomedir'],row[0]),'w') as g:
		i = 0
		for rcd in SeqIO.parse("{}/{}/{}".format(config['input']['source'],row[jobname_index],row[path_index]),'fasta'):
			i += 1
			g.write(">{}|C{:04}\n{}\n".format(row[0],i,rcd.seq))


def sample_from_xdr_genomes(table,sample_size,outf):
	msg("Selecting {} XDR genomes from {}".format(sample_size,table))
	xdr_genomes = []
	with open(table) as f, open(outf,"w") as g:
		rdr = csvrdr(f)
		wtr = csv.writer(g,delimiter="\t")
		header = next(rdr)
		wtr.writerow(header)
		xdr_index = header.index("xdr")
		jobname_index = header.index("jobname")
		path_index = header.index("relative_path")
		for row in rdr:
			if row[xdr_index] == "TRUE":
				xdr_genomes.append(row)
		random.seed("fenyotoboz")
		rnd = random.sample(xdr_genomes,sample_size)
		wtr.writerows(rnd)
	msg("{} is ready".format(outf))
	mkdir_force(config['input']['genomedir'])
	pool = mp.Pool(threads)
	pool.starmap(rename_genome,[(row,jobname_index,path_index,) for row in rnd])

def read_clusters(fname):
	msg("Reading clusters from {}".format(fname))
	with open(fname) as f:
		return sorted({line.strip() for line in f if len(line.strip())})


def is_ready_fastgear(cl):
	return os.path.isfile("{1}/{0}_fastgear/output/lineage_information.txt".format(cl,config['fastgear']['folder']))

def fastgear(cl):
	if is_ready_fastgear(cl):
		msg("FastGEAR is already finished for {}".format(cl))
	else:
		start = time.time()
		mkdir_force("{}/{}_fastgear".format(config['fastgear']['folder'],cl))
		####sysexec("singularity exec -B /node10_R10:/node10_R10 -B /node8_R10:/node8_R10 -B /node8_data:/node8_data /home/vasarhelyib/docker/hgttree_latest.sif cpulimit - -m -l 200 -f -q run_fastGEAR.sh ./results/clusters/{0}.aln ./FastGEAR/{0}_fastgear/{0}.mat fG_input_specs.txt > FastGEAR/{0}.log 2> FastGEAR/{0}.err".format(cl))
		sysexec("run_fastGEAR.sh ./{2}/{0}.aln ./{1}/{0}_fastgear/{0}.mat {3} > {1}/{0}.log 2> {1}/{0}.err".format(cl,config['fastgear']['folder'],config['clusters']['path'],config['fastgear']['settings']))
		if is_ready_fastgear(cl):
			msg("Ready with fastGEAR for {} in {} seconds".format(cl,time.time()-start))
		else:
			msg("ERROR with FastGEAR for {} after {} seconds".format(cl,time.time()-start))


def collect_recombination_data(cluster):
	path = "{}/{}_fastgear/output".format(config['fastgear']['folder'],cluster)
	cluster_path = config['clusters']['path']
	rrname = "{}/recombinations_recent.txt".format(path)
	recombinant_bed = "{}/{}.recombinant_regions.bed".format(cluster_path,cluster)
	recombinant_bed_strict = "{}/{}.recombinant_regions.strict.bed".format(cluster_path,cluster)
	aln = "{}/{}.aln".format(cluster_path,cluster)
	masked_aln = "{}/{}.masked.aln".format(cluster_path,cluster)
	no_recom_aln = "{}/{}.dropped_recombinant.aln".format(cluster_path,cluster)
	trimmed_aln = "{}/{}.masked.trimmed.aln".format(cluster_path,cluster)
	parsimony_aln = "{}/{}.masked.trimmed.parsimony_sites.aln".format(cluster_path,cluster)
	msg('Processing recombination data for {}'.format(cluster))
	if os.path.isfile(rrname):
		outrows = []
		with open(rrname) as f:
			next(f)
			next(f)
			for line in f:
				block = line.strip().split()
				outrows.append([block[5], int(block[0])-1, int(block[1]), "recent"])
		ancestral = []
		with open("{}/recombinations_ancestral.txt".format(path)) as f:
			next(f)
			next(f)
			for line in f:
				block = line.strip().split()
				ancestral.append([block[0], int(block[1])-1, int(block[3])])
		if ancestral:
			with open("{}/lineage_information.txt".format(path)) as f:
				next(f)
				strains_in_lineages = defaultdict(list)
				for line in f:
					block = line.strip().split()
					strains_in_lineages[block[1]].append(block[3])
			for item in ancestral:
				for strain in strains_in_lineages[item[-1]]:
					outrows.append([strain] + item[:2] + ["ancestral"])
		with open(recombinant_bed,"w") as g:
			wtr = csv.writer(g,delimiter="\t")
			wtr.writerows(outrows)
		with open(recombinant_bed_strict, "w") as g:
			wtr = csv.writer(g,delimiter="\t")
			for row in outrows:
				if int(row[2]) - int(row[1]) >= config['fastgear']['min_recomb']:
					wtr.writerow(row)
	else:
		emptyfile(recombinant_bed)
	if not os.path.isfile(recombinant_bed_strict):
		emptyfile(recombinant_bed_strict)
		msg("{} is empty".format(recombinant_bed_strict))
	msg('Ready with {} and {}'.format(recombinant_bed, recombinant_bed_strict))
	'''Mask fasta file with a given bed file using bedtools'''
	sysexec("bedtools maskfasta -fi {} -bed {} -fo {}".format(aln,recombinant_bed_strict,masked_aln))
	sysexec("trimal -in {} -out {} -gappyout -keepheader".format(masked_aln,trimmed_aln,path))
	msg("Ready with aligning {}".format(cluster))
	msg("Filter SNP sites for " + cluster)
	aln = SeqIO.parse(trimmed_aln,"fasta")
	rcd = next(aln)
	totalnum = len(rcd)
	good_coords = []
	unknown_coords = list(range(len(rcd)))
	coord_nuc = defaultdict(lambda:defaultdict(lambda:0))
	rcd.seq = rcd.seq.upper()
	for coord in unknown_coords:
		if rcd.seq[coord] not in ["-","N"]:
			coord_nuc[coord][rcd.seq[coord]] += 1
	ids = [rcd.id]
	sequences = [rcd.seq]
	for rcd in aln:
		rcd.seq = rcd.seq.upper()
		ids.append(rcd.id)
		sequences.append(rcd.seq)
		new_unknown_coords = []
		for coord in unknown_coords:
			if rcd.seq[coord] not in ["-","N"]:
				coord_nuc[coord][rcd.seq[coord]] += 1
			if len(coord_nuc[coord].keys())>=3 or (len(coord_nuc[coord].keys())>=2 and min(coord_nuc[coord].values())>=2):
				good_coords.append(coord)
			else:
				new_unknown_coords.append(coord)
		if len(new_unknown_coords) == 0:
			break
		else:
			unknown_coords = new_unknown_coords
	for rcd in aln:
		ids.append(rcd.id)
		sequences.append(rcd.seq)
	good_coords = sorted(good_coords)
	msg("{}: {} of {} coordinates were kept".format(cluster,len(good_coords),totalnum))
	with open(parsimony_aln, "w") as g:
		for rcdid,seq in zip(ids,sequences):
			g.write(">{}\n{}\n".format(rcdid, "".join([seq[coord] for coord in good_coords])))
	







if __name__ == "__main__":
	logging.basicConfig(format='%(asctime)s. %(levelname)s:%(message)s',filename="{}.{}.log".format(__file__.replace(".py3",""),time.strftime("%Y%m%d%H%M%S")),level=logging.INFO)
	logging.Formatter(fmt='%(asctime)s',datefmt="%Y")
	logging.info("Run started")
	logging.info("Reading " + pargs.configfile)
	config = parse_config()
	mkdir_force(config['resdir'])
	if pargs.n:
		threads = pargs.n
	elif 'threads' in config:
		threads = config['threads']
	else:
		threads = 1
	pool = mp.Pool(threads)
	msg("Using {} threads".format(threads))
	get_versions()
	mkdir_force(config['resdir'])
	'''
	1.
	Collect genome IDs and rename contigs
	'''
	mkdir_force(config['input']['genomedir'])
	sample_from_xdr_genomes(config['input']['aci_tbl'],config['sample']['size'],config['sample']['file'])
	genomes = get_my_genomes()
	'''
	2. Prodigal
	'''
	mkdir_force(config['prodigal']['folder'])
	prodigal_ready = set(map(lambda x: x.split("/")[-1].replace(".proteins.filtered.faa",""), glob.glob("{}/*.proteins.filtered.faa".format(config['prodigal']['folder']))))
	genomes_prodigal_to_do = set(genomes) - prodigal_ready
	msg("I'll do prodigal for {} genomes". format(len(genomes_prodigal_to_do)))
	pool.map(prodigal,genomes_prodigal_to_do)
	concat_prodigal(genomes)
	'''
	3. MMseqs2 (Many against Many sequence searching) is an open-source software suite for very fast, 
	parallelized protein sequence searches and clustering of huge protein sequence data sets.
	M. Steinegger and J. Soding. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi:10.1038/nbt.3988 (2017).
	MMseqs2 Version: 0c4b749805c25f652de6f137ea204155bf78b1d2'''
	mmseq_cmd(config['prodigal']['concat'],config['mmseqs']['prefix'])
	'''4.
	We collect the core sequence identifiers and the sequences and clusters in a separate file'''
	'''Algorithm:
	i. Read cluster details from mmseqs
	ii. Get core clusters (those which prevail in at least 99% of the genomes)
	iii. Drop those core clusters which match resistance or virulence genes
	iv. Check prevalence table: drop those genomes which have paralogs in at least ... clusters and those clusters which have paralogs from more than 1% of the associated genomes - first find distribution
	v. Collect the sequences of the remaining clusters
	vi. Align, trim and find parsimonious sites.
	'''	
	print("Collect cluster members")
	collect_mmseqs_clusters(genomes)
	core_clusters()
	msg("Collecting core cluster sequences")
	mkdir_force(config['clusters']['path'])
	cluster_orflist = read_core_clusters()
	pool.starmap(collect_cluster_sequences,[(cl,orf,) for cl,orf in cluster_orflist]) #Collect core cluster sequences
	'''
	5. 
	Alignment with mafft for core clusters
	
	MAFFT v7.475 (2020/Nov/23)
	MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
	https://mafft.cbrc.jp/alignment/software/
	
	Trimming the alignments:
	
	trimAl v1.4.rev22 build[2015-05-21]. 2009-2015. Salvador Capella-Gutierrez and Toni Gabaldon.
	trimAl webpage: http://trimal.cgenomics.org	
	trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses.
	Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon.
	Bioinformatics 2009, 25:1972-1973.
	'''
	clusters = open(config['clusters']['core_ids']).read().strip().split()
	pool.map(alignment,clusters)
	'''RUN fastgear. '''
	pool2 = mp.Pool(config['fastgear']['threads'])
	mkdir_force(config['fastgear']['folder']
	pool2.map(fastgear,clusters)
	pool2.map(fastgear,clusters) ###ha nem futott le a fastgear, akkor még egyszer próbálkozunk
	'''
	MASK fastgear recombinant sites, TRIMAL and filtering PARSIMONIOUS SITES for all clusters
	'''
	pool.map(collect_recombination_data,clusters)
	'''6. Build the tree with rapidnj and run treeshrink
	https://pure.au.dk/ws/files/19821675/rapidNJ.pdf	
	'''
	concat_snpsites(clusters,genomes)
	build_tree(config['snpsites'],config['tree']['masked_rapidnj'],config['treeshrink']['prefix_masked_rapidnj'],config['treeshrink']['aln_masked_rapidnj'],threads,"iqtree")
