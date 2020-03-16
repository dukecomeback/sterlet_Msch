#!env python

#
# parse_cd-hit_count-occurrences.py
# Copyright (C) 2018 INRA
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Plateforme bioinformatique Midi Pyrenees'
__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'beta'

import configparser
from optparse import *
import datetime, os, sys, re, string
import pysam
from Bio import SeqIO

def version_string ():
	"""
	Return extract_indel_correted_reads_from_bam.py version
	"""
	return "extract_indel_correted_reads_from_bam.py " + __version__

def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def splitHaplotypeFiles(locations, h1File, h2File) :
	"""
	Returns a dictionnary with the reference sequences as a string
	"""
	length = 80

	# load fasta files 
	h1 = {}
	h2 = {}
	len_h1 = 0
	len_h2 = 0
	for seq_record in SeqIO.parse(h1File, "fasta") :
		h1[seq_record.id] = str(seq_record.seq)
		len_h1 = len(str(seq_record.seq))
	for seq_record in SeqIO.parse(h2File, "fasta") :
		h2[seq_record.id] = str(seq_record.seq)
		len_h2 = len(str(seq_record.seq))

	# check if the haplotype files are single sequences 
	if len(h1) > 1 :
		sys.stderr.write('ERROR : the h1 haplotype fasta file contains more than one sequence ' + str(nuc) + "\n")
		exit()
	if len(h2) > 1 :
		sys.stderr.write('ERROR : the h2 haplotype fasta file contains more than one sequence ' + str(nuc) + "\n")
		exit()	

	# create the location files 
	h1pos = []
	h2pos = []
	j = 0
	for i, l in enumerate(locations) :
		if i == 0 :
			h1pos.append([0,int(l)])
		else :
			h1pos.append([j,int(l)])
		j = int(l) + 1
	# write last h1 block positions 
	h1pos.append([j,len_h1])
	for i, l in enumerate(locations) :
		if i == 0 :
			h2pos.append([0,int(l)])
		else :
			h2pos.append([j,int(l)])
		j = int(l) + 1
	# write last h1 block positions 
	h2pos.append([j,len_h2])		
	# print(h1pos)
	# print(h2pos)

	# extract fasta files 
	# rename h1 and h2 files 
	h1_name = os.path.splitext(os.path.basename(h1File))[0]+".split.fasta"
	h1_handle = open(h1_name, "w")
	h2_name = os.path.splitext(os.path.basename(h2File))[0]+".split.fasta"
	h2_handle = open(h2_name, "w")
	# write H1
	for ind, i in enumerate(h1pos) :
		# print(i)
		#print(h1[seq_record.id])
		h1_handle.write(">"+seq_record.id+"_H1_"+str(ind)+"\n")
		for s in split_len(h1[seq_record.id][i[0]:i[1]], length) :
			h1_handle.write(s+"\n")
	h1_handle.close()
	# write H2
	for ind, i in enumerate(h2pos) :
		h2_handle.write(">"+seq_record.id+"_H2_"+str(ind)+"\n")
		for s in split_len(h2[seq_record.id][i[0]:i[1]], length) :
			h2_handle.write(s+"\n")
	h2_handle.close()

	return True

def getVariationLocations(variations, vcfFile) :
	"""
	retrieves the variation location for each of the variation index  
	"""
	locations = []
	vcf_handle = open(vcfFile, "r")
	i = 0
	for line in vcf_handle:
		if line[0:1] != "#" :
			# print(str(i)+"\t"+line)
			i += 1
			for j in variations :
				if int(j) == int(i) :
					# print("HERE")
					s = line.split('\t')
					locations.append(s[1])
	vcf_handle.close()
	#print(locations)
	return locations

def getSplitVariations(hapcut2File) :
	"""
	returns all the split variations of the hapcut2 output file 
	Lines starting with 
	"""
	variations = []
	hapcut2_handle = open(hapcut2File, "r")
	i = 0
	for line in hapcut2_handle:
		i += 1
		#print(line)
		# print("here:\t"+line[0:6]+"\t*"+line[7:19]+"*")
		if line[0:6] == "BLOCK:" and line[7:19] != "(from split)" :
			s = line.split(' ')
			if i != 1 :
				variations.append(s[2])
	hapcut2_handle.close()
	#print("variation "+variations)
	return variations

if __name__ == "__main__":

	parser = OptionParser(usage="Usage: split_haplotype_sequences_using_VCF_and_haplotype_output.py -a HAPCUT2-OUTPUT -v VCF-FILE -1 HAPLOTYPE1 -2 HAPLOTYPE2")

	usage = "usage: %prog -a hapcut2_output_file -v vcf_file -1 haplotype1 -2 haplotype2"
	desc = " produces multi fasta file for haplotype1 and haplotype2 fasta files splitting them using the hapcut2 output file \n"\
		   " ex : esplit_haplotype_sequences_using_VCF_and_haplotype_output.py -a hapcut2_output_file -v vcf_file -1 haplotype1 -2 haplotype2"
	parser = OptionParser(usage = usage, version = version_string(), description = desc)
	
	ogroup = OptionGroup(parser, "Options","")
	ogroup.add_option("-a", "--hapcut2-output-file", dest="hapcut2File",
					  help="Indicate input hapcut2 output file",  type="string")
	ogroup.add_option("-1", "--haplotype1-fasta-file", dest="h1FastaFile",
					  help="Indicate input haplotype 1 fasta file",  type="string")
	ogroup.add_option("-2", "--haplotype2-fasta-file", dest="h2FastaFile",
					  help="Indicate input haplotype 2 fasta file",  type="string")
	ogroup.add_option("-v", "--vcf-file", dest="vcfFile",
					  help="Indicate the input VCF file",  type="string")
	parser.add_option_group(ogroup)
	
	(options, args) = parser.parse_args()
	if options.hapcut2File == None and options.vcfFile == None and  options.h1FastaFile == None and options.h2FastaFile == None:
		parser.print_help()
		sys.exit(1)
	else:
	# process files 
		variations = getSplitVariations(options.hapcut2File)
		variationLocations = getVariationLocations(variations, options.vcfFile)
		splitHaplotypeFiles(variationLocations, options.h1FastaFile, options.h2FastaFile)
