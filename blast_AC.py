#!/usr/bin/env python
#learning how to use blast and biopython stuff


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
N_worst = 5 
#gene_names =[] 

#want to create a list of dictionaries, 1 dictionary for each gene used in MLST
#then when I blast the gene, the keys will be the strains of interest, and the values will be the sequences I get from blasting
# with open("mlst-genes.txt","rU") as infile2_handle:
# 	dictGenes = [dict() for x in SeqIO.parse(infile2_handle,"fasta")]
# 	print dictGenes

counter=0
#with open("adk.txt","rU") as infile_handle: 
with open("mlst-genes.txt","rU") as infile_handle:
	#dictGenes = [dict() for x in SeqIO.parse(infile_handle,"fasta")]
	#can only do SeqIO.parse() once
#The rU means universal readline mode so do not need to worry about differences in mac/unix/pc newline characters
	for search_record in SeqIO.parse(infile_handle, "fasta"):
		dictGenes={}
		#gene_names.append(search_record.name)
		#dictGenes.append(dict())
		print "************************************************************"
		print search_record.name,"\n"
		f_query = search_record.seq
#		fout = open("blast_hit_list", "w")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Note - this is where filtering based on identity cutoff or E_value is set
		#E_VALUE_THRESH = 1.0e-1
		identity_cutoff = 0.9
		print "identity_cutoff=",identity_cutoff,"\n"

#3 required arguments for qblast: blast program, database to search, string with query sequence
		result_handle = NCBIWWW.qblast("blastn", "Complete_Genomes", f_query, hitlist_size = 100, entrez_query=\
'txid511145[orgn] OR txid331112[orgn] OR txid316435[orgn] OR txid1314835[orgn] OR txid574521[orgn] OR txid216592[orgn] OR txid199310[orgn] OR txid754331[orgn] OR txid656449[orgn] OR txid316401[orgn] OR txid155864[orgn] OR txid444450[orgn] OR txid386585[orgn] OR txid588858[orgn] OR txid83334[orgn]')

#use NCBIXML.read() for 1 object and NCBIXML.parse() for multiple objects. NCBIXML.parse() returns an iterator
		hit_strains = " "
		blast_records = NCBIXML.parse(result_handle)

		#save_file = open("my_blast5.xml","w")
		#save_file.write(result_handle.read())
		#save_file.close()
		
		for blast_record in blast_records:
			alignment_list = [] #will hold list of alignments in order of lowest HSP to highest HSP
			for alignment in blast_record.alignments:
				hsp = alignment.hsps[0] #alignment.hsps[0] gives the HSP score
				if len(alignment_list) > 1:
					for i in range(len(alignment_list)):
					#will compare identities (number of nucleotides that match) to that of each alignment in alignment_list
						if hsp.identities < alignment_list[i].hsps[0].identities:
							alignment_list.insert(i,alignment) #adds alignment with lower identity before the next alignment and then break out of loop
							break
					else:
						alignment_list.append(alignment) #if alignment doesn't have lower HSP than any in the list, it gets added to the end of the list
				elif len(alignment_list) == 1:
					if hsp.identities < alignment_list[0].hsps[0].identities:
						alignment_list.insert(0,alignment) # insert the new alignment at position 0 if alignment has lower identity than first item in list
					else:
						alignment_list.append(alignment)
				else:
				#len(alignment_list) has 0 items. adds first alignment to list. 
					alignment_list.append(alignment)
					
#adds alignment to hit_strains if the fraction of identities meets the identity_cutoff
				if float(hsp.identities)/float(len(f_query)) > identity_cutoff:
					hit_strains += alignment.title
					dictGenes[alignment.title]=alignment.hsps[0].sbjct
				else:
					print "This entry has no match with percent identity > ",identity_cutoff,":\n",alignment.title
					print "query length:",len(f_query)," identities:",hsp.identities," percent:",float(hsp.identities)/float(len(f_query))
					print "E value:",hsp.expect,"\n"
		print "------------------------------------------------------------\n"
		
		for i in range(N_worst):
			print alignment_list[i].title
			print alignment_list[i].hsps[0]
			#print alignment_list[i].hsps[0].match
			print alignment_list[i].hsps[0].sbjct
			print alignment_list[i].hsps[0].identities,"identities out of ",alignment_list[i].hsps[0].align_length,"\n"
		counter = counter+1

		gene_name = search_record.name+".txt"
		outFile = open(gene_name,"w")
		for k in dictGenes:
			outFile.write(k)
			outFile.write("\t")
			outFile.write(dictGenes[k])
			outFile.write("\n")
		outFile.close()			