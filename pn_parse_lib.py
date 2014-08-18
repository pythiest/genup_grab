#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pn_parse_lib.py
#  
#  Copyright 2014 <k7n7vi@ivan-UbuVBox>
#  
#  This is a library to basically get a list of NCBI protein GIs and
#  use eUTILS and biopython to automatically retrieve a desired 
#  (upstream) nucleotide sequence from the corresponding CDS in GenBank
#  For now designed to work with "old" RefSeq and GenBank IDs (e.g. YP_)
#  not new WP_ IDs
#  
#  The library contains three main functions:
#  - readGIs reads a list of GIs from file
#  - get_nuc_acc receives a protein GI and retrieves the CDS field
#    position feature (nucleotide seq, plus positions therein)
#    this function is duplicate as SeqIO_get_nuc_acc using the SeqIO
#    functionality to get the same info
#  - parse_nuc_acc gets the CDS field feature and parses it


#imports the Entrez module from the Bio (biopython) library
#Bio is the official name for biopython
#we also set the basic entrez parameters (e.g. email address)
from Bio import Entrez, SeqIO
Entrez.email ="ivan.erill@gmail.com"


#system library for closing the process if on error
import sys

#time library for "sleeping" between NCBI queries
import time


#---------------------------------------------------------------
def read_GIs(file_name = "GI_list.txt"):
    """opens file for reading and gets the GIs into a list, assuming 
    that there is one GI per line.    
    """
    #create file handler for reading file
    try:
        GI_file = open(file_name,"r")
    except (IOError, OSError) as file_open_exception:
        print "The file name provided:", file_name, " does not exist"
        sys.exit()
        
      
    #go through the file, reading GIs into a list
    #we first create an empty GI list
    GIs = []
    
    #we iterate the file as a list of lines
    for line in GI_file:
        #remove \n and empty characters from read line
        #(by default, the iteration reads \n as independent lines)
        #rstrip on a "\n" line will essentially make it an empty string
        line=line.rstrip()
        
        #if line is not empty, then
        if line:
            #if line is numeric
            if (line.isdigit()):
                #append read line to GI list
                GIs.append(line)
    
    GI_file.close()
    #return GIs as list
    return GIs

#---------------------------------------------------------------
def parse_nuc_accession(nuc_acc):
    """parses a nucleotide accession with positions
    basically, two possibilities are contemplated:
    1) forward: ID:pos..pos
    2) reverse: complement(ID:pos..pos) 
    
    returns (using efetch standards):
    - genome accession
    - the strand (1 for forward, 2 for reverse)
    - the sequence start
    - the sequence end
    """
    
    #check for "complement" and assign strand
    #also remove complement from string, and parentheses
    if nuc_acc[0:10]=="complement":
        seq_strand="2"
        nuc_acc=nuc_acc.lstrip("complement(")
        nuc_acc=nuc_acc.rstrip(")")
    else:
        seq_strand="1"
    
    if nuc_acc[0:4]=="join":
        return -1,-1,-1,-1    
        
    #get genome accession and positions
    gen_acc, gen_pos=nuc_acc.split(":")
    
    #get start and stop positions
    seq_start, seq_stop = gen_pos.split("..")
    
    return gen_acc, seq_strand, seq_start, seq_stop   


#---------------------------------------------------------------
def get_nuc_acc(prot_GI):
    """gets a protein GI and uses eUTILS to navigate NCBI and fetch the
    associated CDS nucleotide accession GI with coordinates
    """
    #get handle to DB object (protein DB)
    prot_handle = Entrez.efetch(db="protein", id=prot_GI, retmode="xml")
    #donwload the handler associated data
    prot_records = Entrez.read(prot_handle)
    prot_handle.close()
    time.sleep(4)  #sleep for 5 seconds
    
    #get the feature list from first element in the records list
    features=prot_records[0]["GBSeq_feature-table"]
    
    #for the features, scan until detecting a GBFeature_key equal to 'CDS'
    #we assume that the CDS feature is unique for a protein record
    for a_feature in features:
        if a_feature["GBFeature_key"]=="CDS":
            feature_qualifiers=a_feature["GBFeature_quals"]
            for f_qual in feature_qualifiers:
                if f_qual["GBQualifier_name"]=="coded_by":
                    genome_loc=f_qual["GBQualifier_value"]
                    return genome_loc
                    
    #return "null" if feature not found
    return None
               

#---------------------------------------------------------------
def SeqIO_get_nuc_acc(prot_GI):
    """gets a protein GI and uses eUTILS to navigate NCBI and fetch the
    associated CDS nucleotide accession GI with coordinates
    Uses SeqIO to parse the protein record and extract the CDS + position
    """
    #get handle to DB object (protein DB) in GB/text format
    prot_handle = Entrez.efetch(db="protein", id=prot_GI, rettype="gb", 
                                retmode="text")
    #donwload the handler associated data with GB parsers
    prot_record = SeqIO.read(prot_handle,"gb")
    prot_handle.close()
    time.sleep(5)  #sleep for 5 seconds

    #the SeqIO parser provides us a neat protein record split into 
    #annotations, features, format, sequence... we use the features list
    #to get to the CDS
    features=prot_record.features
        
    #for the features, scan until detecting a feature type equal to 'CDS'
    #we assume that the CDS feature is unique for a protein record
    for a_feature in features:
        if a_feature.type=="CDS":
            feature_qualifiers=a_feature.qualifiers
            #return with "get" the "coded_by" qualifier 
            #so that it returns none if the qualifier does not exist
            genome_loc=feature_qualifiers.get("coded_by")
            return genome_loc[0]        

    #return "null" if feature not found
    return None
            
#---------------------------------------------------------------
def SeqIO_get_nuc_rec(nuc_acc):
    """gets a list with gen_acc, seq_strand, seq_start and seq_stop
       performs a query to Entrez to retrieve the nucleotide record
       at those positions in FASTA format
       Uses SeqIO read to parse the object
    """
    nuc_handle = Entrez.efetch(db="nuccore", id=nuc_acc[0], 
                                retmode="text", rettype="gb",
                                strand=nuc_acc[1], seq_start=nuc_acc[2],
                                seq_stop=nuc_acc[3])
    nuc_records = SeqIO.read(nuc_handle,"gb")
    nuc_handle.close()
    time.sleep(5)  #sleep for 5 seconds
    #return record in fasta format
    return nuc_records.format("fasta")
    
#---------------------------------------------------------------
def get_nuc_rec(nuc_acc):
    """gets a list with gen_acc, seq_strand, seq_start and seq_stop
       performs a query to Entrez to retrieve the nucleotide record
       at those positions in FASTA format
    """
    
    #retrieve record with Entrez.read xml
    nuc_handle = Entrez.efetch(db="nuccore", id=nuc_acc[0], 
                                retmode="xml", strand=nuc_acc[1], 
                                seq_start=nuc_acc[2],seq_stop=nuc_acc[3])
    nuc_records = Entrez.read(nuc_handle)
    nuc_handle.close()
    
    time.sleep(5)  #sleep for 5 seconds
    
    #get the GBSeq_definition
    seq_def=nuc_records[0]["GBSeq_definition"]
    #get the sequence
    seq=nuc_records[0]["GBSeq_sequence"]
    #get strand symbol
    strand=""
    if nuc_acc[1]=="2":
        strand="c"
    
    #define fasta header
    fasta_header=">"+nuc_acc[0]+":"+strand+nuc_acc[2]+"-" \
                 +nuc_acc[3]+" "+seq_def
    
    #define fasta field
    fasta_field=fasta_header + "\n" + seq

    #return record in fasta format
    return fasta_field   
    
"""
some explicit code example for efetch/read testing
nuc_handle = Entrez.efetch(db="nuccore", id="NC_010681.1", retmode="text", strand="2", seq_start="4165997",seq_stop="4166596")
nuc_records = Entrez.read(nuc_handle)
    
time.sleep(10)  #sleep for 10 seconds
"""        

#---------------------------------------------------------------
def SeqIO_grab_ext_nuc_rec(nuc_acc, window_size):
    """gets a list with gen_acc, seq_strand, seq_start and seq_stop
    plus a parameter indicating the size of the upstream window to be
    grabbed
    performs a query to Entrez to retrieve the nucleotide record
    at those positions in genbank format
    Uses SeqIO read to parse the object  
    """
    sstart=0
    sstop=0
    sstrand=nuc_acc[1]
    #if forward gene
    if (sstrand=="1"):
        #start of sequence to grab becomes start-window
        sstart=str(int(nuc_acc[2]) - window_size)
        sstop=nuc_acc[3]
    else:
        #end of sequence to grab becomes start+window
        sstop=str(int(nuc_acc[3]) + window_size)
        sstart=nuc_acc[2]

    #the record is retrieved on the "intended" strand, meaning that the
    #gbk record grabbed is effectively on the "forward" strand
    #with the "leading" gene at the end of the record
    #this means that, from now on, both strand cases can be handled
    #equally
    """
       for -o->, we obtain      --> <-- --> -o->
       for <-o-, we have        <-o- <-- --> <--
           but we obtain        --> <-- --> -o->
    """
    nuc_handle = Entrez.efetch(db="nuccore", id=nuc_acc[0], 
                               retmode="text", rettype="gb", strand=sstrand,
                               seq_start=sstart, seq_stop=sstop)
    nuc_records = SeqIO.read(nuc_handle,"gb")
    nuc_handle.close()
    
    time.sleep(4)  #sleep for 5 seconds
    
    #return the extended grabbed record
    return nuc_records    
       
#---------------------------------------------------------------
def SeqIO_extract_operon_up(nuc_rec, up_d, down_d, int_d):
    """gets a nucleotide record, going "back" window size from an
    original CDS derived from a protein GI provided initially
    also gets the up and down parameters relative to the first operon
    gene TLS to extract a window of sequence
    also gets the maximum intergenic distance within an operon
    traces the nucleotide record "back" through CDS records, locating
    the first gene in the operon as the one for which the intergenic
    distance exceeds teh max size, or the last before switching 
    strand.
    it then returns the FASTA sequence corresponding to that region
    adding to the header the set of genes considered operon
    """
    locuses = []

    print "\n" + "------------------------------------------------"
    print "Initiating upstream scan for: " + nuc_rec.id
    last_start=None
    adist=0
    #iterate the features list in reverse
    for cur_feature in reversed(nuc_rec.features):
        #if a CDS is found
        if (cur_feature.type=="CDS"):
            print "CDS found"
            #if we reached the end of the segment, return -1
            #will save the protein ID for which we did not get 
            #an intergenic region (the segment was too small)
            if (str(cur_feature.location).startswith("[<")):
                print "|---->Run out of upstream space."
                if (last_start is None):
                    print "|---->Likely misannotated start."
                    return(-1,-2)
                else:
                    return(-1,-1)
            else:
                #if in proper orientation
                if (cur_feature.strand==1):
                    #if on first CDS
                    if (last_start is None):
                        print "-- First CDS \n--", cur_feature.location, cur_feature.qualifiers["protein_id"][0], cur_feature.qualifiers["product"][0]
                        #add first CDS
                        locuses.append(cur_feature)
                        #set last_start to start of feature
                        last_start=cur_feature.location.start.position
                    else:
                        print "-- Another CDS \n--", cur_feature.location, cur_feature.qualifiers["protein_id"][0], cur_feature.qualifiers["product"][0]
                        #if intergenic region < max allowed (i.e. operon)
                        adist=last_start-cur_feature.location.end.position
                        print "---- Intergenic distance: ", adist
                        if (adist<int_d):
                            print "---- CDS added"
                            #add CDS & set last_start to start of feat.
                            locuses.append(cur_feature)
                            last_start=cur_feature.location.start.position
                        else:
                            print "|---->Intergenic distance too long; ending here"
                            break
                else:
                    print "|---->Gene in reverse strand; ending here"
                    break
        else:
            print "Non-CDS: ", cur_feature.type
    
    #if we did not manage to get something
    if (last_start is None):
        #assuming internet connection problem, wait it out a bit
        #time.sleep(30)  #sleep for 30 seconds
        print "|---->No data retrieved."
        #return that we were not able to grab anything
        return(-2,-2)
        
    #generate positions to grab
    grab_up=last_start - up_d
    grab_dw=last_start + down_d
    
    #define taxonomy as species + first 4 taxonomy tags (if available)
    taxonomy="|"+nuc_rec.annotations["source"]
    ind = 0
    ll=len(nuc_rec.annotations["taxonomy"])
    while ind < 4:
        if ind>=ll:
            taxonomy+="|"+""
        else:
            taxonomy+="|"+nuc_rec.annotations["taxonomy"][ind]
        ind+=1
    
    #pair locus tag and protein ID
    operon_members=""
    for loc in reversed(locuses):
        operon_members+="|"+loc.qualifiers["locus_tag"][0]+"["+\
                        loc.qualifiers["protein_id"][0]+"]"
    #define fasta header
    fasta_header=">"+nuc_rec.id+taxonomy+operon_members
    

    #assign fasta body
    fasta_body=str(nuc_rec.seq[grab_up:grab_dw])
    
    #define fasta field
    fasta_field=fasta_header + "\n" + fasta_body
    
    #pair locus tag and protein ID
    tab_operon_genes=""
    tab_operon_prots=""
    for loc in reversed(locuses):
        tab_operon_genes+="|"+loc.qualifiers["locus_tag"][0]
        tab_operon_prots+="|"+loc.qualifiers["protein_id"][0]
    
    #define tabbed field
    tabulated=nuc_rec.id + "|" + fasta_body \
              +taxonomy+tab_operon_prots
    
    return fasta_field, tabulated
    
            
