#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  lib_test.py
#  
#  Copyright 2014 <k7n7vi@ivan-UbuVBox>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

#imports the pn_parse lib containing the functions for GI accession
import pn_parse_lib as PN
#imports getopt to handle the cmd line reading
import getopt
import sys

def main():
    """gets list of protein GIs, acccesses their CDS nucleotide records
       crawls back until getting over max intergenic space or finding
       gene in opposite orientation. 
       Returns list of genes and their shared upstream region according
       to the offsets specified in the input (up and downstream of the
       TLS) 
       writes these upstream regions in the specified output file
    """
    
    #set default parameters
    GI_filename="GIs.txt"
    off_up=250
    off_dw=50
    w_size=5000
    int_dist=50
    out_filename="upstream_regs.fas"
    
    #get cmd parameters
    try:
        opts, args=getopt.getopt(sys.argv[1:],"G:Ws:Id:Pup:Pdw:O:")
    except getopt.GetoptError:
        print 'prot_upstream_grab.py -G <GI file> -Ws <window size> -Id \
               <intergenic distance> -Pup <upstream_offset> \
               -Pdw <downstream_offset> -O <output file>'
        sys.exit(2)
   
    #assign parameters
    for opt, arg in opts:
        if opt == '-G':
            GI_filename=arg
        elif opt == '-Ws':
            w_size=int(arg)
        elif opt == '-Id':
            int_dist=int(arg)
        elif opt == '-Pup':
            off_up=int(arg)
        elif opt == '-Pdw':
            off_dw=int(arg)            
        elif opt == '-O':
            out_filename=arg
        elif opt == '-askme':
            GI_filename = raw_input('Enter the GI file name\n')
            
    #verbose
    print "Using: ", GI_filename, " as input"
    print "Writing to: ", out_filename
    print "Logging errors to: ", out_filename + ".err"
    
    #open files for ouput
    try:
        out_file = open(out_filename,"w")
    except (IOError, OSError) as file_open_exception:
        print "Something went wrong while opening the output file"
        sys.exit()    
    try:
        out_file_tab = open(out_filename+".tab","w")
    except (IOError, OSError) as file_open_exception:
        print "Something went wrong while opening the output file"
        sys.exit()  
        
    #open file for error recording
    try:
        err_file = open(out_filename+".err","w")
    except (IOError, OSError) as file_open_exception:
        print "Something went wrong while opening the error file"
        sys.exit()    
        
    #calls function to read GI list from file
    myGIs = PN.read_GIs(GI_filename)

    #for every GI in read list
    for GI in myGIs:
        #get the associated nucleotide accession (+position)
        nuc_acc=PN.get_nuc_acc(GI)
        #if we can obtain nucleotide accession
        if (nuc_acc !=-1) and (nuc_acc != -2):
            #parse the nucleotide accession string to obtain acc, strand...
            acc, strand, sstart, sstop=PN.parse_nuc_accession(nuc_acc)
            if acc==-1:
                print "Error parsing: ", GI
                err_file.write(GI)
                err_file.write("\t"+"Parsing error (join?)!")
                err_file.write("\n")
            else:
                #grab the upstream (window) sequence object
                up_seq = PN.SeqIO_grab_ext_nuc_rec([acc,strand,sstart,
                                               sstop], w_size)   
                #if there was an error fetching
                if (not up_seq):
                    print "Efetch error: ", GI
                    err_file.write(GI)
                    err_file.write("\t"+"Efetch error!")
                    err_file.write("\n")
                    
                #trace operon and grab upstream sequence
                fasta, tabbed = PN.SeqIO_extract_operon_up(up_seq, off_up, 
                                                       off_dw, int_dist)
            
                #if operon could not be completed for this one
                if (fasta==-1):
                    #write the protein ID to the error file
                    err_file.write(GI)
                    if tabbed == -2:
                        err_file.write("\t"+"Likely misannotated start")
                    err_file.write("\n")
                    print  "Logging error. \n"
                elif (fasta==-2):
                    #write to error file
                    print GI
                    err_file.write(GI)
                    err_file.write("\t"+"No data obtained!")
                    err_file.write("\n")
                    print  "Logging error. \n"
                else:
                    #print fasta_seq
                    out_file.write(fasta+"\n")
                    #print tabbed seq
                    out_file_tab.write(tabbed+"\n")

        #handle access to nucleotide record errors
        #either protein has no nucleotide record, or efetch crashed
        else:
            if nuc_acc==-1:
                print "No nucleotide accession available for: ", GI
                err_file.write(GI)
                err_file.write("\t"+"No nucleotide accession available!")
                err_file.write("\n")
            if nuc_acc==-2:
                print "Efetch failed for: ", GI
                err_file.write(GI)
                err_file.write("\t"+"Failed efetch!")
                err_file.write("\n")
            
    out_file.close()
    out_file_tab.close()
    err_file.close()
    return 0
    
if __name__ == '__main__':
    main()

