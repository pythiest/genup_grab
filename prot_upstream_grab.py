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
       and returns their upstream region according to the offsets
       specified in the input (up and downstream of the TLS)
       writes these upstream regions in the specified output file
    """
    
    #set default parameters
    GI_filename="GIs.txt"
    off_up=250
    off_dw=50
    out_filename="upstream_regs.fas"
    
    #get cmd parameters
    try:
        opts, args=getopt.getopt(sys.argv[1:],"G:Pup:Pdw:O:")
    except getopt.GetoptError:
        print 'prot_upstream_grab.py -G <GI file> -Pup <upstream_offset> -Pdw <downstream_offset> -O <output file>'
        sys.exit(2)
   
    #assign parameters
    for opt, arg in opts:
        if opt == '-G':
            GI_filename=arg
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
    
    #open file for ouput
    try:
        out_file = open(out_filename,"w")
    except (IOError, OSError) as file_open_exception:
        print "Something went wrong while opening the output file"
        sys.exit()    

    #calls function to read GI list from file
    myGIs = PN.read_GIs(GI_filename)

    #for every GI in read list
    for GI in myGIs:
        #get the associated nucleotide accession (+position)
        nuc_acc=PN.get_nuc_acc(GI)
        #if we can obtain nucleotide accession
        if nuc_acc:
            #parse the nucleotide accession string to obtain acc, strand...
            acc, strand, sstart, sstop=PN.parse_nuc_accession(nuc_acc)
            #define new grabbing coordinates for upstream region
            #depending on which strand we are operating on
            if (strand=="1"):
                newstart=int(sstart)-off_up
                newstop=int(sstart)+off_dw
            else:
                newstop=int(sstop)+off_up
                newstart=int(sstop)-off_dw
            #grab the corresponding sequence in fasta format
            fasta_seq=PN.get_nuc_rec([acc,strand,str(newstart),str(newstop)])
            #print fasta_seq
            out_file.write(fasta_seq)
            out_file.write("\n")
            
        else:
            print "No nucleotide accession available for: ", GI
        
    
 
    return 0
    
if __name__ == '__main__':
    main()

