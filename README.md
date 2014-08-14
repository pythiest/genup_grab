This is an extension of the PN_parse repo. The PN_parse repo contains a library to identify and grab the CDS for a set of protein GI accession numbers. The GenUp_grab essentially builds upon PN_parse to identify the CDS, and then crawl back upstream until it finds the first, non-intraoperon intergenic region and grab that.

Schematically, the process is as follows:
- For all GIs in file
  - Identify CDS in NCBI
  - Move upstream of CDS, until (new CDS end is farther than X bp from current CDS start) OR (new CDS is in reverse orientation)
  - Report all CDSs (i.e. the operon) identified within the FASTA header for the intergenic sequence grabbed
