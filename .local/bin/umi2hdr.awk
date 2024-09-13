#!/bin/awk -f

###################################################################################################################
# Take as input a FASTQ file with interleaved read-pairs.
# Check successive header names for similarity, and output to stderr read1 
# with read2 sequence in its header line
###################################################################################################################


{
  r[NR%8] = $0 # Put each read pair in its slot
}
NR%8==0{ # Once 8 lines of the fastq where read
  split(r[1], hdr1, " ")
  split(r[5], hdr2, " ")
  if (hdr1[1]!=hdr2[1]) { # compare sequence id
    print "Error at line " NR "header names do not match" > "/dev/stderr";
    exit 1;
  } else {
    print(hdr1[1] "_" r[6]) # Add read2 sequence to header line
    print(r[2]) # read1 sequence
    print(r[3]) # + seperator
    print(r[4]) # read1 quality
  }
}


