# PEmbryo
scripts related to PEmbryo work

Below are examples of command line arguments to run "PEmbryo_align_" .py scripts to evaluate editing efficiency on collected samples/datasets based on the specific target site and target edit. The scripts pull edit details and relevant filenames from 'PEmbryo_sample_sheet_master.xlsx'. Fastq files for the relevant samples, available through SRA, must be stored in a directory called "input_libraries/" or the path updated in the script. 

#tubb5;G>A
#dnmt1;G>T
#dnmt1;G>C
#rnf2;C>G
#chd2;G>A
#chd2;G>C
#hoxd13_1;G>T
#hoxd13_1;G>C
#ctnnb1;+6GtoA
#crygc;+1Gdel
#tspan2;+6GtoC
#col12a1;+2AtoC
#ar2;+1GtoT
#GFP;CT>GC
#dnmt1;GGC>TAG
#dnmt1;GGCAT>TAGTG
#hoxd13_2;G>T
#hoxd13_2;+1nt
#hoxd13_2;+3nt
#hoxd13_2;+8nt
#hoxd13_2;+17nt

Outputs generated by the "PEmbryo_align_" .py scripts (.fastq files containing reads with a determined outcome designation i.e. wild-type, precisely edited, or containing errors) are stored in the directory "output_libraries/" with subdirectories for "wt_reads/", "pe_reads/", and "error_reads/". A .log file is also generated to summarize the results. 
