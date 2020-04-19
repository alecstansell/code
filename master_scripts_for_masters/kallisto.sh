#Kallisto index generation and quantification

#Get index file from BiomaRt

#Generate transcript debruijn index file
kallisto index -i transcripts.idx transcripts.fasta.gz


#Quantify samples

#rather use a loop here but - due to computer constantly crashing I ran them line by line

kallisto quant -i index.idx -o out_SRR3680435 -b 100 SRR3680435_1.fastq.gz SRR3680435_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680442 -b 100 SRR3680442_1.fastq.gz SRR3680442_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680443 -b 100 SRR3680443_1.fastq.gz SRR3680443_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680444 -b 100 SRR3680444_1.fastq.gz SRR3680444_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680445 -b 100 SRR3680445_1.fastq.gz SRR3680445_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680446 -b 100 SRR3680446_1.fastq.gz SRR3680446_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680447 -b 100 SRR3680447_1.fastq.gz SRR3680447_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680448 -b 100 SRR3680448_1.fastq.gz SRR3680448_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680449 -b 100 SRR3680449_1.fastq.gz SRR3680449_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680450 -b 100 SRR3680450_1.fastq.gz SRR3680450_2.fastq.gz 
kallisto quant -i index.idx -o out_SRR3680451 -b 100 SRR3680451_1.fastq.gz SRR3680451_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680452 -b 100 SRR3680452_1.fastq.gz SRR3680452_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680453 -b 100 SRR3680453_1.fastq.gz SRR3680453_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680454 -b 100 SRR3680454_1.fastq.gz SRR3680454_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680455 -b 100 SRR3680455_1.fastq.gz SRR3680455_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680456 -b 100 SRR3680456_1.fastq.gz SRR3680456_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680457 -b 100 SRR3680457_1.fastq.gz SRR3680457_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680458 -b 100 SRR3680458_1.fastq.gz SRR3680458_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680459 -b 100 SRR3680459_1.fastq.gz SRR3680459_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680460 -b 100 SRR3680460_1.fastq.gz SRR3680460_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680461 -b 100 SRR3680461_1.fastq.gz SRR3680461_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680462 -b 100 SRR3680462_1.fastq.gz SRR3680462_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680463 -b 100 SRR3680463_1.fastq.gz SRR3680463_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680464 -b 100 SRR3680464_1.fastq.gz SRR3680464_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680465 -b 100 SRR3680465_1.fastq.gz SRR3680465_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680466 -b 100 SRR3680466_1.fastq.gz SRR3680466_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680467 -b 100 SRR3680467_1.fastq.gz SRR3680467_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680468 -b 100 SRR3680468_1.fastq.gz SRR3680468_2.fastq.gz
kallisto quant -i index.idx -o out_SRR3680469 -b 100 SRR3680469_1.fastq.gz SRR3680469_2.fastq.gz

