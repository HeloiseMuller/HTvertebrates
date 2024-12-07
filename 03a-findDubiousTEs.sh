#To get diamon:
conda activate /opt/anaconda3/envs/Busco5/

#Dowload nr = the non-redundant “nr” protein database of NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

#blastx nr against my TE consensus
diamond makedb --in nr.gz -d nr -p 10
diamond blastx -p 20 --more-sensitive -q  DB_clustered/TE_cat_all.long.clustered_rep_seq.fasta -d nr.dmnd --quiet -k 20 -f 6
        "qseqid" "sseqid" "pident" "length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" "stitle" -o  DB_clustered/findDubious/TTE_cat_all.long.clustered_rep_seq_on_nr.out

    #NOTE: if following error "Error reading input stream at line 106386: Invalid character (u) in sequence",
    #One can change special nucleotides to N :
    #sed '/>/!s/u\|o/N/g' DB_clustered/TE_cat_all.long.clustered_rep_seq.fasta > DB_clustered/findDubious/TE_cat_all.long.clustered_rep_seq_noUO.fasta
    
#blastx RepeatPep against my TE consensus
diamond makedb --in RepeatPeps.lib -d repeatPeps
diamond blastx -p 20 --more-sensitive -q DB_clustered/TE_cat_all.long.clustered_rep_seq.fasta -d repeatPeps.dmnd --quiet -k 20 -f 6 "qseqid" "sseqid" "pident" "length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" "stitle" -o DB_clustered/findDubious/TE_cat_all.long.clustered_rep_seq_on_RepeatPeps.out