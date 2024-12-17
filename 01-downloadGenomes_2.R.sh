# First, download summaries about all available genomes
 # Mammals are there: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt
 # Vertebrates other than mammals are there:  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/assembly_summary.txt
 # Invertebrates are there ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/assembly_summary.txt

# Concatenate all three summaries
cat assembly_summary.txt assembly_summary.txt1 assembly_summary.txt2 > allvertebrates_invertebrates_19012022_assembly_summary.txt 

# Get list of the 247 assemblies we selected
cut -f1 -d '.'metadata.tbl | tail -n +2 > tmp_assemblies #don't keep header
grep -F -f tmp_assemblies allvertebrates_invertebrates_19012022_assembly_summary.txt >  selected_species_assembly_summary.txt
rm tmp_assemblies
cut -f1 selected_species_assembly_summary.txt > selected_species.lst

# Download genomes one by one
while read i
do
    echo $i
        #find the ftp link corresponding to the species & replace the 1st ftp by http.
        #otherwise, wget does not work correctly (now work only with ftp...)
    link=`grep $i selected_species_assembly_summary.txt | cut -f20 | sed 's/ftp/http/'`
   
    #the last element of the link is the begining of the name of the file, we need it for wget (wildcard* does not work)
    file=`echo $link | rev | cut -d'/' -f 1 | rev`

    #download fasta only if does not exist
    if [ ! -f genomes/${file}_genomic.fna ] ; then
        wget  ${link}/${file}_genomic.fna.gz -P genomes/ --no-check-certificate
        #decompress the fasta and not keep the compressed one (by default)
        gzip -d genomes/${file}_genomic.fna.gz 
    fi
done <selected_species.lst