cd ~/Project/Busco/ 
 
#Make a list of the assembly names
     #assembly names is the first colomn,
     # so the 1st element before ''.'' is the assembly without its version number
assembly=`grep GCA ../metadata.tbl | cut -f1 -d '.'`

#remove the output in case already existed
rm buscoNames-ids.txt

#Go through each assembly name
for i in $assembly
do
    #Get ids of the fasta for this assembly, ie busco names
    #the output contains the file in which it found something:the id
    grep ">" ${i}_busco/run*/busco_sequences/single_copy_busco_sequences/*fna >> buscoNames-ids.txt
done

#The first element before '/' is assembly_busco
#The 5th element is the last element. It is buscoName.fna:>bedBusco
cut -f1,5 -d '/' buscoNames-ids.txt >  buscoNames-ids.txt2

#Remove "_busco" so only the assembly name remain
sed 's/_busco\// /' buscoNames-ids.txt2 -i

#Remove ":>" and replace it by a space
sed 's/:>/ /' buscoNames-ids.txt2 -i
#Remove ".fna"
sed 's/.fna//' buscoNames-ids.txt2 -i

#buscoNames-ids.txt2 thus contains 3 colums: assembly name, busco name, bed of the busco

#buscoNames-ids.txt can be removed
