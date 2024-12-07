#script to change ID of fasta
#$1: correspondance between ID and copy names
#$2 input
#$3: output

awk -F'\t' '
    NR==FNR { map[">"$1] = ">"$2; next }
    $0 in map { $0 = map[$0] }
    { print }
' $1 $2 > $3
