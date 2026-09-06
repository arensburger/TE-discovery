Sept 6, 2026 Downloaded known transpoases using the commands below:
```
esearch -db protein -query "transposase[Title] AND txid6656[Organism:exp]" | efetch -format gp > transposases.gp
xtract -input transposases.gp -pattern INSDSeq -element INSDSeq_locus INSDSeq_definition INSDSeq_sequence > all_tnps.txt
cd-hit -i all_tnps.fa -o alltnps_unique.fa
bash rs.sh all_tnps.txt  > all_tnps.fa
```

```
peter@carmen:~/Desktop$ more rs.sh
#!/bin/bash
# usage: ./tsv_to_fasta.sh input.tsv > output.fasta

input="$1"

if [[ -z "$input" ]]; then
    echo "Usage: $0 <input.tsv>" >&2
    exit 1
fi

awk -F'\t' '
{
    accession = $1
    description = $2
    sequence = $3
    print ">" accession " " description
    # wrap sequence at 70 characters per line (standard FASTA width)
    seqlen = length(sequence)
    for (i = 1; i <= seqlen; i += 70) {
        print substr(sequence, i, 70)
    }
}
' "$input"
```

