cat reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa | grep ">" | awk -F" " '{print $1"\t"$7}' | sed -e 's/>//g' -e 's/gene_symbol://g'
