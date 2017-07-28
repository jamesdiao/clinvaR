#!/bin/bash

wget -O- http://www.uniprot.org/docs/pkinfam.txt | egrep -o "\w*_HUMAN" | sed 's/_HUMAN//g' | src/update_symbols.py --hgnc gene_with_protein_product.txt > lists/kinases.tsv