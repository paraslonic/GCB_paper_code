gzip -d *.gz
python3.5 geneGraph/source/orthofinder_parse.py -i Orthogroups.txt -o graph
python3.5 geneGraph/source/start_computing.py -i graph.sif -o LF82_B2_graph --reference 50 --save_db graph.db
grep rRNA lf82_ncbi.gb | grep "\.\." | sed -e 's/complement(//' -e 's/)//' | sed -e 's/.\+\s//' -e 's/\.\./\t/' > lf82.rrna
grep "tRNA" lf82_ncbi.gb | grep "\.\." | sed -e 's/complement(//' -e 's/)//' -e 's/.\+\s//' -e 's/\.\./\t/' > lf82.trna
Rscript plot.R
