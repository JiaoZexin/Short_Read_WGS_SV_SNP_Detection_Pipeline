grep -oE 'ENSSMAG.{11}' ~/data/turbot/snpeff/22537_gtf/rank_HIGH_LOW_MODERATE_MODIFIER/ann_info.txt | sort | uniq >> 12952.22537.gene.txt
grep -oE 'ENSSMAG.{11}' ~/data/turbot/snpeff/8717_new_gtf/rank_HIGH_LOW_MODERATE_MODIFIER/ann_info.txt | sort | uniq >> 7339.8717.gene.txt
ln -s ./12952.22537.gene.txt background_gene.txt 
ln -s ./7339.8717.gene.txt background_gene.txt
for i in $(cat background_gene.txt); do if grep -w -q "${i}" ../KOannotation.tsv; then echo ${i} >> 1.txt ; else echo no >> 1.txt ; fi; done
for i in $(cat background_gene.txt); do if grep -w -q "${i}" ../GOannotation.tsv; then echo ${i} >> 2.txt ; else echo no >> 2.txt ; fi; done
cp 1.txt KO.gene.txt
sed -i '/no/d' KO.gene.txt 
cp 2.txt GO.gene.txt
sed -i '/no/d' GO.gene.txt 
for i in $(cat KO.gene.txt); do grep ${i} ../KOannotation.tsv >> 8717_SVs_KOannotation.tsv; done 
for i in $(cat GO.gene.txt); do grep ${i} ../GOannotation.tsv >> 8717_SVs_GOannotation.tsv; done
for i in $(cat KO.gene.txt); do grep ${i} ../KOannotation.tsv >> 22537_SVs_KOannotation.tsv; done 
for i in $(cat GO.gene.txt); do grep ${i} ../GOannotation.tsv >> 22537_SVs_GOannotation.tsv; done