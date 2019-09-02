##Deconstructing library
##From 1 to 162 are sgRNAs
##different stiles, separator "_" provides trhee kinds, those with sfu sgRNAs, those with 2mm sgRNAs and those specific made by Mohammed. 
##Idea: convert into trheee temp files to then merge them again...
###ALso, after I have to figureout how to add random stuff and the extra made by Christian... lets check!!!

###Main headers of table
#Well	Number	WBID	Locus	Transcript	Type(location_in_gene)	Chr	Start	End	Strand	Sequence	Promoter	cr(1_or_2)	FwdPrimer	RevPrimer
cp OligoDBv2List.csv OligoDBv2List.original.csv

perl -p -e 's/TSSRegion-/TSSRegionbed-/g' OligoDBv2List.original.csv > OligoDBv2List.csv

awk -F"," '{if($1<163){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | grep "WBGene" | awk -F"_" '{if(NF==12){print $0}}' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk -F">" '{print $1""$2}' | awk -F"\t" '{split($3,coso,";"); split($6,toto,"bed-"); split($9,as,":"); OFS="\t"; split(as[4],bs,"-"); split(bs[2],cs,"("); split(cs[2],ds,")"); split($9,ss,";"); split(ss[2],seq,"="); print $1,$2,coso[1],coso[2],coso[3],toto[1],as[3],bs[1],cs[1],ds[1],seq[2],$10,$11,$13,$14}' > tempA

awk -F"," '{if($1<163){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | grep "WBGene" | awk -F"_" '{if(NF==9){print $0}}' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F">" '{print $1""$2}' | awk -F"\t" '{split($3,coso,";"); split($6,toto,"bed-"); split($6,as,":"); OFS="\t"; split(as[6],bs,"-"); split(bs[2],cs,"("); split(cs[2],ds,")"); split(toto[2],seq,";"); print $1,$2,coso[1],coso[2],coso[3],toto[1],as[5],bs[1],cs[1],ds[1],seq[1],$7,$8,$10,$11}' > tempB

awk -F"," '{if($1<163){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | grep "WBGene" | awk -F"_" '{if(NF==7){print $0}}' | perl -p -e 's/\t_/_/g' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7}' | awk -F">" '{print $1""$2}' |  awk -F"\t" '{split($3,coso,";"); split($4,toto,"-"); OFS="\t"; print $1,$2,coso[1],coso[2],coso[3],toto[1]"Design","-","-","-","NULL","NULL",$5,$6,$8,$9}' | perl -p -e 's/CeN50-1/CeN50-/g' | perl -p -e 's/CeN50-/CeN50-1/g' > tempC

awk -F"," '{if($1<163){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | grep -v "WBGene" | awk -F">" '{print $1""$2}' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8}' | awk -F"\t" '{OFS="\t"; if($5=="20M"){bin="RandomBinding"}else{bin="RandomNotBinding"}print $1,$2,"NULL","NULL","NULL",bin,"NULL","NULL","NULL","NULL",$3"-"$4,$6,$7,$9,$10}' > tempD

cat tempA tempB tempC tempD | sort -n -k1 -k2,2 > sgRNAs1-162.tab

awk -F"," '{if($1<187){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | grep -v "WBGene" | awk -F">" '{print $1""$2}' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8}' | awk -F"\t" '{OFS="\t"; if($5=="20M"){bin="RandomBinding"}else{bin="RandomNotBinding"}print $1,$2,"NULL","NULL","NULL",bin,"NULL","NULL","NULL","NULL",$3"-"$4,$6,$7,$9,$10}' > tempD

cat tempA tempB tempC tempD | sort -n -k1 -k2,2 > sgRNAs1-186.tab

awk -F"," '{if(($1<191)&&($1>186)){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv |  awk -F">" '{print $1""$2}' |  awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5}' | awk -F"\t" '{OFS="\t"; split($3,coso,";"); split($4,toto,"bed"); print $1,$2,coso[1],coso[2],coso[3],toto[1],"NULL","NULL","NULL","NULL",$4,$5,"oldcr",$6,$7}' > tempE

awk -F"," '{if(($1==191)){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | awk -F">" '{print $1""$2}' | perl -p -e 's/coding_exon/coding-exon/g' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F"\t" '{split($5,papa,":"); if($6 == "+"){start=papa[2]-16; end=papa[2]+3 } else { start=papa[2]-2; end=papa[2]+17}; OFS="\t"; print $1,$2,"WBGene00006843","unc-119","M142.1a","Expanded",papa[1],start-1,end,$6,"NULL",$8,"oldcr",$10,$11}' > tempF

awk -F"," '{if(($1==192)){print $1"\t"$2"\t"$5}}' OligoDBv2List.csv | awk -F">" '{print $1""$2}' | perl -p -e 's/standard_GFP/standard-GFP/g' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F"\t" '{OFS="\t"; print $1,$2,$4,$4,$4,$3,$5,$6,$6,$6,"NULL",$8,"oldcr",$10,$11}' > tempG


cat tempA tempB tempC tempD tempE tempF tempG| sort -n -k1 -k2,2 > sgRNAs1-192.tab
