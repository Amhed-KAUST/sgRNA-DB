##################All from Scratch
###After multiple intents I realized several errors in gff3:
#Exon and mRNA in WormBase contain the 5'UTR
##Resources: Extracted from
#Caenorhabditis_elegans.WBcel235.95.gff3.gz: Ensemble95 Database
#ce11-appris_data.principal.txt: APPRIS database based on Ensemble 88; grep "PRINCIPAL:1" ce11-appris_data.principal.txt | awk -F"\t" '{print $2}' | sort | uniq | sort > tmp; grep "PRINCIPAL:1" ce11-appris_data.principal.txt | awk -F"\t" '{print $2}' | sort | uniq -d > dupplicated-appris.txt ; comm tmp dupplicated-appris.txt -2 -3 > appris-to-get.txt
#Sorted and combined TSS: awk -F"\t" '{if(array[$1"-"$2"-"$3] != 0){print array[$1"-"$2"-"$3]"\t;"$5 }else{array[$1"-"$2"-"$3]=$0}}' ATAC-TSSs-FWD.BED ATAC-TSSs-REV.BED | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5$7"\t"$6}' - > ATAC-CombTSSs.BED
#EnsePrim-Lists:
#awk -F"\t" '{array[$2]=$0}END{for(i in array){print i"\t"array[i]}}' ce11-appris_data.principal.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - appris-to-get.txt| awk -F"\t" '{print $3";"$4}' > appris-to-getG.txt;
#zcat ../Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="mRNA"){print $0}}}' - | awk -F"\t|;|=|:" '{print $1"\t"$4"\t"$5"\t"$14";"$16"\t.\t"$7}' | awk -F"\t" '{print $4"\t"$0}' | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - appris-to-getG.txt > Ensembl95-Primary.txt;
#awk -F"\t|;" '{print $2}' Ensembl95-Primary.txt > EnsePrim-List.txt;
#zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="mRNA"){print $0}}}' - | awk -F"\t|;|=|:" '{print $1"\t"$4"\t"$5"\t"$14";"$16"\t.\t"$7"\t"$14}' | awk -F"\t" '{print $7"\t"$0}' | awk -F'\t' '{if(array[$1] != 0){array[$1]=array[$1]"lalolanda."$0}else{array[$1]=$0}}END{for(keys in array){print array[keys]}}' - | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - non-PrimGenes.txt | awk -F"lalolanda." '{for(i=1;i<=NF;i++){print $i}}' - | awk -F"\t" '{$1=$8="";print $5"\t"$0}' - > Ensembl95-NON-Primary.txt
#WormBaseANNOT: zcat c_elegans.PRJNA13758.WS268.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="mRNA"){print $0}}}' - | awk -F"\t|;|=|:" '{print $14";"$20}' | sort | uniq | awk -F";" '{print $1"\t"$2}' > WormBaseID-locus.txt
#WS268GenesOperon.BED Done previously by GetOperon LIst
#Obtain CDS GeneBodies
zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=|:" '{print $1"\t"$4"\t"$5"\t"$11"\t.\t"$7}' | awk -F"\t" '{if(array[$4] != 0){if(start[$4] > $2){start[$4]=$2};if(end[$4] < $3){end[$4]=$3};array[$4]=$1"\t"start[$4]"\t"end[$4]"\t"$4"\t"$5"\t"$6}else{start[$4]=$2;end[$4]=$3;array[$4]=$0}}END{for(key in array){print key"\t"array[key]}}' - | sort > GeneBodies.txt
zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=|:" '{print $1"\t"$4"\t"$5"\t"$11"\t.\t"$7}' | awk -F"\t" '{if(array[$4] != 0){if(start[$4] > $2){start[$4]=$2};if(end[$4] < $3){end[$4]=$3};array[$4]=$1"\t"start[$4]"\t"end[$4]"\t"$4"\t"$5"\t"$6}else{start[$4]=$2;end[$4]=$3;array[$4]=$0}}END{for(key in array){print array[key]}}' - | sort > GeneBodies.bed

##Gff3 from Wormbase is 1 based convert to 0 based (bedformat)
awk -F"\t" '{OFS="\t";$3=($3-1);$4=$4;print $0}' GeneBodies.txt | sort -k1,1 -k2,2n > GeneBodies.0.txt
awk -F"\t" '{OFS="\t";$2=($2-1);$3=$3;print $0}' GeneBodies.bed | sort -k1,1 -k2,2n > GeneBodies.0.bed

##Obtain closest and higher ATAC-opening for all Isoforms
bedtools closest -io -id -D "a" -k 2 -a GeneBodies.0.bed -b ATAC-CombTSSs.BED | awk '{if($NF > -600){print $0}}' | awk -F"\t" '{if(array[$4] !=0){ split($11,coso,";"); if($6=="+"){tmp=coso[1]}else{tmp=coso[2]};if(val[$4]<tmp){array[$4]=$0}}else{array[$4]=$0; split($11,coso,";"); if($6=="+"){val[$4]=coso[1]}else{val[$4]=coso[2]}}}END{for(key in array){print array[key]}}' | sort -k1,1 -k2,2n > Isoform-with-closestandhigher-ATAC.txt

#Remove those not found
awk -F"\t" '{if($8 == "-1"){}else{print $0}}' Isoform-with-closestandhigher-ATAC.txt > GoodIsoform-with-closestandhigher-ATAC.txt

#Obtain first set of those with primary form and ATAC
awk -F"\t" '{print $4"\t"$0}' GoodIsoform-with-closestandhigher-ATAC.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - EnsePrim-List.txt > Primary-forms-Promoter.txt
#PutWB AS List
awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print $2"\t"array[$1]}else{array[$1]=$0}}' Primary-forms-Promoter.txt - | awk '{print $1}' - | sort > WBListwithGoodProm.txt
awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print $2"\t"array[$1]}else{array[$1]=$0}}' Primary-forms-Promoter.txt - > WBGoodIsoform-with-closestandhigher-ATAC.txt
#Make subtemporary directory
mkdir almost
cp WBGoodIsoform-with-closestandhigher-ATAC.txt almost/Primary+ATAC.txt
cp GoodIsoform-with-closestandhigher-ATAC.txt Isoform-with-closestandhigher-ATAC.txt



awk '{print $4"\t"$0}' GoodIsoform-with-closestandhigher-ATAC.txt > testes
awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-NON-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]"\t"$0}else{array[$1]=$2}}' - testes >WBGenes-non-prim-but-ATAC.txt

awk '{if(array[$1] != 0){if(val[$1] < $15){array[$1]=$0; val[$1]=$15}}else{array[$1]=$0; val[$1]=$15}}END{for(keys in array){print array[keys]}}' WBGenes-non-prim-but-ATAC.txt > WBGenes-nonP-ATAC-ShortestD.txt
cp WBGenes-nonP-ATAC-ShortestD.txt almost/NonPrimary+ATAC.txt


#Obtain semi Lists
cd almost/
awk '{print $2}' Primary+ATAC.txt | sort > A
sort ../EnsePrim-List.txt > B
comm -2 -3 B A > PrimaryNoAtac-list.txt
cd ..

#Obtain Prim No ATAC
awk -F"\t" '{print $4"\t"$0}' GeneBodies.0.bed | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - almost/PrimaryNoAtac-list.txt | awk -F"\t" '{if($7 =="-"){cen=$4+250}else{cen=$3-250}; print $0"\t"$2"\t"(cen-75)"\t"(cen+75)"\tPutative_Centered250Upstream-150bpRegion\t.\t"$7"\t-175"}' - > PartialPrimNOATAC.txt
awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]"\t"$0}else{array[$1]=$2}}' - PartialPrimNOATAC.txt > almost/Primary-ATAC.txt

#Obtain semi Lists
cd almost/
awk '{print $2}' NonPrimary+ATAC.txt | sort > A
awk -F"\t|;" '{print $2}' ../Ensembl95-NON-Primary.txt | sort > B
comm -2 -3 B A > NONPrimaryNoAtac-list.txt
cd ..

awk -F"\t" '{print $4"\t"$0}' GeneBodies.0.bed | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - almost/PrimaryNoAtac-list.txt | awk -F"\t" '{if($7 =="-"){cen=$4+150}else{cen=$3-150}; print $0"\t"$2"\t"(cen-75)"\t"(cen+75)"\tPutative_Centered150Upstream-150bpRegion\t.\t"$7"\t-75"}' - > PartialPrimNOATAC.txt
awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]"\t"$0}else{array[$1]=$2}}' - PartialPrimNOATAC.txt > almost/Primary-ATAC.txt

awk -F"\t" '{print $4"\t"$0}' GeneBodies.0.bed | awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - almost/NONPrimaryNoAtac-list.txt | awk -F"\t" '{if($7 =="-"){cen=$4+150}else{cen=$3-150}; print $0"\t"$2"\t"(cen-75)"\t"(cen+75)"\tPutative_Centered150Upstream-150bpRegion\t.\t"$7"\t-75"}' - > PartialNONPrimNOATAC.txt

awk -F"\t|;" '{print $2"\t"$1}' Ensembl95-NON-Primary.txt | awk -F"\t" '{if(array[$1] != 0){print array[$1]"\t"$0}else{array[$1]=$2}}' - PartialNONPrimNOATAC.txt | awk -F"\t" '{if(array[$1] != 0){tmp=($5-$4); if(tmp > val[$1]){array[$1]=$0; val[$1]=($5-$4)}}else{array[$1]=$0; val[$1]=($5-$4)}}END{for(keys in array){print array[keys]}}' > almost/NonPrimary-ATAC.txt

cd almost/

rm A
rm B

mkdir tmp
mv PrimaryNoAtac-list.txt tmp/
mv NONPrimaryNoAtac-list.txt tmp/

cat *.txt | awk '{print $1}' | sort | uniq -c | awk '{if($1>1){print $2}}' > Duplicated

awk '{print $1}' NonPrimary-ATAC.txt | sort > A
sort Duplicated > B
comm -2 -3 A B > C

awk -F"\t" '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' NonPrimary-ATAC.txt C > RealNonPrimary-ATAC.txt
mv NonPrimary-ATAC.txt tmp/OldNonPrimary-ATAC.txt
mv RealNonPrimary-ATAC.txt NonPrimary-ATAC.txt
rm A
rm B
rm C
rm Duplicated

for file in `ls *.txt`; do echo $file; awk -F"\t" '{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $file > ${file%.txt}.bed; done
for file in `ls *.txt`; do echo $file; awk -F"\t" '{print $9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' $file > ${file%.txt}.TSS.bed; done
cd ..

##First constitutive Exon (REMEMBER TO CHANGE TO 0 index!!!!!!!)
zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="exon"){print $0}}}' - | grep "constitutive=1" | awk -F"\t|:|;" '{print $10"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' | awk '{if(array[$1] != 0){if($6=="-"){tmp=$4; if(tmp>val[$1]){array[$1]=$0; val[$1]=$4}}else{tmp=$3; if(tmp<val[$1]){array[$1]=$0; val[$1]=$3}}}else{array[$1]=$0; if($6=="-"){val[$1]=$4}else{val[$1]=$3}}}END{for(keys in array){print array[keys]}}' > FirstConstitutiveExon-Ensemble95.txt

##Get CDS
zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|:|;" '{print $1"\t"$4"\t"$5"\t"$10"-CDS\t"$6"\t"$7}' | sort -k1,1 -k2,2n > CDS-all.bed
awk -F"\t" '{OFS="\t";print $2,$3,$4,$1,$5,$6}' FirstConstitutiveExon-Ensemble95.txt | sort -k1,1 -k2,2n > FirstConstitutiveExon-Ensemble95.bed
bedtools intersect -a FirstConstitutiveExon-Ensemble95.bed -b CDS-all.bed | awk '{if(array[$4] != 0){array[$4]=$0}else{array[$4]=$0}}END{for(keys in array){print array[keys]}}' | sort -k1,1 -k2,2n > FirstCDS.bed

awk '{print $1}' GeneBodies.0.txt | sort | uniq > A
awk '{print $4}' FirstCDS.bed | sort | uniq > B
comm -2 -3 A B > List-without-Constitutive.txt

zcat Caenorhabditis_elegans.WBcel235.95.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=|:" '{print $1"\t"$4"\t"$5"\t"$11"\t.\t"$7}' | awk -F"\t" '{if(array[$4] != 0){if(start[$4] > $2){start[$4]=$2; end[$4]=$3};if(end[$4] < $3){end[$4]=$3; start[$4]=$2};array[$4]=$1"\t"start[$4]"\t"end[$4]"\t"$4"\t"$5"\t"$6}else{start[$4]=$2;end[$4]=$3;array[$4]=$0}}END{for(key in array){print key"\t"array[key]}}' - | awk '{if(array[$1] != 0){print array[$1]}else{array[$1]=$0}}' - List-without-Constitutive.txt | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -k1,1 -k2,2n > FirstCDS-IsoformNonConstitutive.bed

#Convert-O indexed
cat FirstCDS.bed FirstCDS-IsoformNonConstitutive.bed | awk -F"\t" '{OFS="\t";$2=($2-1);$3=$3;print $0}' | sort -k1,1 -k2,2n > AlmostConstitutive.bed


cd almost/tmp/

##TOFinish-all-in-loop
for file in `ls ../*.txt`; do
echo $file;
awk -F"\t" '{print $2"\t"$0}' $file > tmp;

###This filter removes dead genes or pseudo genes
awk -F"\t" '{print $4"\t"$0}' ../../AlmostConstitutive.bed | awk -F"\t" '{if(array[$1] != 0){print $2"\t"$0"\t"array[$1]}else{array[$1]=$0}}' - tmp | awk -F"\t" '{if(array[$1] != 0){print array[$1]"\t"$0}else{array[$1]=$2}}' ../../WormBaseID-locus.txt - | awk -F"\t" 'OFS="\t"{print $2,$1,$3,$6,$7,$8,$2";"$1";"$9,$10,$11,$12,$13,$14,$15";WormbaseScore="$16";DistancetoATG="$18,$16,$17,$20,$21,$22,"Exon-"$23,$24,$25}' > ${file%.txt}+CDS.txt;
done

##Run R code
