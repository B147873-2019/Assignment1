#SOME PREPARATION

cd
#change to home directory
mkdir assignment
#build a mew directory 
cd assignment
#change directory again




#QUESTION1
#PROCESSING DATA WITH FASTQC



mkdir fastqc_output
#build a mew directory called xxxxxx

#run fastqc
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/216_L8_1.fq.gz 

echo "216_L8_1 finished"
#then print to tell us that fastqc successfully finished

#run other files too
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/216_L8_2.fq.gz
echo "216_L8_2 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/218_L8_1.fq.gz
echo "218_L8_1 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/218_L8_2.fq.gz
echo "218_L8_2 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/219_L8_1.fq.gz
echo "219_L8_1 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/219_L8_2.fq.gz
echo "219_L8_2 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/220_L8_1.fq.gz
echo "220_L8_1 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/220_L8_2.fq.gz
echo "220_L8_2 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/221_L8_1.fq.gz
echo "221_L8_1 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/221_L8_2.fq.gz
echo "221_L8_2 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/222_L8_1.fq.gz
echo "222_L8_1 finished"
fastqc -o fastqc_output /localdisk/data/BPSM/Assignment1/fastq/222_L8_2.fq.gz
echo "222_L8_2 finished"
#all the sequences files are done 

cd fastqc_output

ls -al
#show us all the quality check results

cd -
#go back to last repository




#QUESTION2
#the analyse of reports can be seen in my PDF file




#QUESTION3
#USE BOWTIE2 TO ALIGN SEQUENCES


#some preparation

mkdir bowtie2_output

cd bowtie2_output

#unzip files and build template before running bowtie2

cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz .
#copy the file to this repository

gunzip Tb927_genome.fasta.gz
#unzip the file

#build index of template
bowtie2-build Tb927_genome.fasta index

#show us all the builded data(6)
ls -al

#run bowtie2

bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/216_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/216_L8_2.fq.gz -S 216.sam
#the sequencing method we used is "paired-end sequencing"
#so "-x" means to find the reference by search the index named"index"
#while  "-1" and "-2" mean the directory of two DNA fragments
#"-S" means the output file name would be "216.sam"

echo "216group aligned"
#tells us that 216 is done aligning

#do the same to other files
bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/218_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/218_L8_2.fq.gz -S 218.sam
echo "218group aligned"
bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/219_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/219_L8_2.fq.gz -S 219.sam
echo "219group aligned"
bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/220_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/220_L8_2.fq.gz -S 220.sam
echo "220group aligned"
bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/221_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/221_L8_2.fq.gz -S 221.sam
echo "221group aligned"
bowtie2 -x index -1 /localdisk/data/BPSM/Assignment1/fastq/222_L8_1.fq.gz -2 /localdisk/data/BPSM/Assignment1/fastq/222_L8_2.fq.gz -S 222.sam
echo "222group aligned"
#all the groups of sequences is finished aligning

ls -al
#show us all the results

cd -
#go back to last repository



#QUESTION4

#some prepration

mkdir Tbbgenes_bed

cd Tbbgenes_bed
#make new repository and go in to it

#copy files 
cp /localdisk/data/BPSM/Assignment1/Tbbgenes.bed .
cp /localdisk/home/s1919406/assignment/bowtie2_output/216.sam .
cp /localdisk/home/s1919406/assignment/bowtie2_output/218.sam .
cp /localdisk/home/s1919406/assignment/bowtie2_output/219.sam .
cp /localdisk/home/s1919406/assignment/bowtie2_output/220.sam .
cp /localdisk/home/s1919406/assignment/bowtie2_output/221.sam .
cp /localdisk/home/s1919406/assignment/bowtie2_output/222.sam .

#other preparations

#change file format

samtools view -bS 216.sam > 216.bam 
#view is used to change format
#"b" means to output files in ".bam" format
#"S" means the input files are in "sam" format
#sort files

samtools sort 216.bam -o 216
#sort means to sort alignments in a proper order by leftmost column
#we have to sort the file before indexing
#build bai index

samtools index -b 216
#"-b" means to create a "bai" index (which is the common format with no special format options)
#we have to build the index before fast random access

#do the same to all sam files
samtools view -bS 218.sam > 218.bam
samtools view -bS 219.sam > 219.bam
samtools view -bS 220.sam > 220.bam
samtools view -bS 221.sam > 221.bam
samtools view -bS 222.sam > 222.bam  
samtools sort 218.bam -o 218
samtools sort 219.bam -o 219
samtools sort 220.bam -o 220
samtools sort 221.bam -o 221
samtools sort 222.bam -o 222
samtools index -b 218
samtools index -b 219
samtools index -b 220
samtools index -b 221
samtools index -b 222
   
#sort template
sort -k1,1 -k2,2n Tbbgenes.bed > sorted
#"-k1,1 -k2,2n" means to sort the file by first column, then output the file by numeric order of column 2
#the first colume is "roughly grouped gene name" and the second colume is the start position of gene

#run bedtools
#use bedtools multicov to count number of reads

bedtools multicov -bams 216 218 219 -bed sorted | grep -w "gene" | grep "+" | grep -v "RNA" > slender
#"bams" means the name of input bam file is "216"
#"bed" means the name of input template bed file is "sorted"
#"grep -w "gene"" means to find in file ahead the word gene, and output the whole line(the gene name contains pseudogene which has uncertain usage, and probably need to be filtered because most of them don't transcript)
#"grep "+"" means to filter antisense("-") which is RNA that regulates gene expression(except "code for genes" in the topic)
#"grep -v "RNA""means to filter "RNA"(because some gene names include "RNA" in it and we only need real DNA)
#"> slender" means to pipe the output to a new file named slender

bedtools multicov -bams 220 221 222 -bed sorted | grep -w "gene" | grep "+" | grep -v "RNA" > stumpy

bedtools multicov -bams 216 218 219 220 221 222 -bed sorted | grep -w "gene" | grep "+" | grep -v "RNA" > grouped

#count total number of reads per group
cat slender | awk 'BEGIN{a=0};{a += $7}END{print a}' > number_of_reads_216_slender
cat slender | awk 'BEGIN{a=0};{a += $8}END{print a}' > number_of_reads_218_slender
cat slender | awk 'BEGIN{a=0};{a += $9}END{print a}' > number_of_reads_219_slender
cat stumpy | awk 'BEGIN{a=0};{a += $7}END{print a}' > number_of_reads_220_stumpy
cat stumpy | awk 'BEGIN{a=0};{a += $8}END{print a}' > number_of_reads_221_stumpy
cat stumpy | awk 'BEGIN{a=0};{a += $9}END{print a}' > number_of_reads_222_stumpy
#BEGIN and END is like loop in awk,and "+=" means to add up the whole column

#finally print those numbers and pipe the file to "number_of_reads_per_group"





#QUESTION5

#use bedtools groupby to count mean of reads by gene name
bedtools groupby -i grouped -g 1 -c 7,8,9,10,11,12 -o mean > grouped_meaned
#groupby is perfect for summary data in several column(column 7,8,9,10,11,12) by another grouped column(column 1) from number_of_reads_per_group(input file), and count the mean value for each group("-o mean", the "-o" means option), then pipe the data to "grouped_meaned"

#count reads per group
cat grouped_meaned | 
while read name one two three four five six;
 do all_slender=$(echo "$one + $two + $three"|bc);
 mean_slender=$(echo "scale=5;$all_slender / 3"|bc);
 all_stumpy=$(echo "$four + $five + $six"|bc);
 mean_stumpy=$(echo "scale=5;$all_stumpy / 3"|bc);
 echo -e "$name\t$mean_stumpy\t$mean_slender";
 done

