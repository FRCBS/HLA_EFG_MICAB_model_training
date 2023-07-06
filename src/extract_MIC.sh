# Extraction of MICAB region reads using Bazam; https://github.com/ssadedin/bazam
# download bazam.jar:
wget https://github.com/ssadedin/bazam/releases/download/1.0.1/bazam.jar
chmod +x /path/to/bazam/bazam.jar
alias bazam="java -jar /path/to/bazam/bazam.jar"
# or add to ~/.bash_aliases: alias bazam="java -jar /path/to/bazam/bazam.jar"

# MICA (GRCh38): 6:31,399,784-31,415,315 
cd ../1000Genomes_MHC/TSI
mkdir FASTQ/MICA
for BAM in `ls -1 BAM/*bam`; do
  NAME=`echo $BAM | sed 's/BAM\///g'`
  bazam -bam $BAM -L chr6:31399684-31415415 -o FASTQ/MICA/$NAME.fastq
done

# MICB (GRCh38): 6:31,494,881-31,511,124 
cd ../1000Genomes_MHC/TSI
mkdir FASTQ/MICB
for BAM in `ls -1 BAM/*bam`; do
  NAME=`echo $BAM | sed 's/BAM\///g'`
  bazam -bam $BAM -L chr6:31494781-31511224 -o FASTQ/MICB/$NAME.fastq
done

