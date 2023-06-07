for i in $(ls *.bam | sed 's/.bam//')
do
samtools coverage -o ./coverage/$i.tsv $i\.bam
done

for i in $(ls *.bam | sed 's/.bam//')
do
samtools flagstat -@ 32 -O tsv $i\.bam > ./flagstats/$i.tsv
done

cd ./coverage/

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $6}' $i.tsv > $i.tmp
done
paste *tmp > coverage.txt
rm -rf *tmp

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $7}' $i.tsv > $i.tmp
done
paste *tmp > meandepth.txt
rm -rf *tmp

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $8}' $i.tsv > $i.tmp
done
paste *tmp > baseQ.txt
rm -rf *tmp

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $9}' $i.tsv > $i.tmp
done
paste *tmp > mapQ.txt
rm -rf *tmp

cd ./../flagstats/

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $1}' $i.tsv > $i.tmp
done
paste *tmp > flagstats.txt
rm -rf *tmp
