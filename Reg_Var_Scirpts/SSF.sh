# This is a sample PBS script. It will request 1 processor on 1 node
# for 120 hours.
#   
#   
#PBS -l nodes=1:ppn=1
#PBS -l walltime=96:00:00
#  
##PBS -q lionxc-makova
#PBS -j oe
#
#   
#
cd $PBS_O_WORKDIR
#
#
echo " "
echo " "
echo "Job started on `hostname` at `date`"
echo " Snoope Single microsat in read Fetch flanking"
##python microsat_snoop.py ${INPUT}.sangerfq --fastq --period=1 --partialmotifs --minlength=3 --prefix=20 --suffix=20 --hamming=0 --multipleruns  >${INPUT}.mono.out
##python microsat_snoop.py ${INPUT}.sangerfq --fastq --period=2 --partialmotifs --minlength=6 --prefix=20 --suffix=20 --hamming=0 --multipleruns  >${INPUT}.di.out
##python microsat_snoop.py ${INPUT}.sangerfq --fastq --period=3 --partialmotifs --minlength=9 --prefix=20 --suffix=20 --hamming=0 --multipleruns  >${INPUT}.tri.out
##python microsat_snoop.py ${INPUT}.sangerfq --fastq --period=4 --partialmotifs --minlength=12 --prefix=20 --suffix=20 --hamming=0 --multipleruns  >${INPUT}.tetra.out
##echo "finish snooping"
##wc -l ${INPUT}.mono.out
##wc -l ${INPUT}.di.out
##wc -l ${INPUT}.tri.out
##wc -l ${INPUT}.tetra.out
##echo " "
##echo "add name column"
##python addreadname_c.py ${INPUT}.mono.out > ${INPUT}.mono.out1
##python addreadname_c.py ${INPUT}.di.out > ${INPUT}.di.out1
##python addreadname_c.py ${INPUT}.tri.out > ${INPUT}.tri.out1
##python addreadname_c.py ${INPUT}.tetra.out > ${INPUT}.tetra.out1
##echo " "
##echo "sort by readname"
##cat ${INPUT}.mono.out1 |  sort -k 9,9 > ${INPUT}.mono.out2
##cat ${INPUT}.di.out1 |  sort -k 9,9 > ${INPUT}.di.out2
##cat ${INPUT}.tri.out1 |  sort -k 9,9 > ${INPUT}.tri.out2
##cat ${INPUT}.tetra.out1 |  sort -k 9,9 > ${INPUT}.tetra.out2
##
##echo " "
##echo "remove crowded"
##python remove_crowd.py ${INPUT}.mono.out2 > ${INPUT}.mono.new
##python remove_crowd.py ${INPUT}.di.out2 > ${INPUT}.di.new 
##python remove_crowd.py ${INPUT}.tri.out2 > ${INPUT}.tri.new 
##python remove_crowd.py ${INPUT}.tetra.out2 > ${INPUT}.tetra.new 
##
##echo "finish remove crowded"
##rm ${INPUT}.mono.out
##rm ${INPUT}.di.out
##rm ${INPUT}.tri.out
##rm ${INPUT}.tetra.out
##rm ${INPUT}.mono.out1
##rm ${INPUT}.di.out1
##rm ${INPUT}.tri.out1
##rm ${INPUT}.tetra.out1
##rm ${INPUT}.mono.out2
##rm ${INPUT}.di.out2
##rm ${INPUT}.tri.out2
##rm ${INPUT}.tetra.out2
##
##wc -l ${INPUT}.mono.new
##wc -l ${INPUT}.di.new
##wc -l ${INPUT}.tri.new
##wc -l ${INPUT}.tetra.new
##echo " "
##python pair_fetch_DNA_ff.py ${INPUT}.mono.new
##python pair_fetch_DNA_ff.py ${INPUT}.di.new
##python pair_fetch_DNA_ff.py ${INPUT}.tri.new
##python pair_fetch_DNA_ff.py ${INPUT}.tetra.new
##echo "finish fetch flanking"
##wc -l ${INPUT}.mono_ff_L.txt
##wc -l ${INPUT}.di_ff_L.txt
##wc -l ${INPUT}.tri_ff_L.txt
##wc -l ${INPUT}.tetra_ff_L.txt
echo " "
echo "Job Ended at `date`"
echo " "
