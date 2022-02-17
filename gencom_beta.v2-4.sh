#bin/bash



#Copyright Alexander A. Martinez C & Gorgas Memorial Institute for Health Studies
#Written by: Alexander A. Martinez C, Genomics and proteomics research unit, Gorgas memorial #Institute For Health Studies.
#Licensed under the Apache License, Version 2.0 (the "License"); you may not use
#this work except in compliance with the License. You may obtain a copy of the
#License at:
#http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software distributed
#under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
#CONDITIONS OF ANY KIND, either express or implied. See the License for the
#specific language governing permissions and limitations under the License.

#This software uses diferent smart programs prepared and collected elsewhere each one retain its particular licence, please comply with them.



#This script takes all the fastq files within a folder and according to the pathogen chose #clean the reads, pileup them and generate several outputs. Currently the pathogen with more #application and outputs availables is SARS-CoV-2

##Start of declaration for positionals arguments.
############ Abreviation of the microorganism to be analyzed ################ 
test="$1"
#read_1="$2"
#read_2="$3"
output_dir="$2"
#echo $test

funct1="$3"

##End of declaration for positionals arguments.

#source /home/drakkar4/anaconda3/etc/profile.d/conda.sh
start=$(date +%s.%N)

#conda  activate genover2 

if  [ -z $test ] 
then
	echo "Virus not specified, for consensus generation use the following ./gencom.sh Virus_abreviature read_1 read_2 consensus_percentage"
	exit 1
fi

################################################################
################ZIKA Virus Analysis#############################

if [ $test = ZIKV ]; then

	echo "starting Zika Virus consensus generation"
	cp /home/drakkar4/data/ZIKV.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/ZIKV.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip -f $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
	#echo $c
	#echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Zika Virus polyprotein coverage" $output_dir/analysis/$d/$d.tsv ZIKV
	echo "finished $d Zika Virus Analysis"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Zika Virus polyprotein coverage" $output_dir/analysis/coverage/depth_file.tsv ZIKV
	

################################################################
################SARS-CoV-2 Analysis#############################


elif [ $test = SARS2 ]; then

	now1="$(date +%c)"
	echo "#####################SARS-CoV-2 pipeline started at $now1#######################"
	echo "################### Reviewing and unzipping files #########################"    
	for read_1 in $output_dir/*fastq  
    do
    #echo ${read_1/.gz/}
    #echo $fastq
    if [[ -r $read_1 ]] ; 
    then
    echo "Using the following files for analysis"
    echo "using $read_1 reads" 
    else
    echo "Not fastq file provided looking for fastq.gz files"
    for fastqgz in $output_dir/*.fastq.gz 
    do
    if [[ -r $fastqgz ]] ; 
    then
    echo "fastq.gz files found unziping them"
    gunzip -f $output_dir/*.gz
    echo "unziping done!!!"    
    else
    echo "Neither gz or fastq files provided for analysis"
    echo "############################################"
    fi
    done
    fi
	done

	echo "                                                                           " 
	echo "                                                                           "     
	echo "############################################"
	echo "starting SARS-Cov-2019 Virus consensus generation"
	cp /home/jovyan/shared/data/SARS2.fas  /opt/conda/lib/python3.9/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/jovyan/shared/data/SARS2.bed /opt/conda/lib/python3.9/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/jovyan/shared/data/SARS2_mutation_db.tsv  /opt/conda/lib/python*/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /opt/conda/lib/python3.9/site-packages/quasitools/data/hxb2_pol.fas /opt/conda/lib/python3.9/site-packages/quasitools/data/hxb2_pol


echo "############################################"
#	echo "unziping files"
#	gunzip -f $output_dir/*.gz
#	echo "unziping done!!!"
	echo "############################################"
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	mkdir $output_dir/analysis/aafiles/  
	mkdir $output_dir/analysis/vcffiles/  
	mkdir $output_dir/analysis/fastas/
	mkdir $output_dir/analysis/nextclade/
    
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
    #echo $c
	d=$(echo $c | awk -F "/" '{print $7}')
    #echo $d
	echo "starting $d SARS-Cov-2019 Analysis"
	echo "############################################"
	now="$(date +%c)"
	echo -e "starting $d SARS-Cov-2019 Analysis at \t$now" >> "$output_dir/analysis/mensajes.log"
#	echo $c
#	echo $d
	echo "starting $d primers removal"
	echo "################################################"
	bwa mem -t 10 /home/jovyan/work/data/SARS2.fas $read_1 $read_2 | python /home/jovyan/shared/data/Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py /home/jovyan/shared/data/Alt_nCov2019_primers/Primers/ver_niid-200325/primerv3.bed $read_1 $read_2
	echo "starting $d quasitools quality control and mapping"
	echo "################################################"
	quasitools hydra -mf 0.15 -vq 70 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	/home/jovyan/shared/data/makeconsensus-illumina2.sh $output_dir $d     
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	python /home/jovyan/shared/data/changevalueincsv.py $output_dir/analysis/$d $d
	python /home/jovyan/shared/data/changevalueinaavf.py $output_dir/analysis/$d $d
	/home/jovyan/shared/data/generatevariantfilesars2.sh $output_dir $d
	/home/jovyan/shared/data/runingpango.sh $output_dir/analysis/$d/$d.fas $output_dir/analysis/$d/ 4
	cp $output_dir/analysis/$d/*_mutation* $output_dir/analysis/aafiles/.
	cp $output_dir/analysis/$d/$d.fas $output_dir/analysis/fastas/.    
    cp $output_dir/analysis/$d/*_final_* $output_dir/analysis/vcffiles/.
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv
	samtools coverage $output_dir/analysis/$d/align.bam >  $output_dir/analysis/$d/$d.coverage
	samtools stats $output_dir/analysis/$d/align.bam | grep ^SN | cut -f 2- | grep 'average length:' >  $output_dir/analysis/$d/$d.readlength
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/opt/conda/bin/Rscript  /home/jovyan/shared/data/Code_for_reference2.R "SARS-CoV-2019 Genome coverage" $output_dir/analysis/$d/$d.tsv SARS-CoV-2019


	echo "finished $d SARS-Cov-2019 Analysis"
	echo "################################################"
	now="$(date +%c)"
	echo -e "finished $d SARS-Cov-2019 Analysis at \t$now" >> "$output_dir/analysis/mensajes.log"
	done
#	python /home/jovyan/shared/data/mergedv2.py  $output_dir/analysis/coverage/
	echo "drawing final plots of the run"
	echo "################################################"
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/opt/conda/bin/Rscript  /home/jovyan/shared/data/Code_for_reference2.R "SARS-CoV-2019 coverage" $output_dir/analysis/coverage/depth_file.tsv SARS-CoV-2019
	python /home/jovyan/shared/data/plotgroupedcoverage.py $output_dir/analysis

	cat $output_dir/analysis/fastas/*.fas > $output_dir/analysis/fastas/consensus.fasta
	/home/jovyan/shared/data/runingpango.sh $output_dir/analysis/fastas/consensus.fasta $output_dir/analysis/fastas/ 15

/home/jovyan/shared/data/nextclade-Linux-x86_64 --in-order --input-fasta $output_dir/analysis/fastas/consensus.fasta --input-dataset /home/jovyan/shared/data/sars-cov-2/data/sars-cov-2 --output-tsv $output_dir/analysis/nextclade/nextclade.tsv --output-tree $output_dir/analysis/nextclade/nextclade.auspice.json --output-dir $output_dir/analysis/nextclade/ --output-basename nextclade



################################################################
################SARS-CoV-2 v3 not usable Analysis################

elif [ $test = SARS3 ]; then

	echo "starting SARS-Cov-2019 Virus consensus generation"
	cp /home/drakkar4/data/SARS2.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/SARS3.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/SARS2_mutation_db.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	echo "############################################"
	echo "unziping files"
	gunzip -f $output_dir/*.gz
	echo "unziping done!!!"
	echo "############################################"
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	mkdir $output_dir/analysis/aafiles/    
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
	echo "starting $d SARS-Cov-2019 Analysis"
	echo "############################################"
	now="$(date +%c)"
	echo -e "starting $d SARS-Cov-2019 Analysis at \t$now" >> "$output_dir/analysis/mensajes.log"
#	echo $c
#	echo $d
	bwa mem -t 10 /home/drakkar4/data/SARS2.fas $read_1 $read_2 | python /home/drakkar4/data/Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py /home/drakkar4/data/Alt_nCov2019_primers/Primers/ver_niid-200325/primerv3.bed $read_1 $read_2
	quasitools hydra -gc -cp 60 -mf 0.15 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	python /home/drakkar4/data/changevalueincsv.py $output_dir/analysis/$d $d
	python /home/drakkar4/data/changevalueinaavf.py $output_dir/analysis/$d $d
	cp $output_dir/analysis/$d/*_mutation* $output_dir/analysis/aafiles/.
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv
	/home/drakkar4/anaconda3/envs/llama/bin/samtools coverage $output_dir/analysis/$d/align.bam >  $output_dir/analysis/$d/$d.coverage
	/home/drakkar4/anaconda3/envs/llama/bin/samtools stats $output_dir/analysis/$d/align.bam | grep ^SN | cut -f 2- | grep 'average length:' >  $output_dir/analysis/$d/$d.readlength
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/opt/conda/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "SARS-CoV-2019 Genome coverage" $output_dir/analysis/$d/$d.tsv SARS-CoV-2019

	echo "finished $d SARS-Cov-2019 Analysis"
	echo "################################################"
	now="$(date +%c)"
	echo -e "finished $d SARS-Cov-2019 Analysis at \t$now" >> "$output_dir/analysis/mensajes.log"
	done
	#python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	echo "drawing final plots of the run"
	echo "################################################"
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/opt/conda/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "SARS-CoV-2019 coverage" $output_dir/analysis/coverage/depth_file.tsv SARS-CoV-2019
#	cp $output_dir/analysis/*/*_mutation* $output_dir/analysis/aafiles/.

################################################################
################Dengue Virus Analysis#############################

elif [ $test = DENV2 ]; then

	echo "starting Dengue 2 consensus generation"
	cp /mnt/nfs/admin/data/DENV2.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/DENV2.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	echo $c
	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	mv $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	echo "finished $d Dengue 2 Analysis"
	 
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/analysis/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Dengue Virus type 2 polyprotein coverage" $output_dir/analysis/coverage/coveragetotal.csv DENV2
	




elif [ $test = CTL22 ]; then

	echo "starting Clamydia L22 coverage consensus generation"
	cp /home/drakkar4/data/CTampli.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/CTampli.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_CT.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Clamydia L22 coverage" $output_dir/analysis/$d/$d.tsv CTamplicon

	echo "finished $d Clamydia L22 mergedcoverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Clamydia L22 coverage" $output_dir/analysis/coverage/depth_file.tsv CTL22
	





elif [ $test = CTL4 ]; then

	echo "starting Clamydia L4 coverage consensus generation"
	cp /home/drakkar4/data/CTL4.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/CTL4.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_CTL4.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Clamydia L4 coverage" $output_dir/analysis/$d/$d.tsv CTL4

	echo "finished $d Clamydia L4 mergedcoverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Clamydia L4 coverage" $output_dir/analysis/coverage/depth_file.tsv CTL4
	





elif [ $test = T4T ]; then

	echo "starting EColi Phage coverage consensus generation"
	cp /home/drakkar4/data/T4T.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/T4T.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "EColi Phage coverage coverage" $output_dir/analysis/$d/$d.tsv EColiPhage

	echo "finished $d EColi Phage coverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "EColi Phage coverage" $output_dir/analysis/coverage/depth_file.tsv EColiPhage
	





elif [ $test = SHIfago ]; then

	echo "starting Shigela Phage coverage consensus generation"
	cp /home/drakkar4/data/SHIfago.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/SHIfago.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Shigela Phage coverage coverage" $output_dir/analysis/$d/$d.tsv EColiPhage

	echo "finished $d Shigela Phage coverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Shigela Phage coverage" $output_dir/analysis/coverage/depth_file.tsv EColiPhage
	





elif [ $test = ECP1 ]; then

	echo "starting Enterobacter phage coverage consensus generation"
	cp /home/drakkar4/data/EcP1.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/EcP1.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Enterobacter phage coverage" $output_dir/analysis/$d/$d.tsv Enterobacterphage

	echo "finished $d Enterobacter Phage coverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Enterobacter phage coverage" $output_dir/analysis/coverage/depth_file.tsv Enterobacterphage
	





elif [ $test = LEPTO ]; then

	echo "starting Leptospira genome coverage consensus generation"
	cp /home/drakkar4/data/lepto.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/lepto.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "leptospira fragment coverage" $output_dir/analysis/$d/$d.tsv Leptospira

	echo "finished $d leptospira fragment coverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "leptospira fragment coverage" $output_dir/analysis/coverage/depth_file.tsv Leptospira
	




elif [ $test = ompL1 ]; then

	echo "starting Leptosporira sp OmpL1 phage coverage consensus generation"
	cp /home/drakkar4/data/ompL1in.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/ompL1in.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
#	echo $c
#	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "leptospira sp fragment coverage" $output_dir/analysis/$d/$d.tsv Leptospirasp

	echo "finished $d leptospira sp fragment coverage"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "leptospira sp fragment coverage" $output_dir/analysis/coverage/depth_file.tsv Leptospirasp
	


elif [ $test = PB2 ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm PB2 fragment 1 files into reservoir"
elif [ $test = PB1 ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm PB1 fragment 2 files into reservoir"
elif [ $test =  PA ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm PA fragment 3 files into reservoir"
elif [ $test =  HA ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm HA fragment 4 files into reservoir"
elif [ $test =  NP ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm NP fragment 5 files into reservoir"
elif [ $test =  INFLUAH1N1 ]; then
	echo "starting Influenza H1N1 pdm consensus generation"
	cp /mnt/nfs/admin/data/NA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/NA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/NA/
	mkdir $output_dir/NA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NA/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/NA/$d/coverage_file.csv 
	mv $output_dir/NA/$d/coverage_file.csv $output_dir/NA/coverage/coverage_file_NA_$d.csv
	echo "finished $d Neuroaminidase  analysis"
	
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Neuroaminidase gene coverage" $output_dir/NA/coverage/coveragetotal.csv NA
	

	cp /mnt/nfs/admin/data/HA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/HA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/HA/	
	mkdir $output_dir/HA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/HA/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/HA/$d/coverage_file.csv 
	mv $output_dir/HA/$d/coverage_file.csv $output_dir/HA/coverage/coverage_file_HA_$d.csv
	echo "finished $d Hemaglutinine analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/HA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Hemaglutinine gene coverage" $output_dir/HA/coverage/coveragetotal.csv HA

	


	cp /mnt/nfs/admin/data/NP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/NP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NP/	
	mkdir $output_dir/NP/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NP/$d 
	
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/NP/$d/coverage_file.csv 
	mv $output_dir/NP/$d/coverage_file.csv $output_dir/NP/coverage/coverage_file_NP_$d.csv
	echo "finished $d Nucleoprotein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Nucleoprotein gene coverage" $output_dir/NP/coverage/coveragetotal.csv NP
	


	cp /mnt/nfs/admin/data/MP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/MP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/MP/
	mkdir $output_dir/MP/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/MP/$d 

	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/MP/$d/coverage_file.csv 
	mv $output_dir/MP/$d/coverage_file.csv $output_dir/MP/coverage/coverage_file_MP_$d.csv
	
	echo "finished $d Matrix protein analysis"
	done

	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/MP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Matrix protein gene coverage" $output_dir/MP/coverage/coveragetotal.csv MP

	

	cp /mnt/nfs/admin/data/PA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/PA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PA/
	mkdir $output_dir/PA/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PA/$d 

	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/PA/$d/coverage_file.csv 
	mv $output_dir/PA/$d/coverage_file.csv $output_dir/NP/coverage/coverage_file_PA_$d.csv
	
	echo "finished $d Polymerase analysis"
	done		
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Polymerase gene coverage" $output_dir/PA/coverage/coveragetotal.csv PA



	



	cp /mnt/nfs/admin/data/NS.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/NS.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NS/
	mkdir $output_dir/NS/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NS/$d 
	
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/NS/$d/coverage_file.csv 
	mv $output_dir/NS/$d/coverage_file.csv $output_dir/NS/coverage/coverage_file_NS_$d.csv
	echo "finished $d Non Structural protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NS/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Non Structural protein gene coverage" $output_dir/NS/coverage/coveragetotal.csv NS


	


	cp /mnt/nfs/admin/data/PB1.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/PB1.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PB1/	
	mkdir $output_dir/PB1/coverage/
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB1/$d 

	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/PB1/$d/coverage_file.csv 

	mv $output_dir/PB1/$d/coverage_file.csv $output_dir/PB1/coverage/coverage_file_PA_$d.csv
	echo "finished $d Non Structural protein 1 analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB1/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Non Structural protein 1  gene coverage" $output_dir/PB1/coverage/coveragetotal.csv PB1

	
	


	cp /mnt/nfs/admin/data/PB2.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/PB2.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PB2/	
	mkdir $output_dir/PB2/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $5}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB2/$d 

	sed -i -r '1{s/frame: 2/gene_position,'$d'/g}' $output_dir/PB2/$d/coverage_file.csv 
	mv $output_dir/PB2/$d/coverage_file.csv $output_dir/PB2/coverage/coverage_file_PB2_$d.csv
	
	echo "finished $d Non Structural 2 protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB2/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza A H1N1 Non Structural protein 2  gene coverage" $output_dir/PB2/coverage/coveragetotal.csv PB2

	




	/root/TblToFasta.sh $output_dir/NA/*/consensus.fasta > $output_dir/NA/consensusfinal_NA.fa

	/root/TblToFasta.sh $output_dir/HA/*/consensus.fasta > $output_dir/HA/consensusfinal_HA.fa

	/root/TblToFasta.sh $output_dir/NP/*/consensus.fasta > $output_dir/NP/consensusfinal_NP.fa	

	/root/TblToFasta.sh $output_dir/MP/*/consensus.fasta > $output_dir/MP/consensusfinal_MP.fa

	/root/TblToFasta.sh $output_dir/PA/*/consensus.fasta > $output_dir/PA/consensusfinal_PA.fa

	/root/TblToFasta.sh $output_dir/NS/*/consensus.fasta > $output_dir/NS/consensusfinal_NS.fa	

	/root/TblToFasta.sh $output_dir/PB1/*/consensus.fasta > $output_dir/PB1/consensusfinal_PB1.fa	

	/root/TblToFasta.sh $output_dir/PB2/*/consensus.fasta > $output_dir/PB2/consensusfinal_PB2.fa	

	cp $output_dir/NA/*/*_coverage_file_NA* $output_dir/NA/.

elif [ $test =  MP ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm MP fragment 7 files into reservoir"
elif [ $test =  NS1 ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm NS1 fragment 8 files into reservoir"
elif [ $test = INFLUH1N1 ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Influenza H1N1 pdm Neuroaminidase files into reservoir"
elif [ $test = CHIKV ]; then

	echo "starting Chikungunya Virus consensus generation"
	cp /home/drakkar4/data/CHIKV.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/CHIKV.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
	echo $c
	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Chikungunya Virus Genome coverage" $output_dir/analysis/$d/$d.tsv SARS-CoV-2019

	echo "finished $d Chikungunya Virus Analysis"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Chikungunya Virus Genome coverage" $output_dir/analysis/coverage/depth_file.tsv CHIKV
	
    
    

elif [ $test = HIV_subtype ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping HIV Subtype ###  into reservoir"
	
	
	
elif [ $test = INFLUBYAM ]; then
	echo "starting Influenza B Yamagata consensus generation"
	cp /mnt/nfs/admin/data/yamNA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamNA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/NA/
	mkdir $output_dir/NA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NA/$d 
	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/NA/$d/coverage_file.csv 
	mv $output_dir/NA/$d/coverage_file.csv $output_dir/NA/coverage/coverage_file_NA_$d.csv
	echo "finished $d Neuroaminidase  analysis"
	
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Yamagata Neuroaminidase gene coverage" $output_dir/NA/coverage/coveragetotal.csv NA



	cp /mnt/nfs/admin/data/yamHA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamHA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/HA/	
	mkdir $output_dir/HA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/HA/$d 
	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/HA/$d/coverage_file.csv 
	mv $output_dir/HA/$d/coverage_file.csv $output_dir/HA/coverage/coverage_file_HA_$d.csv
	echo "finished $d Hemaglutinine analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/HA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza  B Yamagata Hemaglutinine gene coverage" $output_dir/HA/coverage/coveragetotal.csv HA


	cp /mnt/nfs/admin/data/yamNP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamNP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NP/	
	mkdir $output_dir/NP/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NP/$d 
	
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/NP/$d/coverage_file.csv 
	mv $output_dir/NP/$d/coverage_file.csv $output_dir/NP/coverage/coverage_file_NP_$d.csv
	echo "finished $d Nucleoprotein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Yamagata Nucleoprotein gene coverage" $output_dir/NP/coverage/coveragetotal.csv NP
	

	cp /mnt/nfs/admin/data/yamMP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamMP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/MP/
	mkdir $output_dir/MP/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/MP/$d 

	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/MP/$d/coverage_file.csv 
	mv $output_dir/MP/$d/coverage_file.csv $output_dir/MP/coverage/coverage_file_MP_$d.csv
	
	echo "finished $d Matrix protein analysis"
	done

	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/MP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Yamagata Matrix protein gene coverage" $output_dir/MP/coverage/coveragetotal.csv MP



	

	cp /mnt/nfs/admin/data/yamPA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamPA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PA/
	mkdir $output_dir/PA/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PA/$d 

	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/PA/$d/coverage_file.csv 
	mv $output_dir/PA/$d/coverage_file.csv $output_dir/PA/coverage/coverage_file_PA_$d.csv
	
	echo "finished $d Polymerase analysis"
	done		
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Yamagata Polymerase gene coverage" $output_dir/PA/coverage/coveragetotal.csv PA




	cp /mnt/nfs/admin/data/yamNS.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamNS.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NS/
	mkdir $output_dir/NS/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NS/$d 
	
	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/NS/$d/coverage_file.csv 
	mv $output_dir/NS/$d/coverage_file.csv $output_dir/NS/coverage/coverage_file_NS_$d.csv
	echo "finished $d Non Structural protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NS/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B yamagata H1N1 Non Structural protein gene coverage" $output_dir/NS/coverage/coveragetotal.csv NS






	cp /mnt/nfs/admin/data/yamPB1.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamPB1.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PB1/	
	mkdir $output_dir/PB1/coverage/
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB1/$d 

	sed -i -r '1{s/frame: 2/gene_position,'$d'/g}' $output_dir/PB1/$d/coverage_file.csv 

	mv $output_dir/PB1/$d/coverage_file.csv $output_dir/PB1/coverage/coverage_file_PA_$d.csv
	echo "finished $d Non Structural protein 1 analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB1/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B yamagata Non Structural protein 1  gene coverage" $output_dir/PB1/coverage/coveragetotal.csv PB1



	cp /mnt/nfs/admin/data/yamPB2.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/yamPB2.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA


	
	mkdir $output_dir/PB2/	
	mkdir $output_dir/PB2/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB2/$d 

	sed -i -r '1{s/frame: 2/gene_position,'$d'/g}' $output_dir/PB2/$d/coverage_file.csv 
	mv $output_dir/PB2/$d/coverage_file.csv $output_dir/PB2/coverage/coverage_file_PB2_$d.csv
	
	echo "finished $d Non Structural 2 protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB2/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Yamagata Non Structural protein 2  gene coverage" $output_dir/PB2/coverage/coveragetotal.csv PB2

	


	/root/TblToFasta.sh $output_dir/NA/*/consensus.fasta > $output_dir/NA/consensusfinal_NA.fa

	/root/TblToFasta.sh $output_dir/HA/*/consensus.fasta > $output_dir/HA/consensusfinal_HA.fa

	/root/TblToFasta.sh $output_dir/NP/*/consensus.fasta > $output_dir/NP/consensusfinal_NP.fa	

	/root/TblToFasta.sh $output_dir/MP/*/consensus.fasta > $output_dir/MP/consensusfinal_MP.fa

	/root/TblToFasta.sh $output_dir/PA/*/consensus.fasta > $output_dir/PA/consensusfinal_PA.fa

	/root/TblToFasta.sh $output_dir/NS/*/consensus.fasta > $output_dir/NS/consensusfinal_NS.fa	

	/root/TblToFasta.sh $output_dir/PB1/*/consensus.fasta > $output_dir/PB1/consensusfinal_PB1.fa	

	/root/TblToFasta.sh $output_dir/PB2/*/consensus.fasta > $output_dir/PB2/consensusfinal_PB2.fa	

	cp $output_dir/NA/*/*_coverage_file_NA* $output_dir/NA/.



	
	
elif [ $test = INFLUBVIC ]; then
	echo "starting Influenza B Yamagata consensus generation"
	cp /mnt/nfs/admin/data/vicNA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicNA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/NA/
	mkdir $output_dir/NA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NA/$d 
	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/NA/$d/coverage_file.csv 
	mv $output_dir/NA/$d/coverage_file.csv $output_dir/NA/coverage/coverage_file_NA_$d.csv
	echo "finished $d Neuroaminidase  analysis"
	
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Neuroaminidase gene coverage" $output_dir/NA/coverage/coveragetotal.csv NA



	cp /mnt/nfs/admin/data/vicHA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicHA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/HA/	
	mkdir $output_dir/HA/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/HA/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/HA/$d/coverage_file.csv 
	mv $output_dir/HA/$d/coverage_file.csv $output_dir/HA/coverage/coverage_file_HA_$d.csv
	echo "finished $d Hemaglutinine analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/HA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza  B Victoria Hemaglutinine gene coverage" $output_dir/HA/coverage/coveragetotal.csv HA


	cp /mnt/nfs/admin/data/vicNP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicNP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NP/	
	mkdir $output_dir/NP/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NP/$d 
	
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/NP/$d/coverage_file.csv 
	mv $output_dir/NP/$d/coverage_file.csv $output_dir/NP/coverage/coverage_file_NP_$d.csv
	echo "finished $d Nucleoprotein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Nucleoprotein gene coverage" $output_dir/NP/coverage/coveragetotal.csv NP
	

	cp /mnt/nfs/admin/data/vicMP.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicMP.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/MP/
	mkdir $output_dir/MP/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/MP/$d 

	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/MP/$d/coverage_file.csv 
	mv $output_dir/MP/$d/coverage_file.csv $output_dir/MP/coverage/coverage_file_MP_$d.csv
	
	echo "finished $d Matrix protein analysis"
	done

	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/MP/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Matrix protein gene coverage" $output_dir/MP/coverage/coveragetotal.csv MP



	

	cp /mnt/nfs/admin/data/vicPA.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicPA.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PA/
	mkdir $output_dir/PA/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PA/$d 

	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/PA/$d/coverage_file.csv 
	mv $output_dir/PA/$d/coverage_file.csv $output_dir/PA/coverage/coverage_file_PA_$d.csv
	
	echo "finished $d Polymerase analysis"
	done		
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PA/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Polymerase gene coverage" $output_dir/PA/coverage/coveragetotal.csv PA




	cp /mnt/nfs/admin/data/vicNS.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicNS.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/NS/
	mkdir $output_dir/NS/coverage/	
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/NS/$d 
	
	sed -i -r '1{s/frame: 1/gene_position,'$d'/g}' $output_dir/NS/$d/coverage_file.csv 
	mv $output_dir/NS/$d/coverage_file.csv $output_dir/NS/coverage/coverage_file_NS_$d.csv
	echo "finished $d Non Structural protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/NS/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria H1N1 Non Structural protein gene coverage" $output_dir/NS/coverage/coveragetotal.csv NS






	cp /mnt/nfs/admin/data/vicPB1.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicPB1.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA
	
	mkdir $output_dir/PB1/	
	mkdir $output_dir/PB1/coverage/
	
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB1/$d 

	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/PB1/$d/coverage_file.csv 

	mv $output_dir/PB1/$d/coverage_file.csv $output_dir/PB1/coverage/coverage_file_PA_$d.csv
	echo "finished $d Non Structural protein 1 analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB1/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Non Structural protein 1  gene coverage" $output_dir/PB1/coverage/coveragetotal.csv PB1



	cp /mnt/nfs/admin/data/vicPB2.fas  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas
	cp /mnt/nfs/admin/data/vicPB2.bed  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.bed
	cp /mnt/nfs/admin/data/mutation_db_influ.tsv  /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA.fas /home/drakkar1/miniconda3/envs/influ/lib/python3.7/site-packages/quasitools/data/NA


	
	mkdir $output_dir/PB2/	
	mkdir $output_dir/PB2/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $6}')
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/PB2/$d 

	sed -i -r '1{s/frame: 2/gene_position,'$d'/g}' $output_dir/PB2/$d/coverage_file.csv 
	mv $output_dir/PB2/$d/coverage_file.csv $output_dir/PB2/coverage/coverage_file_PB2_$d.csv
	
	echo "finished $d Non Structural 2 protein analysis"
	done
	python /mnt/nfs/admin/Influ2019/mergedv2.py  $output_dir/PB2/coverage/
	Rscript /mnt/nfs/admin/Influ2019/Code_for_reference.R "Influenza B Victoria Non Structural protein 2  gene coverage" $output_dir/PB2/coverage/coveragetotal.csv PB2

	


	/root/TblToFasta.sh $output_dir/NA/*/consensus.fasta > $output_dir/NA/consensusfinal_NA.fa

	/root/TblToFasta.sh $output_dir/HA/*/consensus.fasta > $output_dir/HA/consensusfinal_HA.fa

	/root/TblToFasta.sh $output_dir/NP/*/consensus.fasta > $output_dir/NP/consensusfinal_NP.fa	

	/root/TblToFasta.sh $output_dir/MP/*/consensus.fasta > $output_dir/MP/consensusfinal_MP.fa

	/root/TblToFasta.sh $output_dir/PA/*/consensus.fasta > $output_dir/PA/consensusfinal_PA.fa

	/root/TblToFasta.sh $output_dir/NS/*/consensus.fasta > $output_dir/NS/consensusfinal_NS.fa	

	/root/TblToFasta.sh $output_dir/PB1/*/consensus.fasta > $output_dir/PB1/consensusfinal_PB1.fa	

	/root/TblToFasta.sh $output_dir/PB2/*/consensus.fasta > $output_dir/PB2/consensusfinal_PB2.fa	

	cp $output_dir/NA/*/*_coverage_file_NA* $output_dir/NA/.




elif [ $test = HBV ]; then
	cp /Users/AlexanderMartinez/testfile.txt /Users/AlexanderMartinez/HINIfile.txt
	echo "coping Hepatitis B Virus files into reservoir"

elif [ $test = SALMO ]; then

	echo "starting Zika Virus consensus generation"
	cp /home/drakkar4/data/senterica.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/SE.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
	echo $c
	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Salmonella genome coverage" $output_dir/analysis/$d/$d.tsv Salmonella
	echo "finished $d Salmonella Analysis"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Salmonella genome coverage" $output_dir/analysis/coverage/depth_file.tsv Salmonella
	


elif [ $test = SALMOPLASMID ]; then

	echo "starting Zika Virus consensus generation"
	cp /home/drakkar4/data/SEplasmid.fas  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas
	cp /home/drakkar4/data/SEplasmid.bed /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.bed
	cp /home/drakkar4/data/mutation_db_influ.tsv  /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/mutation_db.tsv
	bowtie2-build /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol.fas /home/drakkar4/anaconda3/envs/genover2/lib/python3.7/site-packages/quasitools/data/hxb2_pol
	
	gunzip $output_dir/*.gz
	
	mkdir $output_dir/analysis/
	mkdir $output_dir/analysis/coverage/
	for read_1 in $output_dir/*R1*
	do #echo $read_1
	read_2=${read_1/R1/R2}
	c=$(echo $read_1 | awk -F "_S" '{print $1}')
	d=$(echo $c | awk -F "/" '{print $8}')
	echo $c
	echo $d
	quasitools hydra -gc -cp 20 -i $d  $read_1 $read_2 -o $output_dir/analysis/$d 
	sed -i -r '1{s/frame: 0/gene_position,'$d'/g}' $output_dir/analysis/$d/coverage_file.csv 
	cp $output_dir/analysis/$d/coverage_file.csv $output_dir/analysis/coverage/coverage_file_$d.csv
	samtools depth -aa $output_dir/analysis/$d/align.bam | awk -v sample=$d '{$1=sample ; print;}' >  $output_dir/analysis/$d/$d.tsv 
	awk -v sample=$d '{$1=sample ; print;}' $output_dir/analysis/$d/$d.tsv > $output_dir/analysis/coverage/coverage_file_$d.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Salmonella plasmid coverage" $output_dir/analysis/$d/$d.tsv Salmonella_Plasmid
	echo "finished $d Salmonella plasmid Analysis"
	 
	done
	python /home/drakkar4/data/mergedv2.py  $output_dir/analysis/coverage/
	cat $output_dir/analysis/coverage/coverage_file_*.tsv >  $output_dir/analysis/coverage/depth_file.tsv
	/usr/bin/Rscript  /home/drakkar4/data/Code_for_reference2.R "Salmonella plasmid coverage" $output_dir/analysis/coverage/depth_file.tsv Salmonella_Plasmid


else echo "organism $test is not currently supported"

exit 1
fi

if  [ $funct1 = move ]; then
echo "################################################"
	echo "moving fastq files and zipping them"
echo "################################################"
	mkdir $output_dir/analyzedfiles     
	mv $output_dir/*fastq $output_dir/analyzedfiles/.
	pigz $output_dir/analyzedfiles/*fastq
end=$(date +%s.%N)    
runtime=$(python -c "print((${end} - ${start})/60)")

now="$(date +%c)"


echo "################################################"
echo "gencom_beta.v2 Analysis has finished"
echo "################################################"
echo "Finalizado a las $end"
echo "Tiempo total de corrida  $runtime minutos"
echo "################################################"
echo -e "finalizing $test at \t$now \tTotal runtime: \t $runtime minutos" >> "$output_dir/analysis/mensajes.log"
	exit 1
fi




#gunzip $2/*.gz

#echo $read_1 $read_2 $output_dir

#for read_1 in $output_dir/*R1*
#do #echo $read_1
#read_2=${read_1/R1/R2}
#c=$(echo $read_1 | awk -F "_S" '{print $1}')
#quasitools hydra -gc -cp 20 -i $c.consensus  $read_1 $read_2 -o $c 
#echo "finished $read_1 analysis"
#done

#conda deactivate

