set -xe
dataset="Staph"
asms=("Bambus2")
# crashed  and worked Bambus2, MSR-CA, SOAPdenovo, Velvet, ABySS, SGA, Allpaths-LG
#asms=("ABySS" "ABySS2" "Allpaths-LG" "Bambus2" "MSR-CA" "SGA" "SOAPdenovo" "Velvet")
count=2
#====================================================
str="_"
folder_name=""
folder_name="$folder_name$dataset"
f=$HOME"/GAGE_Data/"
f1=$HOME"/Results/"
sl="/"
folder_name="$f1$folder_name$sl$asm"
rm -rf Gaps
mkdir -p Gaps
mkdir -p "$folder_name"

r="/Data/original/"

read1="frag_1.fastq"
read2="frag_2.fastq"

#read1="shortjump_1.fastq"
#read2="shortjump_2.fastq"

read_dir1=$f$dataset$r$read1
read_dir2=$f$dataset$r$read2

#i="z_index/indexFile"
#index_dir=$f$i

ref="/Data/original/genome.fasta"
ref_genome=$f$dataset$ref

asm_dir="/Assembly/"
scaffold="/genome.scf.fasta"
for asm in "${asms[@]}"
do
  rm -rf Gaps
  mkdir Gaps
  gapped_genome=$dataset$asm_dir$asm$scaffold
  rm -rf indexFile
  rm -rf result.sam
  $bowtie2_build $gapped_genome indexFile
  $bowtie2 -x indexFile -1   $read_dir1  -2  $read_dir2  -S result.sam

  echo "running main.cpp"
  g++ main.cpp && ./a.out $gapped_genome

  echo "running figbird"
  g++ figbird.cpp && ./a.out $gapped_genome 1

  echo "measuing distance"
  mv gapout.txt gapOut2.txt

  arg1="gapOut2.txt"

  g++ distance.cpp
  res2="$(./a.out $arg1)"
  res3=$res2

  #Iterative method====================

  echo Count = $count ED = $res1
  new_file="filledContigs.fa"

  while [ $count -lt 2 ]
  do
    if [ $res2 -eq 0  ];then
        break
    else
        $bowtie2_build $new_file $index_dir
        $bowtie2 -x $index_dir -1 $read_dir1 -2 $read_dir2 -S result.sam

        g++ main.cpp && ./a.out $new_file
        #echo "While loop main ends..."

        g++ figbird.cpp && ./a.out $new_file 1
        #echo "While loop figbird ends..."

        #g++ new_gap.cpp && ./a.out
        #mv new.txt gapOut2.txt

        #arg1="gapOut2.txt"
        #g++ distance.cpp
        res2="$(./a.out $arg1)"
    fi
    count=$((count+1))
    #echo Count = $count ED = $res1
  done
  # format the output

  #pyenv shell 2.7.6
  mkdir -p "EM/"$dataset
  python reference.py $new_file "EM/"$dataset"/filledContigs_"$asm".fa"
  #Evaluate using quast

  #python quast.py --strict-NA -R {input.ref} -o output --no-plots {input.scaffolds}
  #python $QUAST --strict-NA -R $ref_genome -o output --no-plots $dataset$asm_dir$asm"/filledContigs.fa"
  # QUAST quast_correction
  #python $HOME"/eval_Gap2Seq/scripts/correct_quast.py" --N 4000 {input.quast_stdout {input.quast_misassm_file} {input.quast_report} {input.scaffolds}",iterable=True))[0]
  #result=$(python $HOME"/eval_Gap2Seq/scripts/correct_quast.py" --N 4000 output/contigs_reports/contigs_report_*.stdout output/report.txt output/contigs_reports/misassemblies_report.txt  output/report.txt)
  #printf $result >> out.txt
done
