#! /bin/bash
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=20gb
#PBS -j oe
#PBS -m bea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N create_subbasin

# get number of cpus
#export NCPUS=`cat ${PBS_NODEFILE} | wc -l`
NCPUS=10

# OMP Setting
export OMP_NUM_THREADS=$NCPUS            #   OpenMP cpu num

### goto path ##
cd $PBS_O_WORKDIR

### SET CaMa-Flood directory
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v395b_20191030/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396_20191225/"
# CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CAMADIR="/cluster/data6/menaka/CaMa-Flood_v4"
#CAMADIR="/work/a01/modi/cama_v395"

###SET map name
# map="glb_06min"
# map="glb_15min"
#map="reg_06min_srtm"
map="amz_06min"

####output
outdir="./output"

### threshold
threshold=100000 #1.0e11 #

### precent
perct=0.30

echo ./scripts/subbasin $map $CAMADIR $outdir $threshold $perct
time ./scripts/subbasin $map $CAMADIR $outdir $threshold $perct

#echo ./get_subbasin $map $CAMADIR $outdir $threshold
#time ./get_subbasin $map $CAMADIR $outdir $threshold