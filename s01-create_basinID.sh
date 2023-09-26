# /bin/sh
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -j oe
#PBS -m bea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N create_basinID

# get number of cpus
# export NCPUS=`cat ${PBS_NODEFILE} | wc -l`
NCPUS=20

# OMP Setting
export OMP_NUM_THREADS=$NCPUS            #   OpenMP cpu num

### goto path ##
#cd $path
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/CaMa_subbasin"

### SET CaMa-Flood directory
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v395b_20191030/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396_20191225/"
# CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# CAMADIR="/cluster/data6/menaka/CaMa-Flood_v4"
#CAMADIR="/work/a01/modi/cama_v395"
# CAMADIR="/work/a07/uddin/Cama-Flood/CaMa-Flood_v4"
CAMADIR="/cluster/data6/menaka/CaMa-Flood_v410"

###SET map name
# map="glb_06min"
# map="glb_15min"
#map="reg_06min_srtm"
# map="amz_06min"
# map="conus_06min"
# map="GBM_15min"
map="yangtze_05min"

####output
outdir="./output"

# make out dir
mkdir -p $outdir

###
echo ./scripts/basinID $CAMADIR $map $outdir
time ./scripts/basinID $CAMADIR $map $outdir