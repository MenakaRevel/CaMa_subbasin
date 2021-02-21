# /bin/sh
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -j oe
#PBS -m bea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N create_basinID

# get number of cpus
export NCPUS=`cat ${PBS_NODEFILE} | wc -l`

# OMP Setting
export OMP_NUM_THREADS=$NCPUS            #   OpenMP cpu num

### goto path ##
#cd $path
cd $PBS_O_WORKDIR

### SET CaMa-Flood directory
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa_package_v394_20190508/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v395b_20191030/"
#CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396_20191225/"
CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
#CAMADIR="/work/a01/modi/cama_v395"

###SET map name
#map="glb_06min"
map="amz_06min"
#map="reg_06min_srtm"

####output
outdir="../output"

# make out dir
mkdir $outdir

###
echo ./basinID $CAMADIR $map $outdir
time ./basinID $CAMADIR $map $outdir