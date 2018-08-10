#$ -V
#$ -cwd
#$ -P fedorov-iles.prj
#$ -q serial-low.q
#$ -j y
#$ -o out.$JOB_ID
#$ -ac runtime="1h"

date
./run3drismgaff.sh benzyl_alcohol
date
