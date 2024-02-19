date

conda_path=/nird/home/owul/miniforge3

if [ ! -e /projects/NS9873K/owul/projects/discharge_forecast/data/regular_downloads/senorge/smaakraft/seNorge_$(date +%Y)$(date +%m)$(date +%d).csv ]
then
source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast
# update latest seNorge data:
ipython /projects/NS9873K/owul/projects/discharge_forecast/discharge_forecast/update_senorge.py
# create merged files that include latest senorge data and forecast intialized on day of running this script:
ipython /projects/NS9873K/owul/projects/discharge_forecast/discharge_forecast/update_merge_file.py
else
echo seNorge already updated. $date
fi