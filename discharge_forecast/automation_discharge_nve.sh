date

conda_path=/nird/home/owul/miniforge3

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/data/regular_downloads/sildre_nve/nve/nve_$(date +%Y)-$(date +%m)-$(date +%d).csv ]
then
source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast
# run plotting routine:
ipython /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/update_sildre.py
else
echo NVE discharge already updated. $date
fi