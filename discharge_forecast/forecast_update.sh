date

conda_path=/nird/home/owul/miniforge3

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/data/regular_downloads/metno/smaakraft/metno_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv ]
then
source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast
# run forecast update routine:
ipython /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/update_forecast.py

else
echo not running forecast update, forecast already exists.
date
fi