date

conda_path=/nird/home/owul/miniforge3

source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/smaakraft/daily_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv ]
then
# make HBV model run for smaakraft:
Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_forecast_smaakraft.R
# push to web-exposed directory:
rsync -avz /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/smaakraft/ /projects/NS9873K/www/smaakraft/discharge_forecast/
else
echo HBV forecast runs for smaakraft already done. $date
fi

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/nve/daily_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv ]
then
# make HBV model run for nve catchments:
Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_forecast_nve.R
else
echo HBV forecast runs for nve already done. $date
fi