date

source /nird/home/owul/.bashrc

conda_path=/nird/home/owul/miniforge3

source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/smaakraft/daily_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv ]
then
# make HBV model run for smaakraft:
Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_forecast_smaakraft.R
# push to web-exposed directory:
rsync -avz /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/smaakraft/*.csv /projects/NS9873K/www/smaakraft/discharge_forecast/
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

ensemble_members=$(seq 1 31)

if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/smaakraft/daily21d_$(date +%Y)-$(date +%m)-$(date +%d) ]
then
# make HBV model ensemble run for smaakraft catchments:
# Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_ensemble_forecast_smaakraft.R
# run the ensemble members in parallel
echo "$ensemble_members" | xargs -I {} -P $(echo "$ensemble_members" | wc -w) bash -c "Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_singleensemble_forecast_smaakraft.R {}"
else
echo HBV ensemble forecast runs for smaakraft already done. $date
fi


if [ ! -e /projects/NS9001K/owul/projects/discharge_forecast/results/discharge_forecast/nve/daily21d_$(date +%Y)-$(date +%m)-$(date +%d) ]
then
# make HBV model ensemble run for nve catchments:
# Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_ensemble_forecast_nve.R
# run the ensemble members in parallel:
echo "$ensemble_members" | xargs -I {} -P $(echo "$ensemble_members" | wc -w) bash -c "Rscript /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/RHBV/scripts/make_runoff_singleensemble_forecast_nve.R {}"
ipython /projects/NS9001K/owul/projects/discharge_forecast/discharge_forecast/compute_discharge_percentiles.py
else
echo HBV ensemble forecast runs for nve already done. $date
fi