date

conda_path=/nird/home/owul/miniforge3

if [ ! -e /projects/NS9873K/owul/projects/discharge_forecast/data/regular_downloads/metno_21d/smaakraft/metno_$(date +%Y)-$(date +%m)-$(date +%d).csv ]
then
source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast
# run forecast update routine:
ipython /projects/NS9873K/owul/projects/discharge_forecast/discharge_forecast/update_forecast_21d.py

else
echo not running forecast update, forecast already exists.
fi

# if the latest 10d forecast exists (which serves the missing 18h between senorge and the 21d forecast),
# perform the copula coupling
fc10d_smaakraft=/projects/NS9873K/owul/projects/discharge_forecast/data/regular_downloads/metno/smaakraft/metno_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv
fc21d_ppens_smaakraft=/projects/NS9873K/owul/projects/discharge_forecast/results/ens_forecast_input/smaakraft/fc_init_$(date +%Y)-$(date +%m)-$(date +%d).csv
fc10d_nve=/projects/NS9873K/owul/projects/discharge_forecast/data/regular_downloads/metno/nve/metno_$(date +%Y)-$(date +%m)-$(date +%d)T06:00:00Z.csv
fc21d_ppens_nve=/projects/NS9873K/owul/projects/discharge_forecast/results/ens_forecast_input/nve/fc_init_$(date +%Y)-$(date +%m)-$(date +%d).csv
if [ -e "$fc10d_smaakraft" ] && [ ! -e "$fc21d_ppens_smaakraft" ] && [ -e "$fc10d_nve" ] && [ ! -e "$fc21d_ppens_nve" ]
then
source $conda_path/etc/profile.d/conda.sh # look for <base_path>/etc/profile.d/conda.sh in <base_path> given by conda info | grep -i 'base environment'
# activate the environment:
conda activate discharge_forecast
# run forecast update routine:
ipython /projects/NS9873K/owul/projects/discharge_forecast/discharge_forecast/update_merge_21d_forecast.py

else
echo "not post-processing 21d forecast, 10d forecast doesn't exist yet or the 21d file is already there."
date
fi