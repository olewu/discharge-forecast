Workflow:
* calibrate_everything.R : Calibrate model for NVE catchments.
* make_runoff_simulations_smaakraft.R : Find donors for smaakraft catchments. Make simulations.


Evaluate for NVE catchments:
* donor_search: Find donors for the NVE catchments.
* test_our_donors: Evaluate performance for NVE catchments (when donors are used).
* make_runoff_simulations_nve: More evaluations


Data preparation:
* nve_api: Download streamflow data from NVE api.
* nve_api_feltparam: Download catchment properties from NVE api.
* get_seNorge_nve: Download seNorge temp/precip for NVE catchments.
* get_seNorge_smaakraft: Download seNorge temp/precip for Smaakraft catchments.
* scripts/DataDownload/catchment_data_from_pdf_thea.py: Read data from pdfs with Ole's python script
* get_feltparam_smaakraft.R: A lot of data cleaning ++ for preparing smaakraft "feltparam" on the same format as NVE catchment properties.


Smaakraft/Results/historical_simulations/date: Here the final streamflow predictions are stored.