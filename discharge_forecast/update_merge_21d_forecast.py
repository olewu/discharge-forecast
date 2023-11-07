import pandas as pd
from datetime import date, timedelta
import numpy as np
from scipy.optimize import minimize
from scipy.stats import gamma, skewnorm

from discharge_forecast.config import *

# functions for fitting gamma and skewnorm distributions:
def gamma_objective(params,percs,obs_val):
    shape, scale = params
    expected_values = gamma.ppf(percs, a=shape, scale=scale)
    return np.sum((obs_val - expected_values) ** 2)

def skewnorm_objective(params,percs,obs_val):
    a, loc, scale = params
    expected_values = skewnorm.ppf(percs, a, loc = loc, scale=scale)
    return np.sum((obs_val - expected_values) ** 2)

def derive_gamma_percs(obs, percentiles=np.array([.1,.5,.9]), zero_thresh=.1,  des_percs = np.arange(0.05,1,.05),initial_guess = [1,2]):
    new_ens_quant_eqspace = []
    for observed_values in obs.T:
        # minimize the SSE for the given percentiles 
        result = minimize(gamma_objective, initial_guess, args = (percentiles, observed_values), method='Nelder-Mead')
        # extract distribution parameters:
        shape, scale = result.x
        # draw equally spaced percentiles fromt the estimated distribution:
        new_ens_quant_eqspace.append(gamma.ppf(des_percs, a=shape, scale=scale))
    
    new_ens_quant_eqspace = np.array(new_ens_quant_eqspace).T
    
    # censor by setting values close to 0 exactly zero (arbitrary threshold of 0.1 mm):
    new_ens_quant_eqspace[new_ens_quant_eqspace < zero_thresh] = 0

    return new_ens_quant_eqspace

def derive_skewnorm_percs(obs, percentiles=np.array([.1,.5,.9]), des_percs = np.arange(0.05,1,.05),initial_guess = [0,0,1]):
    new_ens_quant_eqspace_st = []
    for observed_values in obs.T:
        # minimize the SSE for the given percentiles 
        result = minimize(skewnorm_objective, initial_guess, args = (percentiles, observed_values), bounds = [(-50,50),(-50,50),(0,20)], method='L-BFGS-B')
        # extract distribution parameters:
        a, loc, scale = result.x
        # draw equally spaced percentiles fromt the estimated distribution:
        new_ens_quant_eqspace_st.append(skewnorm.ppf(des_percs, a, loc = loc, scale=scale))

    return np.array(new_ens_quant_eqspace_st).T


def pp_21d_ens(catchments_from, rdate = date.today(), des_percs = np.arange(0.05,1,0.05), back_extend = False):

    tmrrw = rdate  + timedelta(days=1)

    # name of file to write:
    outf = proj_base + '/results/ens_forecast_input/{1:s}/fc_init_{0:}.csv'.format(rdate,catchments_from)

    # load corresponding 21-day forecast:
    fc21 = pd.read_csv(proj_base + '/data/regular_downloads/metno_21d/{1:s}/metno_{0:}.csv'.format(rdate,catchments_from))

    # load 10-day forecast for the first forecast day:
    fc10 = pd.read_csv(proj_base + '/data/regular_downloads/metno/{1:s}/metno_{0:}T06:00:00Z.csv'.format(rdate,catchments_from))
    fc10_init = fc10[fc10.date == tmrrw.strftime('%Y-%m-%d')]

    #--------create a probabilistic forecast out of the three percentiles of the 21d forecast by fitting parametric distributions--------#
    #-------- then draw a number (20) of equally spaced percentiles from the fitted distribution--------#
    # reorganize data:
    prec_fc = np.concatenate([[fc21.prec_10.to_list()],[fc21.prec_50.to_list()],[fc21.prec_90.to_list()]]).round(2)
    st_fc = np.concatenate([[fc21.st_10.to_list()],[fc21.st_50.to_list()],[fc21.st_90.to_list()]]).round(2)

    prec_resampled = derive_gamma_percs(prec_fc, des_percs=des_percs)
    st_resampled = derive_skewnorm_percs(st_fc, des_percs=des_percs)

    prec_df = pd.concat([fc21.date,fc21.catchname,pd.DataFrame(prec_resampled.T)],axis=1)
    st_df = pd.concat([fc21.date,fc21.catchname,pd.DataFrame(st_resampled.T)],axis=1)

    #--------re-order using inter-dependence structure from observations ('Schaake-shuffle')--------#
    # split by station
    # load historical senorge data back to 1960 (two seperate files):
    sn_hist = pd.read_csv(proj_base + '/data/historical_data/senorge/{0:s}/seNorge_daily_{0:s}.csv'.format(catchments_from)).rename(columns={'catchid':'catchname'}).drop_duplicates()
    if back_extend: # not possible at the moment for nve data
        sn_hist60 = pd.read_csv(proj_base + '/data/historical_data/senorge/{0:s}/seNorge_daily_{0:s}_1960-1989.csv'.format(catchments_from))
        # concatenate:
        sn_hist = pd.concat([sn_hist60,sn_hist])

    # subset historical data to only include same calendar period that the forecast covers:
    calday_list = [dt[5:] for dt in fc21.date.unique()]

    datesel = [ii for ii,dt in enumerate(sn_hist.date) if dt[5:] in calday_list]

    sn_cal_day_sel = sn_hist.iloc[datesel].copy()
    sn_cal_day_sel['year'] = pd.to_datetime(sn_cal_day_sel.date).dt.year

    y_unique = sn_cal_day_sel.year.unique()
    randy = np.random.choice(y_unique, size=len(prec_resampled), replace=False)

    for ii,station in enumerate(sn_cal_day_sel.catchname.unique()):

        if station not in st_df.catchname.unique():
            continue

        sn_sel = sn_cal_day_sel[sn_cal_day_sel.catchname == station]

        TEMPL = np.array([np.concatenate([sn_sel[sn_sel.year == YY].st.to_numpy(), sn_sel[sn_sel.year == YY].prec.to_numpy()]) for YY in randy])
        RANKS = TEMPL.argsort(0).argsort(0)

        stat_sel_st = st_df[st_df.catchname == station]
        stat_sel_st_np = stat_sel_st.drop(columns=['date','catchname']).to_numpy()

        stat_sel_prec = prec_df[prec_df.catchname == station]
        stat_sel_prec_np = stat_sel_prec.drop(columns=['date','catchname']).to_numpy()

        fc_col = np.concatenate([stat_sel_st_np,stat_sel_prec_np],axis=0).T

        ens_reshuffled = fc_col[RANKS,np.arange(fc_col.shape[1])]

        colnames_st = ['st_{:d}'.format(ensn) for ensn in range(len(des_percs))]
        colnames_prec = ['prec_{:d}'.format(ensn) for ensn in range(len(des_percs))]
        
        st_init = fc10_init[fc10_init.catchname == station].st.values.repeat(len(des_percs))[np.newaxis,:]
        prec_init = fc10_init[fc10_init.catchname == station].prec.values.repeat(len(des_percs))[np.newaxis,:]

        st_DF = pd.DataFrame(np.concatenate([st_init,ens_reshuffled[:,:len(stat_sel_st.date)].T],axis=0),columns=colnames_st)
        prec_DF = pd.DataFrame(np.concatenate([prec_init,ens_reshuffled[:,len(stat_sel_st.date):].T],axis=0),columns=colnames_prec)
        
        head = pd.concat([fc21[fc21.catchname==station].date,fc21[fc21.catchname==station].catchname],axis=1)
        head = pd.concat([fc10_init[fc10_init.catchname == station][['date','catchname']],head],axis=0).reset_index(drop=True)
        # append to output file:
        pd.concat([head,st_DF,prec_DF],axis=1).to_csv(outf,mode='a',index=False,header=not os.path.exists(outf))

        # if ii == 0:
        #     pp_fc = pd.concat([head,st_DF,prec_DF],axis=1)
        # else:
        #     pp_fc = pd.concat([pp_fc,pd.concat([head,st_DF,prec_DF],axis=1)],axis=0)


if __name__ == '__main__':
    for catchmnts in ['smaakraft','nve']:
        pp_21d_ens(catchmnts, des_percs=np.linspace(0,1,33)[1:-1])