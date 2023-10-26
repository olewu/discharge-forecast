#!/usr/bin/env python

import properscoring as ps
from sklearn.metrics import mean_pinball_loss
from scipy.stats import spearmanr, pearsonr
import xarray as xr
import numpy as np
import os
import pandas as pd
import sys

from discharge_forecast.config import proj_base

obs_dir = os.path.join(proj_base,'/data/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv')

class verification():
    def __init__(self,find_dir,lead_time,prob_dim='number',ndec=3,obs_dir=obs_dir):

        self.prob_dim = prob_dim
        self.ndec = ndec
        self.lead_time = lead_time

        #-------------- Observations --------------#
        self.obs_dir = obs_dir

        # load senorge data in catchment
        seno_catch_df = pd.read_csv(obs_dir)

        seno_catch_df.date = pd.to_datetime(seno_catch_df.date)   

        # convert to xarray.dataset
        seno_catch = seno_catch_df.set_index(['date','catchname']).to_xarray()
        seno_catch = seno_catch.rename({'date':'valid_time','catchname':'station'})
        seno_catch = seno_catch.assign_coords(valid_time = seno_catch.valid_time + pd.Timedelta(hours=6))

        
        #-------------- Forecasts --------------#
        self.find_dir = find_dir

        self.t2m_file = '{0:s}t2m/t2m_{1:0>2d}d.nc'.format(find_dir,lead_time)
        with xr.open_dataarray(self.t2m_file).transpose('station','valid_time',self.prob_dim) as da1:
            da_t2m  = da1

        self.tp_file = self.t2m_file.replace('t2m','tp')
        with xr.open_dataarray(self.tp_file).transpose('station','valid_time',self.prob_dim)  as da2:
            da_tp  = da2

        
        #-------------- Merge into xr.Dataset --------------#
        # subselect matching dates (and stations) from senorge:
        seno_catch_sel = seno_catch.sel(valid_time=da_t2m.valid_time,station=da_t2m.station).transpose('station','valid_time')

        # merge into one xr.Dataset
        self.verification_dataset = xr.Dataset(
            dict(
                obs_tp = seno_catch_sel.prec,
                obs_t2m = seno_catch_sel.st,
                fc_tp = da_tp,
                fc_t2m = da_t2m,
            )
        )
        
        self.nstations = len(self.verification_dataset.station)
        self.nobs = len(self.verification_dataset.valid_time)
        self.nens = len(self.verification_dataset[self.prob_dim])

        # initialize an xr.Dataset to collect verification metrics in
        #TODO: convert valid_date axis to initialization date instead!
        self.vf_collection = xr.Dataset(coords=self.verification_dataset.coords).drop_vars(self.prob_dim)

        return

    def crps(self):

        if self.mode == 'ensemble':
            self.vf_collection['crps_t2m'] = (
                ('station','valid_time'),
                ps.crps_ensemble(
                    self.verification_dataset.obs_t2m.values.flatten(),
                    self.verification_dataset.fc_t2m.values.reshape(-1,len(self.verification_dataset.number))
                ).reshape(self.nstations,self.nobs)
            )
            self.vf_collection['crps_tp'] = (
                ('station','valid_time'),
                ps.crps_ensemble(
                    self.verification_dataset.obs_tp.values.flatten(),
                    self.verification_dataset.fc_tp.values.reshape(-1,len(self.verification_dataset.number))
                ).reshape(self.nstations,self.nobs)
            )

        return
    
    def crps_fair(self):
        crps_fair_t2m = CRPS_ensemble(
            self.verification_dataset.obs_t2m.values.flatten(),
            self.verification_dataset.fc_t2m.values.reshape(-1,len(self.verification_dataset.number)).transpose(1,0)
            )
        
        self.vf_collection['crps_fair_t2m'] = (
            ('station','valid_time'),
            crps_fair_t2m.reshape(self.nstations,self.nobs)
        )

        crps_fair_tp = CRPS_ensemble(
            self.verification_dataset.obs_tp.values.flatten(),
            self.verification_dataset.fc_tp.values.reshape(-1,len(self.verification_dataset.number)).transpose(1,0)
            )
        
        self.vf_collection['crps_fair_tp'] = (
            ('station','valid_time'),
            crps_fair_tp.reshape(self.nstations,self.nobs)
        )

        return
    
    def crps_decomposition(self,by=None):

        components = ['REL','RES','UNC']
        ix = [1,2,3]

        # reshape:
        if by == None:
            obs_t2m_ = self.verification_dataset.obs_t2m.values.flatten()[:,np.newaxis]
            fc_t2m_ = self.verification_dataset.fc_t2m.values.reshape(-1,len(self.verification_dataset.number)).transpose(1,0)[:,:,np.newaxis]
            obs_tp_ = self.verification_dataset.obs_tp.values.flatten()[:,np.newaxis]
            fc_tp_ = self.verification_dataset.fc_tp.values.reshape(-1,len(self.verification_dataset.number)).transpose(1,0)[:,:,np.newaxis]
            dims = ()
        elif by == 'station':
            obs_t2m_ = self.verification_dataset.obs_t2m.values.transpose(1,0)
            fc_t2m_ = self.verification_dataset.fc_t2m.values.transpose(2,1,0)
            obs_tp_ = self.verification_dataset.obs_tp.values.transpose(1,0)
            fc_tp_ = self.verification_dataset.fc_tp.values.transpose(2,1,0)
            dims = ('station')
        # elif by == 'season': # last dimension 4, -2 dimension nstations * nvalid_times/4
        #     obs_t2m_ = self.verification_dataset.obs_t2m.values.transpose(1,0)
        #     fc_t2m_ = self.verification_dataset.fc_t2m.values.transpose(2,1,0)
        #     obs_tp_ = self.verification_dataset.obs_tp.values.transpose(1,0)
        #     fc_tp_ = self.verification_dataset.fc_tp.values.transpose(2,1,0)
        # elif by == 'stationseason': # last dimension 4 * nstations, -2 dimension nvalid_times/4
        # #     obs_t2m_ = self.verification_dataset.obs_t2m.values.transpose(1,0)
        # #     fc_t2m_ = self.verification_dataset.fc_t2m.values.transpose(2,1,0)
        # #     obs_tp_ = self.verification_dataset.obs_tp.values.transpose(1,0)
        # #     fc_tp_ = self.verification_dataset.fc_tp.values.transpose(2,1,0)
        

        crps_dcmp_t2m = CRPS_dcmpstn(obs_t2m_,fc_t2m_)
        
        for cmpnt,index in zip(components,ix):
            self.vf_collection['crps_{:}_t2m'.format(cmpnt)] = (
                dims,
                crps_dcmp_t2m[index]
            )

        crps_dcmp_tp = CRPS_dcmpstn(obs_tp_,fc_tp_)
        
        for cmpnt,index in zip(components,ix):
            self.vf_collection['crps_{:}_tp'.format(cmpnt)] = (
                dims,
                crps_dcmp_tp[index]
            )


        return
    
    def pinball_loss(self,alpha = [.1,.25,.33,.5,.66,.75,.9]):

        self.alpha = alpha

        if 'alpha' not in self.vf_collection.dims:
            self.vf_collection = self.vf_collection.expand_dims(alpha=xr.DataArray(self.alpha,[('alpha',self.alpha)]))

        mploss_t2m = np.array([mean_pinball_loss(self.verification_dataset.obs_t2m.values.flatten()[np.newaxis,:],np.quantile(self.verification_dataset.fc_t2m.values.reshape(-1, self.nens)[np.newaxis,:],alph,axis=-1),alpha=alph,multioutput='raw_values').round(self.ndec) for alph in self.alpha])
        self.vf_collection['pinball_loss_t2m'] = (
            ('alpha','station','valid_time'),
            mploss_t2m.reshape(len(self.alpha),self.nstations,self.nobs)
        )

        mploss_tp= np.array([mean_pinball_loss(self.verification_dataset.obs_tp.values.flatten()[np.newaxis,:],np.quantile(self.verification_dataset.fc_tp.values.reshape(-1, self.nens)[np.newaxis,:],alph,axis=-1),alpha=alph,multioutput='raw_values').round(self.ndec) for alph in self.alpha])
        self.vf_collection['pinball_loss_tp'] = (
            ('alpha','station','valid_time'),
            mploss_tp.reshape(len(self.alpha),self.nstations,self.nobs)
        )

        return
    
    def brier_score(self,thresholds=[0,0.5,1,5,10,50]):

        for thresh in thresholds:
            
            self.vf_collection['bs_precipyn_{:}'.format(thresh)] = (
                ('station','valid_time'),
                ps.brier_score(
                    self.verification_dataset.obs_tp > thresh,
                    (self.verification_dataset.fc_tp > thresh).sum('number')/self.nens
                )
            )

        return



    def corr(self,corrdim='valid_time',statistic='mean'):

        if statistic == 'mean':
            self.vf_collection['corr_t2m'] = xr.corr(self.verification_dataset.obs_t2m,self.verification_dataset.fc_t2m.mean('number'),dim=[corrdim])
            self.vf_collection['corr_tp'] = xr.corr(self.verification_dataset.obs_tp,self.verification_dataset.fc_tp.mean('number'),dim=[corrdim])

        return
    
    def rank_corr(self,corrdim='valid_time',statistic='mean'):

        if statistic == 'mean':
            spr_t2m = np.array([spearmanr(self.verification_dataset.obs_t2m.values[i],self.verification_dataset.fc_t2m.mean('number').values[i]).statistic for i in range(self.nstations)])
            spr_tp = np.array([spearmanr(self.verification_dataset.obs_tp.values[i],self.verification_dataset.fc_tp.mean('number').values[i]).statistic for i in range(self.nstations)])

        self.vf_collection['rankcorr_t2m'] = (
            ('station'),
            spr_t2m
        )
        self.vf_collection['rankcorr_tp'] = (
            ('station'),
            spr_tp
        )

        return


    def save_coll(self,extension=''):
        savename = 'verification_{0:0>2d}{1:s}d.nc'.format(self.lead_time,extension)
        savedir = os.path.join(self.find_dir,'verification')
        os.makedirs(savedir,exist_ok=True)
        self.savepath = os.path.join(savedir,savename)
        self.vf_collection.to_netcdf(self.savepath)

    def full_eval(self,extension=''):
        self.crps()
        self.pinball_loss()
        self.crps_fair()
        try:
            self.crps_decomposition(by='station')
        except:
            print('CRPS decomposition failed for lead time {:}d'.format(self.lead_time))
            pass
        self.corr()
        self.rank_corr()
        self.brier_score()

        self.save_coll(extension)

    def load_eval(self):



        return



def CRPS_ensemble(obs,fc,fair=True,axis=0):
    """
    @author: Ole Wulff
    @date: 2020-07-08
    
    implementation of fair (adjusted) CRPS based on equation (6) from Leutbecher (2018, QJRMS, https://doi.org/10.1002/qj.3387)
    version with fair=False tested against properscoring implementation crps_ensemble (see https://pypi.org/project/properscoring/)
    
    INPUT:
        obs: observations as n-dimensional array
        fc: forecast ensemble as (n+1)-dimensional where the extra dimension (axis) carries the ensemble members
        fair: if True returns the fair version of the CRPS accounting for the limited ensemble size (see Leutbecher, 2018)
              if False returns the normal CRPS
        axis: axis of fc array that contains the ensemble members, defaults to 0
    OUTPUT:
        CRPS: n-dimensional array
    TODO:
        implement weights for ensemble member weighting
    """
    odims = obs.shape
    M = fc.shape[axis]
    if axis != 0:
        fc = np.swapaxes(fc,axis,0)
    
    # flatten all dimensions except for the ensemble member dimension:
    fc_flat = fc.reshape([M,-1])
    obs_flat = obs.reshape([-1])
    
    dsum = np.array([abs(fc_flat[jj] - fc_flat[kk]) for kk in range(M) for jj in range(M)]).sum(axis=axis)
    if fair:
        CRPS = 1/M * (abs(fc_flat - obs_flat)).sum(axis=axis) - 1/(2*M*(M-1)) * dsum
    else:
        CRPS = 1/M * (abs(fc_flat - obs_flat)).sum(axis=axis) - 1/(2*M**2) * dsum
        
    # is this necessary or even a good idea at all?
#     del dsum, fc_flat, obs_flat
    
    return CRPS.reshape([*odims])


def CRPS_dcmpstn(obs,fc,ensax=0,smpax=1):
    """
    estimation and decomposition of the CRPS after Hersbach (2000)
    
    fc must have 3 dimensions:
        one ensemble member axis, one sample axis (initializations) and one iteration axis only the last dim will remain
    obs must have one less dimension
    
    only works for ensax = 0 and smpax = 1 at the moment...
    
    updated 10/22/2020
        
    """
    ninit = fc.shape[smpax]
    fc = fc.reshape([fc.shape[ensax],ninit,-1])
    nsmp = fc.shape[-1]

    CRPS,REL,CRPS_pot,UNC = np.ones([4,nsmp]) #

    N = fc.shape[0]
    p = np.arange(N+1)/N
    
    for ISMP in range(nsmp):
        # compute the necessaryy terms for the decomposition and the CRPS estimator:
        # initialize alpha, beta, o, p and g:
        alpha,beta = np.zeros([2,N+1,ninit])

        alpha[1:-1] = np.diff(np.sort(fc[:,:,ISMP],axis=0),axis=0)
        beta[1:-1] = np.diff(np.sort(fc[:,:,ISMP],axis=0),axis=0)
        
        ob,oa = 0,0
        
        for k in range(ninit):
            x_a = obs[k,ISMP]
            x = np.sort(fc[:,k,ISMP])

            # find the "bin" (two members) that x_a lies in
            # where is the distance smallest?
            ix_close = (abs(x-x_a)).argmin()
            # check whether x_a lies above or below that index to determine "bin"
            if x[ix_close] < x_a:
                obs_bin = ix_close + 1
            elif x[ix_close] > x_a:
                obs_bin = ix_close

            if obs_bin != 0:
                alpha[obs_bin,k] = x_a - x[obs_bin-1]
            else:
                # count instances of x_a lower than lowest ens member
                ob += 1
            if obs_bin != N:
                beta[obs_bin,k] = x[obs_bin] - x_a
            else:
                # count instances of x_a  lower than highest ens member
                oa += 1

            alpha[obs_bin+1:,k] = 0
            beta[:obs_bin,k] = 0
        
        # average over the sample dimension (initializations)
        ALPHA = np.mean(alpha,axis=1); BETA = np.mean(beta,axis=1)
        g,o = np.zeros_like(ALPHA),np.zeros_like(ALPHA)
        g[1:-1] = ALPHA[1:-1] + BETA[1:-1]
        o[1:-1] = BETA[1:-1]/g[1:-1]

        # treat outliers:
        o[0] = ob/ninit
        o[N] = oa/ninit

        if ob != 0:
            g[0] = BETA[0]/o[0]
        if oa != 0:
            g[N] = ALPHA[N]/(1-o[N])
        
        
        # reliability component:
        REL[ISMP] = np.nansum(g*(o-p)**2) # sum over index N ensemble members
        # potential CRPS
        CRPS_pot[ISMP] = np.nansum(g*o*(1-o))
        # these two should add up to give the full CRPS:
        CRPS[ISMP] = np.sum(ALPHA*p**2 + BETA*(1-p)**2)

        # compute uncertainty component (only a function of the observations):
        pk = np.zeros(ninit)
        pk[0] = 0
        for k in range(ninit):
            pk[k] = pk[k-1] + 1/ninit
        UNC[ISMP] = np.sum(pk[1:]*(1-pk[1:])*np.diff(np.sort(obs[:,ISMP])))

    # resolution component as residual:
    RES = UNC - CRPS_pot
    
    return CRPS,REL,RES,UNC


if __name__ == '__main__':

    # can be run from shell inside its own directory using ./evaluate_univpp <path_to_hindcast>

    assert len(sys.argv) > 1, 'pass directory to forecasts as input!'

    dir = sys.argv[1]

    if len(sys.argv) > 2:
        mode = sys.argv[2]
    else:
        mode = 'ensemble'

    for LT in np.arange(1,25):
        VF = verification(find_dir=dir,lead_time=LT,mode=mode)
        VF.full_eval()