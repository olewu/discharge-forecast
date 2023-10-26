import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.interpolate import interp1d
from properscoring import crps_ensemble

class Qmap():
    def __init__(self,observations,hindcast):
        
        self.observations = observations
        self.hindcast = hindcast.flatten()
        # check that observations and forecasts have consistent sizes

        self.hindcast_oshape = hindcast.shape

        return
    
    def fit(self,fill_val=None,interpolation='linear'):

        self.ecdf_obs = ECDF(self.observations)
        # 0-quantile will get value of -inf. For temperature, it is not clear what the lowest possible value is
        # we arbitrarily set it to the minimum - 0.5 degrees
        if fill_val is not None:
            self.ecdf_obs.x[0] = fill_val
        self.ecdf_hc = ECDF(self.hindcast)

        self.obs_quantiles = self.ecdf_obs(self.observations)
        self.hc_quantiles = self.ecdf_hc(self.hindcast)

        self.inverse_ecdf = interp1d(self.ecdf_obs.y,self.ecdf_obs.x,kind=interpolation)

        return
    
    def predict(self,forecast):

        self.forecast = forecast.flatten()
        self.forecast_oshape = forecast.shape

        self.fc_quantiles = self.ecdf_hc(self.forecast)

        self.prediction = self.inverse_ecdf(self.fc_quantiles).reshape(self.forecast_oshape)

        return self.prediction
    
    def evaluate_prediction(self,verification):

        self.verification = verification

        print(
            '\nMean\nobs:\t{0:1.2f}\nraw:\t{1:1.2f}\nqmap:\t{2:1.2f}'.format(
                self.verification.mean(),
                (self.forecast).mean(),
                (self.prediction).mean()
            )
        )
        print(
            '\nVariance\nobs:\t{0:1.2f}\nraw:\t{1:1.2f}\nqmap:\t{2:1.2f}'.format(
                self.verification.var(),
                (self.forecast).var(),
                (self.prediction).var()
            )
        )

        print(
            '\nCRPS\nraw:\t{0:1.3f}\nqmap:\t{1:1.3f}'.format(
                crps_ensemble(self.verification,self.forecast.reshape(self.forecast_oshape).T).mean(),
                crps_ensemble(self.verification,self.prediction.T).mean()
            )
        )

        return