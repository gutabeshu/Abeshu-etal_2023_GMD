import os
import spotpy
import numpy as np
import pandas as pd
import scipy.io as scio
import matplotlib.pylab as plt
from scipy.stats import qmc
from itertools import product

from xanthos.runoff.abcd_managed import AbcdManaged
from xanthos.data_reader.data_abcd import DataAbcd
from xanthos.data_reader.data_reference import DataReference

#warning
import warnings
warnings.filterwarnings('ignore')


class Calibrate_runoff:
    def __init__(self,
                 pet,
                 precip,
                 tmin,
                 SM,
                 ts_bsn_obs,
                 basin_ids,
                 basin_idx,
                 nmonths,
                 runoff_spinup,
                 calib_algorithm,
                 params_ro
                 ):

        self.SM = SM
        self.pet = pet
        self.precip = precip
        self.tmin = tmin
        self.ts_bsn_obs = ts_bsn_obs
        self.basin_ids = basin_ids
        self.basin_idx = basin_idx
        self.nmonths = nmonths
        self.runoff_spinup = runoff_spinup
        self.calib_algorithm = calib_algorithm
        self.params_ro = params_ro
        #Runoff
        self.LB = 1e-4
        self.UB = 1 - self.LB	
		
     # set up parameters
    def parameters(self):
        params = self.params_ro        
        return spotpy.parameter.generate(params)
			
    def simulation(self, pars):
        """ABCD model and mrtm routing model : this function provides simulated streamflow"""
        he = AbcdManaged(pars=pars,
                        soil_water_initial=self.SM[self.basin_idx],
                        pet=self.pet[self.basin_idx,:],
                        precip=self.precip[self.basin_idx,:],
                        tmin=self.tmin[self.basin_idx,:],
                        basin_ids=self.basin_ids[self.basin_idx],
                        process_steps=self.nmonths,
                        spinup_steps=self.runoff_spinup,
                        method="dist")					 
        he.emulate()
        self.rsim =  np.nanmean(he.rsim, 1)
        annual_rsim = timeseries_coverter(self.rsim , start_yr=1971, ending_yr=2001)

        return annual_rsim[0:20]
        

    #@staticmethod
    def objectivefunction(self, simulation, evaluation):
        """Calculates Model Performance.
        Objective function to be minimized (if sceua /NSGAII is used) and maximized (all others)
        """
        # sceua requires minimization which will result in a negative KGE
        method = self.calib_algorithm
        if (method == 'sceua') | (method == 'NSGAII'):
            multiplier = -1
        else:
            multiplier = 1
		
        obj1 = spotpy.objectivefunctions.rrmse(evaluation, simulation) #* multiplier

        return obj1


    def evaluation(self):
        """observed streamflow data"""
        annual_obs = timeseries_coverter(self.ts_bsn_obs , start_yr=1971, ending_yr=1990)

        return annual_obs

    def save(self, objectivefunctions, parameter, simulations):
        line = str(objectivefunctions) + ',' + str(parameter).strip('[]') + ',' + str(simulations).strip('[]') + '\n'
        self.database.write(line)

# calibration set up
def calibrate_basin(pet,
                    precip,
                    tmin,
                    SM,
                    ts_bsn_obs,
                    basin_ids,
                    basin_idx,
                    nmonths,
                    runoff_spinup,
                    repetitions,
                    calib_algorithm,
                    dbname_dir,
                    params_runoff):
    runoff_model_spot_setup = Calibrate_runoff(pet,
                                        precip,
                                        tmin,
                                        SM,
                                        ts_bsn_obs,
                                        basin_ids,
                                        basin_idx,
                                        nmonths,
                                        runoff_spinup,
                                        calib_algorithm,
                                        params_runoff,
                                        )
    # parallel ='seq' # Runs everthing in sequential mode
    np.random.seed(2000) # Makes the results reproduceable
    skip_duplicates = True
    if calib_algorithm == 'sceua':	          
        sampler = spotpy.algorithms.sceua(runoff_model_spot_setup,dbname= dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                          
        sampler.sample(repetitions, ngs=10, kstop=10, peps=1e-7, pcento=1e-7)

    elif calib_algorithm == 'NSGAII':	    
        n_pop = 10
        repetitions_nsgaii = int(repetitions / n_pop)         
        sampler = spotpy.algorithms.NSGAII(runoff_model_spot_setup, dbname=dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                                
        sampler.sample(repetitions_nsgaii, n_obj= 1, n_pop = n_pop)

    elif calib_algorithm == 'mcmc':	          
        sampler = spotpy.algorithms.mcmc(runoff_model_spot_setup,dbname=dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                          
        sampler.sample(repetitions)


    elif calib_algorithm == 'demcz':	          
        sampler = spotpy.algorithms.demcz(runoff_model_spot_setup,dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                            
        sampler.sample(repetitions)

    elif calib_algorithm == 'dream':	          
        sampler = spotpy.algorithms.dream(runoff_model_spot_setup,dbname=dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                            
        sampler.sample(repetitions)
    elif calib_algorithm == 'abc':	          
        sampler = spotpy.algorithms.abc(runoff_model_spot_setup,dbname=dbname_dir + calib_algorithm + '_Runoff',
                                    dbformat="csv",dbappend=False,save_sim=False)#,
                                    #parallel='mpi' )                                            
        sampler.sample(repetitions)

    results = pd.read_csv(dbname_dir + calib_algorithm + '_Runoff' + '.csv')
    results_sorted = results.sort_values(by = 'like1', ascending=True).reset_index(drop=True)
    ro_params_selected = np.array(results_sorted.loc[0:100][['para',	'parb',	'parc',	'pard',	'parm']])

    return ro_params_selected
    



# converts monthly to annual
def timeseries_coverter(data_array, start_yr, ending_yr):
    from datetime import date, timedelta
    sdate = date(start_yr,1,1)
    edate = date(ending_yr, 12, 31)  
    data_ts = pd.DataFrame(data_array)
    
    data_ts.index = pd.date_range(start=sdate, end=edate, freq='M')
    mean_annual_data = np.squeeze(np.array(data_ts.resample('A').sum()))

    return mean_annual_data