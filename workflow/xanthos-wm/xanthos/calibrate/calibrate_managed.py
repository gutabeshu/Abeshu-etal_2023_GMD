"""
Calibrate the ABCD model with and without water management.
@authors:Guta Abeshu : gwabeshu@uh.edu  
         University of Houston
         HongYi Li : hli57@uh.edu,
         University of Houston
License:  BSD 2-Clause, see LICENSE and DISCLAIMER files
Copyright (c) 2018, Battelle Memorial Institute
"""

import os
import spotpy
import numpy as np
import pandas as pd
import scipy.io as scio
import matplotlib.pylab as plt
from scipy.stats import qmc
from itertools import product

from xanthos.runoff.abcd_managed import AbcdManaged
from xanthos.data_reader.data_mrtm_managed import DataMrtmManaged
from xanthos.data_reader.data_calibration_managed import DataCalibrationManaged
from xanthos.data_reader.data_reference import DataReference
from xanthos.data_reader.data_abcd import DataAbcd
import xanthos.calibrate.calibrate_abcdm_firststage as calibrate_runoff

import xanthos.routing.mrtm_managed as routing_mod
import xanthos.utils.general as helper
import xanthos.utils.math as umth

#warning
import warnings
warnings.filterwarnings('ignore')

class CalibrateManaged:
    """Calibrate the ABCD model with and without water management.
    :params repetitions:                    maximum number of function evaluations allowed during optimization
    :type repetitions:                      int
    """

    ncell = 67420
    ngridrow = 360
    ngridcol = 720
     
	#Runoff
    LB = 1e-4
    UB = 1 - LB
    
	#WM
    LBwm1 = 0.1
    UBwm1 = 10.0
    LBwm2 = 0.5
    UBwm2 = 1
    PARAMETER_NAMES = ['a', 'b', 'd', 'm', 'beta', 'alpha']

    def __init__(self,
                 basin_num,
                 basin_ids,
                 basin_areas,
                 precip,
                 pet,
                 obs,
                 tmin,
                 runoff_spinup,
                 set_calibrate,
                 obs_unit,
                 out_dir,
                 nmonths=None,
                 router_func=None,
                 config_obj=None,
                 cal_observed=None,
                 purpose_file=None,
                 capacity_file=None,
                 hp_release_file=None,
                 water_consumption_file=None,
                 instream_flow_natural_file=None,
                 initial_channel_storage_natural_file=None,
                 sm_file=None,
                 mtif_natural_file=None,
                 maxtif_natural_file=None,
                 total_demand_cumecs_file=None,
                 grdc_coord_index_file=None,
                 start_year=None,
                 end_year=None,
                 repetitions=None,
                 calib_algorithm=None,
                 flow_distance_file=None,
                 flow_direction_file=None,
                 stream_velocity_file=None,
                 historical_mode="True",
                 hist_channel_storage_file=None,
                 hist_channel_storage_varname=None,
                 routing_spinup=None):
        """Initialize calibration data and parameters.
        :param basin_num:      basin number as an integer
        :param basin_ids:      an array of basin ids per grid cell that are associated with the basin
        :param basin_areas:    an array of basin areas per grid cell that are associated with the basin
        :param precip:         precipitation in mm/month
        :param pet:            PET in mm/month
        :param obs:            Observed runoff in mm/month
        :param tmin:           minimum temperature in degrees C
        :param runoff_spinup:  the number of months from time 0 in the input data to use as spin up
        :param set_calibrate:  0 to calibrate to observed runoff, 1 to calibrate to observed streamflow
        :param obs_unit:       the unit of the input data
        :param out_dir:        calibrated parameters output directory
        :param router_func:    objective function for calibrating routing
        """

        if config_obj is None:
            self.start_year = start_year
            self.end_year = end_year
            
            if nmonths is None:
                self.nmonths = (self.end_year - self.start_year + 1) * 12
            else:
                self.nmonths = nmonths

            # load reference data
            self.reference_data = DataReference(nmonths=self.nmonths)

            self.flow_distance_file = flow_distance_file
            self.flow_direction_file = flow_direction_file
            self.stream_velocity_file = stream_velocity_file
            self.historical_mode = historical_mode
            self.hist_channel_storage_file = hist_channel_storage_file
            self.hist_channel_storage_varname = hist_channel_storage_varname
            self.routing_spinup = routing_spinup
            self.repetitions = repetitions
            self.calib_algorithm = calib_algorithm

        else:
            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths
            self.reference_data = DataReference(config=config_obj)
            self.routing_spinup = config_obj.routing_spinup
            self.repetitions = config_obj.repetitions
            self.calib_algorithm = config_obj.cal_algorithm

            self.flow_distance_file = config_obj.flow_distance
            self.flow_direction_file = config_obj.flow_direction
            self.stream_velocity_file = config_obj.strm_veloc
            self.historical_mode = config_obj.HistFlag

            self.hist_channel_storage_file = config_obj.ChStorageFile
            self.hist_channel_storage_varname = config_obj.ChStorageVarName

        self.basin_num = basin_num
        self.basin_ids = basin_ids
        self.basin_areas = basin_areas
        self.precip = precip
        self.pet = pet
        self.obs = obs	
        self.tmin = tmin
        self.runoff_spinup = runoff_spinup
        self.router_func = router_func
        self.set_calibrate = set_calibrate
        self.obs_unit = obs_unit
        self.out_dir = out_dir


        self.parameternames = CalibrateManaged.PARAMETER_NAMES

        # load calibration data
        self.calib_data = DataCalibrationManaged(config_obj=config_obj,
                                                 cal_observed=cal_observed,
                                                 purpose_file=purpose_file,
                                                 capacity_file=capacity_file,
                                                 hp_release_file=hp_release_file,
                                                 water_consumption_file=water_consumption_file,
                                                 instream_flow_natural_file=instream_flow_natural_file,
                                                 initial_channel_storage_natural_file=initial_channel_storage_natural_file,
                                                 sm_file=sm_file,
                                                 mtif_natural_file=mtif_natural_file,
                                                 maxtif_natural_file=maxtif_natural_file,
                                                 total_demand_cumecs_file=total_demand_cumecs_file,
                                                 grdc_coord_index_file=grdc_coord_index_file,
                                                 start_year=self.start_year,
                                                 end_year=self.end_year)

        # index for gauge station locations
        self.grdcData_info  = self.calib_data.grdc_coord_index_file
        self.basin_grdc_indx = np.where(self.grdcData_info['basins'].values == self.basin_num)[0][0]
        self.grdc_xanthosID = self.grdcData_info.loc[self.basin_grdc_indx]['xanthosID_new_adj'] - 1
        self.grdc_drainage_area = self.grdcData_info.loc[self.basin_grdc_indx]['area']*1e6 #m2

        # routing inputs: wdirr, irrmean, tifl, ppose, cpa,dscells
        self.wdirr = np.copy(self.calib_data.total_demand_mmpermonth)
        self.irrmean = np.mean(self.wdirr, axis=1)  # mean demand
        self.ppose = np.array(self.calib_data.purpose)
        self.cpa  = np.array(self.calib_data.capacity) * 1e6  # m3
        self.installed_cap = np.array(self.calib_data.installed_cap) # MW	
        self.q_max= np.array(self.calib_data.maxtif_natural) # m3/s	
        self.surface_area = np.array(self.calib_data.surface_area) # square km	
        self.max_depth = np.array(self.calib_data.max_depth) # square km
		
        self.WConsumption = np.array(self.calib_data.water_consumption)
        self.chs_ini =np.array(self.calib_data.ini_channel_storage)
        self.Initial_instream_flow = np.array(self.calib_data.instream_flow_natural)
        self.SM = np.squeeze(np.array(self.calib_data.sm))
        self.res_data = pd.DataFrame([self.cpa *1e-6, 
		                              self.installed_cap, 
		                              self.q_max, 
		                              self.surface_area, 
		                              self.max_depth]).transpose()
        self.res_data.columns = ["CAP",	"ECAP", "FLOW_M3S", "AREA", "DAM_HGT"]							   
        # routing data
        self.yr_imth_dys = helper.set_month_arrays(self.nmonths, self.start_year, self.end_year)
        self.map_index = umth.sub2ind([CalibrateManaged.ngridrow, CalibrateManaged.ngridcol],
                                      self.reference_data.coords[:, 4].astype(int) - 1,
                                      self.reference_data.coords[:, 3].astype(int) - 1)

        self.routing_data = DataMrtmManaged(start_year=self.start_year,
                                            end_year=self.end_year,
                                            flow_distance_file=self.flow_distance_file,
                                            flow_direction_file=self.flow_direction_file,
                                            stream_velocity_file=self.stream_velocity_file,
                                            historical_mode=self.historical_mode,
                                            hist_channel_storage_file=self.hist_channel_storage_file,
                                            hist_channel_storage_varname=self.hist_channel_storage_varname)

 		
        ## contributing grids for grdc  station
        self.dsid = routing_mod.downstream(self.reference_data.coords, self.routing_data.flow_dir, CalibrateManaged)
        self.upid = routing_mod.upstream(self.reference_data.coords, self.dsid, CalibrateManaged)
        self.basin_idx = routing_mod.grdc_stations_upstreamgrids(self.upid, self.grdc_xanthosID) # grids upstream of grdc
        # upstream genmatrix matrix 			
        self.um, self.up = routing_mod.upstream_genmatrix(self.upid) 

        #transpose data for use in the ABCD model
        #self.basin_idx = np.where(self.basin_ids == self.basin_num)[0]
        self.bsn_areas = self.basin_areas[self.basin_idx]        
        self.bsn_PET = self.pet[self.basin_idx]
        self.bsn_P = self.precip[self.basin_idx]
        # initial soil moisture
        self.bsn_SM = self.SM[self.basin_idx]

		
        
        # Minimum temperature is optional; if not provided, the snow components
        # of the model is effectively removed, so remove the model parameter
        # if no tmin provided, just ensure it is larger than the rain threshold
        self.nosnow = self.tmin is None
        if self.nosnow:
            self.bsn_TMIN = None
            self.tminx = None
        else:
            self.bsn_TMIN = self.tmin[self.basin_idx]
            self.tminx = self.tmin

        # Unit conversion for runoff case
        #if self.obs_unit == "km3_per_mth":
        self.conversion = self.bsn_areas * 1e-6
        ##  Observation data for calibration 	
        self.indx_ro = 0		
        if self.set_calibrate == 0:
            self.bsn_vic_runoff = np.squeeze(self.obs[np.where(self.obs[:,0]==self.basin_num)[0],3])
        else:
            self.bsn_obs = self.obs[np.where(self.obs[:, 0] == self.basin_num)][: self.nmonths, 1]
            # calibration and validation data
            self.bsn_Robs_calib = self.bsn_obs[0:120]
            self.bsn_Robs_valid = self.bsn_obs[121:240]

            # basin observed runoff 
            self.bsn_obs_runoff = np.multiply(self.bsn_obs, self.yr_imth_dys[0:240, 2])*(1e3*24*3600) / self.grdc_drainage_area  #mm/month        

        
        # residence time in hr
        Lst = self.routing_data.flow_dist
        self.Vst = self.routing_data.str_velocity
        self.Vst[self.Vst < 0.01] = 0.01
        grid_size = np.sqrt(self.reference_data.area) * 1000
        nn = np.where(Lst < grid_size)[0]
        Lst[nn] = grid_size[nn]
        self.Tr = np.divide(Lst, self.Vst) / 3600
        self.beta_local = 0.95*self.Tr/3   

        self.flow_dist = self.routing_data.flow_dist
        grid_size = np.sqrt(self.reference_data.area) * 1e3 # grid size in m
        nn_grids = np.where(self.flow_dist < grid_size)[0]
        self.flow_dist[nn_grids] = grid_size[nn_grids] # update flow distance
        self.flow_dir = self.routing_data.flow_dir

        # routing initializations
        self.chs_prev = np.squeeze(self.chs_ini)  # initial channel storage in m3
        self.instream_flow = np.squeeze(self.Initial_instream_flow)  # initial channel flow in m3 /s
        self.mtifl_natural= np.squeeze(self.Initial_instream_flow)
		
       # reservoir params
        usgrid_ppose = np.squeeze(np.array(pd.DataFrame(self.ppose[self.basin_idx])))
        grdc_us_grids_ = np.squeeze(np.array(pd.DataFrame(self.basin_idx)))
        self.HPindex= grdc_us_grids_[(usgrid_ppose==1)] # Hydropower reservoir cells	
        self.us_resrv_idx= grdc_us_grids_[(usgrid_ppose>0)]       

        self.Main_UseHP = (self.ppose==1) & (self.basin_ids == self.basin_num)
        self.halfDegree_global_dfG = self.res_data[self.Main_UseHP].reset_index(drop=True)  
        self.Main_UseHPXP = np.where((self.ppose==1) & (self.basin_ids == self.basin_num))[0]

        # best parameter from preceeding simulations
        self.best_params = None
        self.initial_cond = 0 # to run reservoir initialization the first time only
        self.dir_storage = self.out_dir + '/reservoir_initial_storage/Sini_'

        # calibration result output path
        self.ModelPerformance = os.path.join(self.out_dir, f"basin_calibration_{self.basin_num}")

        # set up parameters  for first stage or runoff  
        # set number of parameter combinations with latin hyper cube sampling
        l_bounds = [CalibrateManaged.LB, CalibrateManaged.LB,   CalibrateManaged.LB, CalibrateManaged.LB, CalibrateManaged.LB]
        u_bounds = [CalibrateManaged.UB, 8- CalibrateManaged.LB,CalibrateManaged.UB, CalibrateManaged.UB, CalibrateManaged.UB]                    
        sampler_lhc = qmc.LatinHypercube(d=5, seed=42)
        sample_params = sampler_lhc.random(n=10000)
        self.sample_params_set = qmc.scale(sample_params, l_bounds, u_bounds)
        self.params_ro = [spotpy.parameter.List('a',list(self.sample_params_set[:,0])),        
                            spotpy.parameter.List('b',list(self.sample_params_set[:,1])),  
                            spotpy.parameter.List('c',list(self.sample_params_set[:,2])),  
                            spotpy.parameter.List('d',list(self.sample_params_set[:,3])),  
                            spotpy.parameter.List('m',list(self.sample_params_set[:,4])),             
                            ]         
        # two stage calibration case
        if self.set_calibrate == 1:	        
            print("\tStarting The First Stage Parameter Selection : Runoff Parameters Selection")
            self.ro_params = calibrate_runoff.calibrate_basin(self.pet,
                                                            self.precip,
                                                            self.tmin,
                                                            self.SM,
                                                            self.bsn_obs_runoff,
                                                            self.basin_ids,
                                                            self.basin_idx,
                                                            self.nmonths,
                                                            self.runoff_spinup,
                                                            self.repetitions,
                                                            self.calib_algorithm,
                                                            self.ModelPerformance,
                                                            self.params_ro)
    


            
            # list of parameters values for second stage calibration
            print("\tStarting The Second Stage Parameter Selection: Runoff + Routing Parameter Selection")   
            # parameter setup for second stage calibration
            wmp_beta  = [0.1 ,0.2 ,0.3 ,0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,        
                            1 ,2 ,3 ,4, 5 , 6, 7, 8, 9, 10] #Give possible beta values as a List
            wmp_alpha = [0.85]                              #Give possible alpha values as a List
            wm_params = pd.DataFrame(product(wmp_beta, wmp_alpha))
            for ii in range(self.ro_params.shape[0] + 1):
                if ii==0:
                    abcdm_parameters = np.mean(self.ro_params,0)  
                else:
                    abcdm_parameters = self.ro_params[ii-1,:]  
                wm_params[['a','b','c','d','m']] =  np.array(abcdm_parameters)  
                if (ii==0):
                    self.wm_abcdm_parameters = wm_params.copy()
                else:
                    self.wm_abcdm_parameters = pd.concat([self.wm_abcdm_parameters, wm_params], 0).reset_index(drop=True)
                
            params_all = np.array(self.wm_abcdm_parameters)
            self.params_wm = [spotpy.parameter.List('a',list(params_all[:,2])),        
                                spotpy.parameter.List('b',list(params_all[:,3])),  
                                spotpy.parameter.List('c',list(params_all[:,4])),  
                                spotpy.parameter.List('d',list(params_all[:,5])),  
                                spotpy.parameter.List('m',list(params_all[:,6])), 
                                spotpy.parameter.List('wmbeta', list(params_all[:,0])),
                                spotpy.parameter.List('wmalpha',list(params_all[:,1])) ]    

    # parameter set up
    def parameters(self):
        '''Returns ABCD Params'''
        if self.set_calibrate == 0:
            params = self.params_ro
        else:       
            params = self.params_wm
            
        return spotpy.parameter.generate(params)

    # model simulation set up			
    def simulation(self, pars):
        """ABCD model and mrtm routing model : this function provides simulated streamflow"""
        if self.set_calibrate == 0:	
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
            self.rsim =  np.nansum(he.rsim * self.conversion, 1)

            return self.rsim
			
        else:
            # runoff model
            he = AbcdManaged(pars=pars,
							soil_water_initial=self.SM,
							pet=self.pet,
							precip=self.precip,
							tmin=self.tmin,
							basin_ids=self.basin_ids,
							process_steps=self.nmonths,
							spinup_steps=self.runoff_spinup,
							method="dist")					 
            he.emulate()
			
            # runoff
            self.runoff = he.rsim.T

            # load routing data
            beta_basin = np.copy(self.beta_local)			
            beta_basin[(self.Tr/pars[5]) > 3] = pars[5]	
                       
            # adjusted velocity
            self.str_velocity = np.multiply(np.copy(self.Vst), beta_basin)
      
			
            # Preallocation for variables to be returned from routing module
            self.Avg_ChFlow = np.zeros_like(self.precip)    # average channel flow m3/s
            self.ChStorage = np.zeros_like(self.precip)     # channel storage   m3
            self.ResStorage = np.zeros_like(self.precip)    # reservoir storage m3
            self.Qout_res_avg = np.zeros_like(self.precip)  # reservoir outflow m3/s
            self.Qin_res_avg = np.zeros_like(self.precip)   # reservoir inflow  m3/s
            self.HP_Release = np.zeros([len(self.basin_ids),1001,12]) #1001 = S discretization & 12 = nmonths/year
            

            self.routing_timestep = 3*3600  # seconds	
            for ii in range(1):
                if ii==0:
                    ###################### WITHOUT RESERVOIR ###########################           
                    # Reservoir flag
                    self.res_flag = 0  # 1 if with reservoir, 0 without reservoir
                    # routing time step               
                    self.routing_length = 120#self.routing_spinup		      
                    #place holders
                    self.Sini = np.zeros_like(self.mtifl_natural)  # initial reservoir storage at time t in m3
                    self.res_prev = np.zeros_like(self.mtifl_natural)  # reservoir storage at beggining of the year in m3           
                else:
                    ###################### WITH RESERVOIR ###########################            
                    ## Reservoir flag
                    self.res_flag = 1  # 1 if with reservoir, 0 without reservoir	                     
                    ## initial values
            
                    self.mtifl_natural = np.mean(self.Avg_ChFlow,1)
                    if ii == 1 :      
                        # routing time step               
                        self.routing_length = self.routing_spinup                        
                        self.Sini = 0.5 * np.squeeze(self.cpa)  # initial reservoir storage at time t in m3
                        self.res_prev = 0.5 * np.squeeze(self.cpa)  # reservoir storage at beggining of the year in m3
                    else: 
                        # routing time step               
                        self.routing_length = 240                      
                        Sini_read = np.load(self.dir_storage + str(self.basin_num) +'.npy')
                        self.Sini =  np.squeeze(Sini_read)  # initial reservoir storage at time t in m3
                        self.res_prev = np.squeeze(Sini_read)  # reservoir storage at beggining of the year in m3
                        
                    #Hydropower reservoirs operation policy                        
                    if len(self.HPindex) > 0:		
                        # hydropower_res = np.sort(self.HPindex)		
                        # self.HP_Release[hydropower_res,:,:] = np.array([routing_mod.optimal_release_policy(
                        #                                                 self.res_data, self.Avg_ChFlow, res, pars[6]) 
                        #                                                 for res in hydropower_res])	
                        # np.save('D:/XanthosDev/Figures/RunResult-02022022/HP_Release' + str(self.basin_num)+ '.npy',
                        #          self.HP_Release[hydropower_res,:,:])
                        dataflowD = self.Avg_ChFlow[self.Main_UseHPXP, :]
                        self.HP_Release = np.zeros([dataflowD.shape[0],1001,12]) #1001 = S discretization & 12 = nmonths/year
                        for res in range(dataflowD.shape[0]):
                            self.HP_Release[res,:,:] = routing_mod.release_functions(res, dataflowD, self.halfDegree_global_dfG, alpha=0.85)


                ###################### ROUTING ###########################            
                for nm in range(self.routing_length):
                    sr = routing_mod.streamrouting(self.flow_dist,#L
                                                self.chs_prev,#S0
                                                self.instream_flow,#F0
                                                self.str_velocity,#ChV
                                                self.runoff[:, nm],#q
                                                self.basin_areas,#area
                                                self.yr_imth_dys[nm, 2],#nday
                                                self.routing_timestep,#dt
                                                self.um,#um
                                                self.up,#up
                                                self.Sini,#Sini
                                                self.wdirr[:, nm],#wdirr
                                                self.irrmean,#irrmean
                                                self.mtifl_natural,#mtifl
                                                self.ppose,#ppose
                                                self.cpa,#cpa
                                                self.HP_Release[:,  :, np.mod(nm, 12)],#self.HP_Release
                                                self.q_max,#maxTurbineFlow
                                                self.WConsumption[:, nm],#WConsumption
                                                pars[6]	,#alpha
                                                self.res_prev,#Sini_resv
                                                self.res_flag,
                                                self.basin_idx) #res_flag

                    (self.ChStorage[:, nm],
                    self.Avg_ChFlow[:, nm],
                    self.instream_flow,
                    self.Qin_Channel_avg,
                    self.Qout_channel_avg,
                    self.Qin_res_avg,
                    self.Qout_res_avg,
                    self.ResStorage[:, nm]) = sr


                    
                    # update data
                    #self.chs_prev = self.ChStorage[:, nm]
                    self.res_prev = self.ResStorage[:, nm]
                    # update the reservoir storage at beginning of year
                    if np.mod(nm, 12) == 11:
                        self.Sini = self.ResStorage[:, nm]

                    # update data channels storage with reservoir effect
                    self.DsDt_channel = (self.Qin_res_avg - self.Qout_res_avg)*self.yr_imth_dys[nm, 2]*24*3600
                    self.chs_prev = self.ChStorage[:, nm] + self.DsDt_channel
                    self.chs_prev[self.chs_prev < 0] = 0

                # storage out if first run
                if ii == 1 :
                    np.save(self.dir_storage + str(self.basin_num) +'.npy',self.res_prev)
                    self.initial_cond += 1
                    np.save('D:/XanthosDev/Figures/RunResult-02022022/Avg_ChFlow.npy', self.Avg_ChFlow[:, 0:240])

        #np.save('D:/XanthosDev/Figures/RunResult-02022022/CH_OutFlow' + str(self.basin_num)+ '.npy', self.Qout_channel_avg[self.us_resrv_idx, 0:240])
        #np.save('D:/XanthosDev/Figures/RunResult-02022022/ResStorage_' + str(self.basin_num)+ '.npy', self.ResStorage[self.us_resrv_idx, 0:240])
        return self.Avg_ChFlow[self.grdc_xanthosID, 0:120]#, list(self.Avg_ChFlow[self.grdc_xanthosID, 121:240])]


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
				
        obj1 = spotpy.objectivefunctions.kge(evaluation, simulation) * multiplier

 
        return obj1


    def evaluation(self):
        """observed streamflow data"""
        if self.set_calibrate==1:
            self.obs_eval = self.bsn_obs[0:120]
        else:
            self.obs_eval = self.bsn_vic_runoff

        return self.obs_eval

    def save(self, objectivefunctions, parameter, simulations):
        line = str(objectivefunctions) + ',' + str(parameter).strip('[]') + ',' + str(simulations).strip('[]') + '\n'
        self.database.write(line)

    # calibration set up
    def calibrate_basin(self):
        """This function is to calibrate the distributed ABCD + water management model against the GRDC to
        obtain optimized parameters of ABCD(a, b, c, d, m) and Water management (beta and c)
        """
        # parallel ='seq' # Runs everthing in sequential mode
        np.random.seed(2000) # Makes the results reproduceable
        skip_duplicates = True
        if self.set_calibrate == 0:
            if self.calib_algorithm == 'sceua':	          
                sampler = spotpy.algorithms.sceua(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow' ,
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                          
                sampler.sample(self.repetitions, ngs=10, kstop=10, peps=1e-7, pcento=1e-7)

            elif self.calib_algorithm == 'NSGAII':	    
                n_pop = 10
                self.repetitions_nsgaii = int(self.repetitions / n_pop)         
                sampler = spotpy.algorithms.NSGAII(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                                
                sampler.sample(self.repetitions_nsgaii, n_obj= 1, n_pop = n_pop)

            elif self.calib_algorithm == 'mcmc':	          
                sampler = spotpy.algorithms.mcmc(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                          
                sampler.sample(self.repetitions)


            elif self.calib_algorithm == 'demcz':	          
                sampler = spotpy.algorithms.demcz(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                            
                sampler.sample(self.repetitions)

            elif self.calib_algorithm == 'dream':	          
                sampler = spotpy.algorithms.dream(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                            
                sampler.sample(self.repetitions)
            elif self.calib_algorithm == 'abc':	          
                sampler = spotpy.algorithms.abc(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                            
                sampler.sample(self.repetitions)


        elif self.set_calibrate == 1:    
            self.repetitions = self.wm_abcdm_parameters.shape[0]  
            if self.calib_algorithm == 'sceua':	          
                sampler = spotpy.algorithms.sceua(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                          
                sampler.sample(self.repetitions, ngs=10, kstop=10, peps=1e-7, pcento=1e-7)

            elif self.calib_algorithm == 'NSGAII':	                 
                n_pop = 10                
                self.repetitions_nsgaii = int(self.repetitions / n_pop)         
                sampler = spotpy.algorithms.NSGAII(self,dbname=self.ModelPerformance + self.calib_algorithm + '_Flow',
                                            dbformat="csv",dbappend=False,save_sim=False)#,
                                            #parallel='mpi' )                                                
                sampler.sample(self.repetitions_nsgaii, n_obj= 1, n_pop = n_pop)

    def calibration_run(self, x):
        qsim_cal = self.simulation(x)
        qsimulated = self.sim_with_vald

        if self.set_calibrate == 0:
            # KGE of the calibration period
            kge_cal = spotpy.objectivefunctions.kge(self.bsn_obs[0:240], self.sim_with_vald[0:240])
            # KGE of the validation period
            kge_val = spotpy.objectivefunctions.kge(self.bsn_obs[241:480], self.sim_with_vald[241:480])
            print("Calibration KGE: {}, Validation KGE: {}".format(kge_cal, kge_val))                   
            ## NSE of the calibration period
            nse_cal = spotpy.objectivefunctions.nashsutcliffe(self.bsn_obs[0:240], self.sim_with_vald[0:240])
            ## NSE of the validation period
            nse_val = spotpy.objectivefunctions.nashsutcliffe(self.bsn_obs[241:480], self.sim_with_vald[241:480])
            print("Calibration NSE: {}, Validation NSE: {}".format(nse_cal, nse_val))
            
            
        else:
            # KGE of the calibration period
            kge_cal = spotpy.objectivefunctions.kge(self.bsn_obs[0:120], self.sim_with_vald[0:120])
            # KGE of the validation period
            kge_val = spotpy.objectivefunctions.kge(self.bsn_obs[121:240], self.sim_with_vald[121:240])
            print("Calibration KGE: {}, Validation KGE: {}".format(kge_cal, kge_val))                   
            ## NSE of the calibration period
            nse_cal = spotpy.objectivefunctions.nashsutcliffe(self.bsn_obs[0:120], self.sim_with_vald[0:120])
            ## NSE of the validation period
            nse_val = spotpy.objectivefunctions.nashsutcliffe(self.bsn_obs[121:240], self.sim_with_vald[121:240])
            print("Calibration NSE: {}, Validation NSE: {}".format(nse_cal, nse_val))


  
        return qsimulated, self.bsn_obs, kge_cal, kge_val


def process_basin(config_obj, calibration_data, pet, router_function=None):
    """Process single basin."""

    # load ABCD runoff module data
    data_abcd = DataAbcd(config=config_obj)

    cal = CalibrateManaged(config_obj=config_obj,
                           basin_num=config_obj.basin_num,
                           set_calibrate=config_obj.set_calibrate,
                           obs_unit=config_obj.obs_unit,
                           basin_ids=calibration_data.basin_ids,
                           basin_areas=calibration_data.area,
                           precip=data_abcd.precip,
                           pet=pet,
                           obs=calibration_data.cal_obs,
                           tmin=data_abcd.tmin,
                           nmonths=config_obj.nmonths,
                           runoff_spinup=config_obj.runoff_spinup,
                           router_func=router_function,
                           out_dir=config_obj.calib_out_dir,
                           start_year=config_obj.StartYear,
                           end_year=config_obj.EndYear
                           )
    cal.calibrate_basin()


def calibrate_all(settings, calibration_data, pet, router_fn=None):
    """Run calibration for ABCD model for a basins."""

    print("\tCalibrating Basin:  {}".format(settings.basin_num))

    process_basin(settings, calibration_data, pet, router_function=router_fn)

			
def plot_kge(calibration_result_file, output_file_name, dpi=300, figsize=(9, 5)):
    """Plot the KGE result of a calibrated basin"""

    split_extension = os.path.splitext(calibration_result_file)
    if len(split_extension) > 1:
        f = split_extension[0]
    else:
        f = calibration_result_file

    # load results
    results = spotpy.analyser.load_csv_results(f)

    # create plot
    fig = plt.figure(1, figsize=figsize)
    plt.plot(results['like1'] * -1)
    plt.ylabel('KGE')
    plt.xlabel('Repetition')
    fig.savefig(output_file_name, dpi=dpi)

    return plt
	
 
# converts monthly to annual
def timeseries_coverter(data_array, start_yr, ending_yr):
    from datetime import date, timedelta
    sdate = date(start_yr,1,1)
    edate = date(ending_yr, 12, 31)  
    data_ts = pd.DataFrame(data_array)
    
    data_ts.index = pd.date_range(start=sdate, end=edate, freq='M')
    mean_annual_data = np.squeeze(np.array(data_ts.resample('A').sum()))

    return mean_annual_data