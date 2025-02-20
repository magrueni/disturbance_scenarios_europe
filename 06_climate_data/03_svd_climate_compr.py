# --------------------------------------------------------------
# Climate Data 03: Climate data compression
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# In this script we use the climate compressor dnn to reduce dimensionality
# of the climate data. 
# 
# Notes:
# - An example of the compressed data is stored in /svd/svd_example/climate/
# - The whole dataset is too large to store here but will
# be shared upon request.
# --------------------------------------------------------------



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 09:28:15 2022

@author: mgruenig
"""

import os
import tensorflow as tf
import gc

GPU = '2'
os.environ['CUDA_VISIBLE_DEVICES'] = GPU # Use 3rd GPU
gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)


# read meta data and examples from simulation db
import sqlite3
import pandas as pd
import numpy as np
import glob
import pyarrow.feather
import pyarrow as pa
import pyarrow.parquet as pq
import re
import csv


# load the climate compressor
from tensorflow.keras.models import Sequential, save_model, load_model


filepath = "/.../models/climcompressor_v4.h5"
model_tas = load_model(filepath, compile = True)


# scenarios
scenarios = ["ICHEC-EC-EARTH_rcp_8_5","ICHEC-EC-EARTH_rcp_4_5", "ICHEC-EC-EARTH_rcp_2_6" \
              "NCC-NorESM1-M_rcp_8_5", "NCC-NorESM1-M_rcp_4_5", "NCC-NorESM1-M_rcp_2_6", \
               "MPI-M-MPI-ESM-LR_rcp_2_6", "MPI-M-MPI-ESM-LR_rcp_4_5", "MPI-M-MPI-ESM-LR_rcp_8_5"]
                

timesteps = ["2011-2020" , "2021-2030", "2031-2040", "2041-2050", "2051-2060", \
              "2061-2070", "2071-2080", "2081-2090", "2091-2100"]

# timesteps = ["1981-1990", "1991-2000", "2001-2010"]
# scenarios = ["ICHEC-EC-EARTH_historical", "NCC-NorESM1-M_historical", "MPI-M-MPI-ESM-LR_historical"]

for s in scenarios:
    
    print(s)

    for timestep in timesteps:
            
            print(timestep)
            
            path = "/clim_data/daily_climate/"
            pred_files_all = glob.glob(path + f"all_grids_{timestep}_{s}.dat", recursive = True)
            
            #  read in the file
            fname = pred_files_all[0]
            df = pyarrow.feather.read_feather(fname)
            df = df.sort_values("year")
            
            # separate climate data from other columns and add soil info for climate compressor
            df_base = df[list(df.loc[:, "point_id":"year"])]
            df_aux = df[list(df.loc[:, "bbgen":"summervpd"])]
            # print(df.columns.tolist())
            df_clim = df[list(df.loc[:, "tas_1":"vpd_365"])]
            # add constant values for site conditions - taken from bins npp basic layers: nitrogen, depth, pctSand, pctSilt, pctClay
            df_user_all_examples_site = pd.DataFrame({'nitrogen' : len(df) * [70], 'depth' : len(df) * [90], 'pctSand' : len(df) * [45], \
                                                      'pctSilt' : len(df) * [30], 'pctClay' : len(df) * [25]})
                     
            # the climate data frame is reformated and then fed to the climate compressor        
            df_clim_np = df_clim.to_numpy()
            clim = np.reshape(df_clim_np, (-1, 4, 365))
            scale_clim = np.array( (10, 3, 20, 1) )
            scale_clim = scale_clim[ np.newaxis, :, np.newaxis]
            clim = clim / scale_clim
            clim = np.moveaxis(clim, 1, 2)
            clim.astype(np.float32)
            clim = tf.convert_to_tensor(clim, dtype = tf.float32)
                      
            site = df_user_all_examples_site.to_numpy()
            site = site / [100, 100, 100, 100, 100]
            site.astype(np.float32)
            env = tf.convert_to_tensor(site, dtype = tf.float32)
            dat = clim, env
                      
            # climate compressor
            tas_layer_output = model_tas.predict(dat, verbose = 0)
            tas_df = pd.DataFrame(tas_layer_output[1])
            tas_df = tas_df.set_axis([f'comp_clim_{i+1}' for i in range(14)], axis = 1, copy = False)
            
            npp_df = pd.DataFrame(tas_layer_output[0])
            npp_df = npp_df.set_axis([f'npp_{i+1}' for i in range(10)], axis = 1, copy = False)
            
            df_clim_averages = pd.DataFrame()
            df_clim_averages["MAT"] = df_clim.loc[:, 'tas_1':'tas_365'].mean(axis = 1)
            df_clim_averages["ANP"] = df_clim.loc[:, 'prec_1':'prec_365'].mean(axis = 1)
            df_clim_averages["Mrad"] = df_clim.loc[:, 'rad_1':'rad_365'].mean(axis = 1)
            df_clim_averages["Mvpd"] = df_clim.loc[:, 'vpd_1':'vpd_365'].mean(axis = 1)
            
    
            days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            day_in_year = 0

            for n, days in enumerate(days_in_month):
                df_clim_averages["tas_month" + str(n+1)] = \
                df_clim.iloc[:, day_in_year:day_in_year + days].mean(axis = 1)
                
                day_in_year += days
                
            day_in_year = 0
            
            for n, days in enumerate(days_in_month):
                df_clim_averages["prec_month" + str(n+1)] = \
                df_clim.iloc[:, 365 + day_in_year: 365 + day_in_year + days].mean(axis = 1)
                day_in_year += days
                
            df_clim_averages = df_clim_averages.reset_index(drop=True)
            df_base = df_base.reset_index(drop=True)   
            df_aux = df_aux.reset_index(drop=True)   
         
            # put everything together 
            df_fin = pd.concat([df_base.reset_index(drop=True), df_clim_averages.reset_index(drop=True), \
                                tas_df.reset_index(drop=True), npp_df.reset_index(drop=True), \
                                df_aux.reset_index(drop=True)], axis = 1)
            
            df_fin = df_fin.rename(columns={"point_id": "climateId"})


            if df_fin.isnull().values.any():
                 print("HAS NAN: " + s + " chnk:" + str(timestep))
                 df_fin = df_fin.dropna().reset_index(drop=True)
            
            filename = f"/clim_data/annual_climate/svd_clim_{s}.csv"

            # Check if file exists
            if os.path.isfile(filename):
                # If file exists, append data without writing the header
                df_fin.to_csv(filename, mode='a', header=False, index=False, float_format='%.3f')
            else:
                # If file does not exist, write the header and data
                df_fin.to_csv(filename, mode='w', header=True, index=False, float_format='%.3f')
          
            
            del clim, env, dat, tas_layer_output, df_fin
            gc.collect()
            tf.keras.backend.clear_session()
