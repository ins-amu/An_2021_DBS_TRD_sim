import numpy as np
from tvb.simulator.lab import *
from tvb.datatypes.cortex import Cortex
from tvb.datatypes.region_mapping import RegionMapping
from tvb.datatypes.local_connectivity import LocalConnectivity
from tvb.datatypes.projections import ProjectionMatrix, ProjectionSurfaceEEG
from tvb.datatypes.sensors import SensorsEEG
import os.path
from tvb.simulator.plot.tools import *
from scipy import signal
LOG = get_logger('demo')
import scipy.io
from scipy.sparse import csc_matrix
from copy import copy, deepcopy



def tic():    
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"


# neural mass model - generic 2D oscillator

def osc_gen2D_diff(v,w,c_str,Stim):
    
    tau=4.0; a=-1.5; b=-15.0; c=0.0; d=0.015; e=3.0; f=1.0; g=0.0;
    v_diff=(d*tau*(-f*v**3 + e*v**2 + g*v + w + c_str + Stim))
    w_diff=(d/tau*(c*v**2 + b*v - w + a))

    return v_diff, w_diff

def osc_gen_2D(v,w,c_str,Stim,h,eta_v,eta_w):
    
    (v_diff1,w_diff1)=osc_gen2D_diff_v2(v,w,c_str,Stim)
    v_diff_t=v + h*v_diff1 + eta_v
    w_diff_t=w + h*w_diff1 + eta_w
    (v_diff2,w_diff2)=osc_gen2D_diff_v2(v_diff_t,w_diff_t,c_str,Stim)
    v_out=v + h/2.0 * (v_diff1 + v_diff2) + eta_v
    w_out=w + h/2.0 * (w_diff1 + w_diff2) + eta_w
    
    return v_out, w_out


# connectivity - transmission delay

str_tck=scipy.io.loadmat("patient folder path.../tract_lengths.mat")
str_con_speed=6.0
str_delay=np.round(str_tck['leng_st_save_11']/str_con_speed) 
str_delay_rnd=np.round(str_delay)*5.0 # 1ms = 5 time steps 

# stimulation

stim_tck_mat=loadtxt("patient folder path.../stim_mat_R0.txt")
sigma=0.5
stim_mat=sigma*stim_tck_mat


con_sw=0.0045 # global coupling parameter G
con_sw_num=45


for time in range(0,5): # for parameter sweep for G
    
    # connectivity - connection strength
    
    str_con=scipy.io.loadmat("patient folder path.../weights.mat")
    str_con_sc=con_sw*(str_con['wei_st_save_11'])
    
    # preparation of simulation
    
    sim_leng=300 #ms
    h=0.2;
    node_num=np.size(str_con_sc,0)

    vtg_node_v = (np.random.uniform(-0.098, -0.098, size=node_num)).reshape(node_num,1) 
    vtg_node_w = (np.random.uniform(-0.03, -0.03, size=node_num)).reshape(node_num,1) 
    vtg_node_local = vtg_node_v 

    mean = 0.0; std = 0.0; # noise parameter

    vtg_node_int=np.zeros((len(vtg_node_v),1))
    vtg_node_int2=np.zeros((len(vtg_node_v),1))
    
    # simulation
    # tic()
    
    time2=0; time3=0;
    for time in range(0,sim_leng*10):

        # noise
        eta_v = (np.random.normal(mean, std, size=len(vtg_node_v))).reshape(len(vtg_node_v),1)
        eta_w = (np.random.normal(mean, std, size=len(vtg_node_v))).reshape(len(vtg_node_v),1)

        # stimulus  
        if time >= 50*10 and time < (50*10 + np.size(stim_mat,1)): 
            stim=1.0*stim_mat[:,time3].reshape(len(vtg_node_v),1)
            time3=time3+1
        else:
            stim=np.zeros((len(vtg_node_v),1))           

        # coupling input  
        if time > 350:
            d_ind1 = deepcopy(str_delay_rnd)
            d_ind1.data = -1.0*(d_ind1.data)-1.0
            str_field_delay = deepcopy(d_ind1)
            str_field_delay.data = vtg_node_int[d_ind1.indices,map(int,d_ind1.data)]
            str_field = (asarray(sum(str_field_delay.multiply(str_con_sc),0)).reshape(node_num,1))           
        else:
            str_field=np.zeros((node_num,1))    

        # calculation of v, w of each node  
        (vtg_node_v,vtg_node_w) = osc_gen_2D(vtg_node_v,vtg_node_w,str_field,stim,h,eta_v,eta_w)        
        vtg_node_local = vtg_node_v + 0.098
        
        vtg_node_int=np.append(vtg_node_int,vtg_node_local,axis=1) # to speed up
        if np.size(vtg_node_int,1)>=700:        
            vtg_node_int2=np.append(vtg_node_int2,vtg_node_int[:,0:350],axis=1)
            vtg_node_int=vtg_node_int[:,350:700]        

        time2=time2+1   
        # print(time,con_sw_num)
    vtg_node_int2=np.append(vtg_node_int2,vtg_node_int,axis=1)        

    #toc()
    
    vtg_node_int_down=np.zeros((len(vtg_node_v),1))
    time_d=0

    for ind_d in range(0,np.size(vtg_node_int2,1)):
        if ind_d%10==0: # downsampling (to match the sampling rate of 500Hz)
            vtg_node_int_down=np.append(vtg_node_int_down,vtg_node_int2[:,ind_d].reshape(len(vtg_node_v),1),axis=1)        
            time_d=time_d+1
    
    con_sw_str=('patient folder path.../sim/erp_R0stim_G_%d_sigma_0_5' %(con_sw_num))
    np.savetxt(con_sw_str, vtg_node_int_down)
    
    con_sw=con_sw+0.001
    con_sw_num=con_sw_num+10

