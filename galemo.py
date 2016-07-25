#! /usr/bin/env python
# Python script for 2-D galaxy disk evolution model output initial analysis

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import subprocess as subp
import os
import random
from ModelIO import Model_IO


#global variables
CMD_COADDED_OUTPUT=True
CLEAN_UP=True
    
def MakePlot(ax,X,Y,PLOT_PRM):
    if 'LINE' in PLOT_PRM.keys():
        ax.plot(X,Y,c='k',lw=1)
    else:
        if 'POINTS_COLOR' in PLOT_PRM.keys():
            ax.scatter(X,Y,edgecolor='none',facecolor=PLOT_PRM['POINTS_COLOR'],s=5, **PLOT_PRM['kwdict'])
        else:
            ax.scatter(X,Y,edgecolor='none',facecolor='k',s=5, **PLOT_PRM['kwdict'])
    for key in PLOT_PRM.keys():
        if key == 'XLIM':
            ax.set_xlim(PLOT_PRM[key])
        elif key == 'YLIM':
            if 'YINVERSE' in PLOT_PRM.keys():
                ax.set_ylim(PLOT_PRM[key][::-1])
            else:
                ax.set_ylim(PLOT_PRM[key])
        elif key == 'YLOG':
            ax.set_yscale('log')
        elif key == 'XLABEL':
            if 'CMD' in PLOT_PRM.keys():
                ax.set_xlabel(PLOT_PRM[key]+' - '+PLOT_PRM['YLABEL'], size=20)
            else:
                ax.set_xlabel(PLOT_PRM[key], size=20)
        elif key == 'YLABEL':
            if 'YLABEL_COLOR' in PLOT_PRM.keys():
                ax.set_ylabel(PLOT_PRM[key], size=20, color=PLOT_PRM['YLABEL_COLOR'])
            else:
                ax.set_ylabel(PLOT_PRM[key], size=20)
        elif key == 'YTICKS':
            ax.set_yticks(PLOT_PRM[key])
        elif key == 'XTICKS':
            ax.set_xticks(PLOT_PRM[key])
        ax.tick_params(axis='both', which='major', labelsize=15,length=10,width=1)

def SetBasicPlotParams(X,Y, Keys):
    XMIN=np.floor( np.amin(X) )
    XMAX=np.ceil( np.amax(X) )
    if 'YLOG' in Keys:
        finite=np.isfinite(Y)
        y_finite=Y[finite]
        finite=np.nonzero(y_finite)
        ylog=np.log10(y_finite[finite])
        YMIN=np.floor( np.amin(ylog) )
        YMAX=np.ceil( np.amax(ylog) )
    else:
        YMIN=np.floor( np.amin(Y) )
        YMAX=np.ceil( np.amax(Y) )
    
    if YMIN==YMAX:
        YMIN=YMIN-0.5
        YMAX=YMAX-0.5

    if np.floor(XMAX+0.1)/XMIN <= 5:
        XTICKS=np.arange(XMIN,XMAX+0.1,0.5)
    else:
        XTICKS=np.linspace(XMIN,XMAX,5)
        
    if np.floor(YMAX+0.1)/YMIN <= 5:
        YTICKS=np.arange(YMIN,YMAX+0.1,0.5)
    else:
        YTICKS=np.linspace(YMIN,YMAX,5)
        
    return {'XLIM':[XMIN,XMAX],'YLIM':[YMIN,YMAX],\
            'XLABEL':Keys[0].replace('_',''),'YLABEL':Keys[1].replace('_',''),\
           'YTICKS':YTICKS,'XTICKS':XTICKS, 'kwdict':{}}
           
def PlotCmds(Model_cmd, Model_0d, Ages, Filters, Iterations, File, Path,  t_axis='t-myr', \
                        tsfr_axis='TSFR-msol-yr',  gas_axis='GAS-msol',  stars_axis='STARS-msol',  zgas_axis='ZGAS'):
    CMD={}
    SHARED_PRM=0
    FIRST_TIME=True
    for idx in range(Iterations):
        CMD[idx]={}
    for akey in Ages:
        for idx in range(Iterations):
            fig=plt.figure(figsize=(12,8))
            gs = gridspec.GridSpec(2,2)
            
            cmd=Model_cmd[idx][akey]
            ax = plt.subplot(gs[:, 0]) #cmd
            PLOT_PRM=SetBasicPlotParams(cmd[Filters[0]]-cmd[Filters[1]],cmd[Filters[1]], Filters)
            if FIRST_TIME:
                SHARED_PRM=PLOT_PRM
                FIRST_TIME=False
            else:
                PLOT_PRM=SHARED_PRM
            PLOT_PRM['YINVERSE']=True
            PLOT_PRM['CMD']=True
            MakePlot(ax,cmd[Filters[0]]-cmd[Filters[1]],cmd[Filters[1]],PLOT_PRM)
            
            global CMD_COADDED_OUTPUT
            if CMD_COADDED_OUTPUT == False:
                del CMD[idx][akey], cmd
            
            ax = plt.subplot(gs[0, 1]) #sfh
            akey_0d=Ages[-1]
            m0d=Model_0d[idx][akey_0d]
            PLOT_PRM=SetBasicPlotParams(m0d[t_axis]*1e-3,m0d[tsfr_axis], [t_axis,tsfr_axis])
            PLOT_PRM['kwdict']['alpha']=0.2
            MakePlot(ax,m0d[t_axis]*1e-3,m0d[tsfr_axis],PLOT_PRM)
            tmin=np.amin(m0d[t_axis]*1e-3)
            tmax=np.amax(m0d[t_axis]*1e-3)
            tstep=(m0d[t_axis][1]-m0d[t_axis][0])
            bins=np.arange(tmin,tmax+0.1,0.5)
            ax.hist(m0d[t_axis]*1e-3,weights=m0d[tsfr_axis]/(500./tstep),bins=bins,histtype='step',lw=2,zorder=1.,log=True, color='k')
            ax.set_ylim(bottom=1e-5)
            
            ax2=ax.twinx()
            PLOT_PRM=SetBasicPlotParams(m0d[t_axis]*1e-3,m0d[gas_axis], [t_axis,gas_axis])
            PLOT_PRM['YLOG']=True
            PLOT_PRM['YLABEL_COLOR']='red'
            PLOT_PRM['POINTS_COLOR']='red'
            MakePlot(ax2,m0d[t_axis]*1e-3,m0d[gas_axis],PLOT_PRM)
            GasMax=np.ceil( np.log10(max(m0d[gas_axis])) )
            ax2.set_ylim(bottom=1e4, top=10**GasMax)
            ax2.set_yticks(np.logspace(4, GasMax, int(GasMax-4)+1))
            
            ax = plt.subplot(gs[1, 1]) #mh
            PLOT_PRM=SetBasicPlotParams(m0d[t_axis]*1e-3,m0d[zgas_axis], [t_axis,zgas_axis])
            PLOT_PRM['YLOG']=True
            MakePlot(ax,m0d[t_axis]*1e-3,m0d[zgas_axis],PLOT_PRM)
            ax.set_ylim(1e-5, 0.2)
            ax.set_yticks(np.logspace(-5, -1, 5))
            ax2=ax.twinx()
            PLOT_PRM=SetBasicPlotParams(m0d[t_axis]*1e-3,m0d[stars_axis], [t_axis,stars_axis])
            PLOT_PRM['YLOG']=True
            PLOT_PRM['YLABEL_COLOR']='red'
            PLOT_PRM['POINTS_COLOR']='red'
            MakePlot(ax2,m0d[t_axis]*1e-3,m0d[stars_axis],PLOT_PRM)
            StarsMax=np.ceil( np.log10(max(m0d[stars_axis])) )
            ax2.set_ylim(bottom=100, top=10**StarsMax)
            ax2.set_yticks(np.logspace(2, StarsMax, int(StarsMax-2)+1))
            
            fig.savefig(Path+File+'_CMD_'+str(akey)+'_ITERATION-'+str(idx)+'.png', dpi=150)
            plt.close()
            
    if CMD_COADDED_OUTPUT == True:
        return CMD
    else:
        return 0

def PlotGenericType(OutPrm, Model, Pairs, Iterations, File, Path, PlotColumns=3):
    for akey in OutPrm:
        fig, ax = plt.subplots(len(Pairs.keys())/PlotColumns,PlotColumns,figsize=(16,16))
        fig.subplots_adjust(hspace=0.25,wspace=0.3)
        for idx in range(Iterations):
            for plkey in Pairs.keys():
                PLOT_PRM=SetBasicPlotParams(Model[idx][akey][Pairs[plkey][0]],Model[idx][akey][Pairs[plkey][1]], Pairs[plkey])
                for idx in range(Iterations):
                    if 'YLOG' in Pairs[plkey]:
                        finite=np.isfinite(Model[idx][akey][Pairs[plkey][1]])
                        x_finite=Model[idx][akey][Pairs[plkey][0]][finite]
                        y_finite=Model[idx][akey][Pairs[plkey][1]][finite]
                        finite=np.nonzero( y_finite )
                        
                        finite=np.nonzero( Model[idx][akey][Pairs[plkey][1]] )
                        MakePlot(ax.flat[plkey],x_finite[finite], np.log10(y_finite[finite]),PLOT_PRM)
                    else:
                        finite=np.isfinite(Model[idx][akey][Pairs[plkey][1]])
                        MakePlot(ax.flat[plkey],Model[idx][akey][Pairs[plkey][0]][finite],\
                                        Model[idx][akey][Pairs[plkey][1]][finite],PLOT_PRM)
        fig.savefig(Path+File+'AGE_'+str(akey)+'.png', dpi=150)
        plt.close()
    
def MainPlots(File, Model, Iterations=15, CmdPhotSys='UBV'):
    subp.call('mkdir '+'Dir_'+File,shell=True,executable='/bin/sh', cwd=os.getcwd())
    PATH='Dir_'+File+"/"
#    OutPrm=ModelOutputFiles(Params)
#    MODEL=ReadModelOutput(os.getcwd()+'/',File,Iterations,OutPrm)
    
    try:
        PAIRS={0:['r-kpc','mgas-msol-pc2','YLOG'],1:['r-kpc','mstr-msol-pc2','YLOG'],2:['r-kpc','zgas','YLOG'],3:['r-kpc','SFR-msol-pc2-tstep','YLOG'],\
           4:['r-kpc','SFR100-msol-pc2-10tsteps','YLOG'],5:['r-kpc','Ogas_tot-msol','YLOG'],6:['r-kpc','Ogas_cur-msol-pc2-10tsteps','YLOG'],\
           7:['r-kpc','Ometals_tot-msol','YLOG'], 8:['r-kpc','Ometals_cur-msol-pc2-10tsteps','YLOG'],\
           9:['r-kpc','SF_events-num','YLOG'],10:['r-kpc','SP_events-num','YLOG'],11:['r-kpc','Tgas-msol-pc2','YLOG']}
        
        PlotGenericType(Model.OUTPUT_FILES['1d'], Model.MODEL['1d'], PAIRS, Iterations, File+'_1d_', PATH, PlotColumns=3)
    except Exception as e:
        print '1d ouput failed'
        print repr(e)
        
    try:
        PAIRS={0:['r-kpc','mgas-msol-pc2','YLOG'],   1:['r-kpc','mstr-msol-pc2','YLOG'], 2:['r-kpc','zgas','YLOG'],\
                   3:['r-kpc','sfe-ratio','YLOG'],4:['r-kpc','last_mstr-msol','YLOG'],5:['r-kpc','ref_t-tstep','YLOG'],\
                   6:['r-kpc','sfr_t','YLOG'],  7:['r-kpc','TVEL']}
        
        PlotGenericType(Model.OUTPUT_FILES['2d'], Model.MODEL['2d'], PAIRS, Iterations, File+'_2d_', PATH, PlotColumns=2)
    except Exception as e:
        print '2d ouput failed'
        print repr(e)

    try:
        PAIRS={0:['t-myr','TSFR-msol-yr','YLOG'],   1:['t-myr','ACC-msol-yr','YLOG'], 2:['t-myr','SP_E-num','YLOG'],\
                   3:['t-myr','TR_E-num','YLOG'],4:['t-myr','OTFL-msol-yr','YLOG'],5:['t-myr','ST_GAS_ACC-msol-yr','YLOG']}
        
        PlotGenericType(Model.OUTPUT_FILES['0d'], Model.MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except Exception as e:
        print '0d ouput failed'
        print repr(e)

    try:
        if CmdPhotSys == 'UBV':
            PlotCmds(Model.MODEL['cmd'], Model.MODEL['0d'], Model.OUTPUT_FILES['cmd'], ['o_B', 'o_I'], Iterations, File, PATH)
        elif CmdPhotSys == 'WFC-ACS':
            PlotCmds(Model.MODEL['cmd'], Model.MODEL['0d'], Model.OUTPUT_FILES['cmd'], ['o_F606W', 'o_F814W'], Iterations, File, PATH)
        PlotGenericType(Model.OUTPUT_FILES['0d'], Model.MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except Exception as e:
        print 'cmd ouput failed'
        print repr(e)
    
    SUCCESS=False
    for key in Model.OUTPUT_FILES.keys():
        try:
            SUCCESS=Model.WriteCsvFiles(key, PATH, Compression='gzip')
        except Exception as e:
            print 'Coaded files ouput for ' +str(key)+' failed'
            print repr(e)
    #make clean up only if coaded output operation was successfull
    global CLEAN_UP
    if CLEAN_UP==True and SUCCESS:
        Model.clean_up(PATH)

    

def MainRun(File, Models, Iterations, CmdPhotSys='UBV'):
    for i in range(Iterations):
        subp.call('cp '+File+' '+File+'-'+str(i),shell=True,executable='/bin/sh')
        SEED=int(random.SystemRandom().random()*1e8)
        subp.call('sed -i \"s/SEED.*/SEED '+str(SEED)+'/\" '+File+'-'+str(i)+'',shell=True,executable='/bin/sh')
        subp.call('./galaxy_2.0 '+File+'-'+str(i),shell=True,executable='/bin/sh')
        
        try:
            for akey in Models.OUTPUT_FILES['cmd']:
                cmd_out=File+'-'+str(i)+'_cmd_'+str(akey)+'.dat'
                
                lines=subp.Popen('cat '+cmd_out +' | wc -l', executable='/bin/sh', shell=True, stdout=subp.PIPE)
                lines=lines.communicate()[0].rstrip('\n')
                lines=int(lines)-1 # to account for the header
                template=-1
                if CmdPhotSys=='UBV':
                    template=open('template-ubv').read()
                elif CmdPhotSys=='WFC-ACS':
                    template=open('template-wfc_acs').read()
                SEED=int(random.SystemRandom().random()*1e8)
                tmp=template.replace('SEED', 'seed '+str(SEED))
                template.close()
                tmp=tmp.replace('GALEMO_RESULTS', 'galemo_results '+str(lines)+' '+cmd_out)
                tmp=tmp.replace('OUT', 'out '+cmd_out.replace('dat', 'cmd'))
                subp.call('echo \"'+tmp+'\" >'+'cmd_'+File+'-'+str(i)+'_'+str(akey),shell=True,executable='/bin/sh')
                
                if CmdPhotSys=='UBV':
                    subp.call('./gCMD_0.21.5_ubv cmd_'+File+'-'+str(i)+'_'+str(akey),shell=True,executable='/bin/sh')
                elif CmdPhotSys=='WFC-ACS':
                    subp.call('./gCMD_0.21.5_acs cmd_'+File+'-'+str(i)+'_'+str(akey),shell=True,executable='/bin/sh')
                
        except:
            print 'CMDs missing in the output'
    print "Calculations of the models complete!" 
    
def Main(File, Iterations, CmdPhotSys):
    MODELS=Model_IO(os.getcwd()+'/', File, Iterations)
    
    MainRun(File, MODELS, Iterations, CmdPhotSys)
    
    MODELS.ReadModelOutput()

    MainPlots(File, MODELS, Iterations, CmdPhotSys)
    
if __name__=='__main__':
    CMD_PHOT_SYS=''
    try:
        if sys.argv[3]=='cmd-acs':
            CMD_PHOT_SYS='WFC-ACS'
        elif sys.argv[3]=='cmd-ubv':
            CMD_PHOT_SYS='UBV'
        else:
            print 'Passed 3rd parameter:', sys.argv[3], ' doesn\'t look like a defintion of phot sys for the cmd [cmd-acs or cmd-ubv]'
    except:
        print 'UBVRIJHK phot system will be used or cmd\'s '
        CMD_PHOT_SYS='UBV'
    Main(File=sys.argv[1], Iterations=int(sys.argv[2]), CmdPhotSys=CMD_PHOT_SYS)
    
    
    
