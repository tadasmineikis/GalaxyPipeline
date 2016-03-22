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

#def is_number(s):
#    try:
#        float(s)
#        return True
#    except ValueError:
#        return False


    
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
    
    XTICKS=np.arange(XMIN,XMAX+0.1,0.5)
    if XTICKS.size>5:
        XTICKS=np.linspace(XMIN,XMAX,5)
    
    YTICKS=np.arange(YMIN,YMAX+0.1,0.5)
    if YTICKS.size>5:
        YTICKS=np.linspace(YMIN,YMAX,5)
        
    return {'XLIM':[XMIN,XMAX],'YLIM':[YMIN,YMAX],\
            'XLABEL':Keys[0].replace('_',''),'YLABEL':Keys[1].replace('_',''),\
           'YTICKS':YTICKS,'XTICKS':XTICKS, 'kwdict':{}}
         

    
    
def PlotCmds(Model_cmd, Model_0d, Ages, Filters, Iterations, File, Path):
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
            PLOT_PRM=SetBasicPlotParams(m0d['t']*1e-3,m0d['TSFR'], ['t','TSFR'])
            PLOT_PRM['kwdict']['alpha']=0.2
            MakePlot(ax,m0d['t']*1e-3,m0d['TSFR'],PLOT_PRM)
            tmin=np.amin(m0d['t']*1e-3)
            tmax=np.amax(m0d['t']*1e-3)
            tstep=(m0d['t'][1]-m0d['t'][0])
            bins=np.arange(tmin,tmax+0.1,0.5)
            ax.hist(m0d['t']*1e-3,weights=m0d['TSFR']/(500./tstep),bins=bins,histtype='step',lw=2,zorder=1.,log=True, color='k')
            ax.set_ylim(bottom=1e-5)
            
            ax2=ax.twinx()
            PLOT_PRM=SetBasicPlotParams(m0d['t']*1e-3,m0d['GAS'], ['t','GAS'])
            PLOT_PRM['YLOG']=True
            PLOT_PRM['YLABEL_COLOR']='red'
            PLOT_PRM['POINTS_COLOR']='red'
            MakePlot(ax2,m0d['t']*1e-3,m0d['GAS'],PLOT_PRM)
            GasMax=np.ceil( np.log10(max(m0d['GAS'])) )
            ax2.set_ylim(bottom=1e4, top=10**GasMax)
            ax2.set_yticks(np.logspace(4, GasMax, int(GasMax-4)+1))
            
            ax = plt.subplot(gs[1, 1]) #mh
            PLOT_PRM=SetBasicPlotParams(m0d['t']*1e-3,m0d['ZGAS'], ['t','ZGAS'])
            PLOT_PRM['YLOG']=True
            MakePlot(ax,m0d['t']*1e-3,m0d['ZGAS'],PLOT_PRM)
            ax.set_ylim(1e-5, 0.2)
            ax.set_yticks(np.logspace(-5, -1, 5))
            ax2=ax.twinx()
            PLOT_PRM=SetBasicPlotParams(m0d['t']*1e-3,m0d['STARS'], ['t','STARS'])
            PLOT_PRM['YLOG']=True
            PLOT_PRM['YLABEL_COLOR']='red'
            PLOT_PRM['POINTS_COLOR']='red'
            MakePlot(ax2,m0d['t']*1e-3,m0d['STARS'],PLOT_PRM)
            StarsMax=np.ceil( np.log10(max(m0d['STARS'])) )
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
    
def MainPlots(File, Model, Iterations=15):
    subp.call('mkdir '+'Dir_'+File,shell=True,executable='/bin/sh', cwd=os.getcwd())
    PATH='Dir_'+File+"/"
#    OutPrm=ModelOutputFiles(Params)
#    MODEL=ReadModelOutput(os.getcwd()+'/',File,Iterations,OutPrm)
    
    try:
        PAIRS={0:['r','mgas','YLOG'],1:['r','mstr','YLOG'],2:['r','zgas','YLOG'],3:['r','SFR','YLOG'],\
           4:['r','SFR100','YLOG'],5:['r','Ogas_tot','YLOG'],6:['r','Ogas_cur','YLOG'],\
           7:['r','Ometals_tot','YLOG'], 8:['r','Ometals_cur','YLOG'],\
           9:['r','SF_events','YLOG'],10:['r','SP_events','YLOG'],11:['r','Tgas','YLOG']}
        
        PlotGenericType(Model.OUTPUT_FILES['1d'], Model.MODEL['1d'], PAIRS, Iterations, File+'_1d_', PATH, PlotColumns=3)
    except Exception as e:
        print '1d ouput failed'
        print repr(e)
        
    try:
        PAIRS={0:['r','mgas','YLOG'],   1:['r','mstr','YLOG'], 2:['r','zgas','YLOG'],\
                   3:['r','sfe','YLOG'],4:['r','last_mstr','YLOG'],5:['r','ref_t','YLOG'],\
                   6:['r','sfr_t','YLOG'],  7:['r','TVEL']}
        
        PlotGenericType(Model.OUTPUT_FILES['2d'], Model.MODEL['2d'], PAIRS, Iterations, File+'_2d_', PATH, PlotColumns=2)
    except Exception as e:
        print '2d ouput failed'
        print repr(e)

    try:
        PAIRS={0:['t','TSFR','YLOG'],   1:['t','ACC','YLOG'], 2:['t','SP_E','YLOG'],\
                   3:['t','TR_E','YLOG'],4:['t','OTFL','YLOG'],5:['t','ST_GAS_ACC','YLOG']}
        
        PlotGenericType(Model.OUTPUT_FILES['0d'], Model.MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except Exception as e:
        print '0d ouput failed'
        print repr(e)
    Cmd=0
    try:
        Cmd=PlotCmds(Model.MODEL['cmd'], Model.MODEL['0d'], Model.OUTPUT_FILES['cmd'], ['o_B', 'o_I'], Iterations, File, PATH)
        
        PlotGenericType(Model.OUTPUT_FILES['0d'], Model.MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except Exception as e:
        print 'cmd ouput failed'
        print repr(e)
    
    for key in Model.OUTPUT_FILES.keys():
        try:
            Model.WriteVOFiles(key, PATH)
        except:
            pass
#    
#    if Cmd !=0 :
#        WriteVOFiles_CMD_Only(PATH, File, Cmd,Iterations, OutPrm)
#    
    global CLEAN_UP
    if CLEAN_UP==True:
        Model.clean_up(PATH)

def MainRun(File, Models, Iterations):
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
                template=open('template').read()
                SEED=int(random.SystemRandom().random()*1e8)
                tmp=template.replace('SEED', 'seed '+str(SEED))
                tmp=tmp.replace('GALEMO_RESULTS', 'galemo_results '+str(lines)+' '+cmd_out)
                tmp=tmp.replace('OUT', 'out '+cmd_out.replace('dat', 'cmd'))
                subp.call('echo \"'+tmp+'\" >'+'cmd_'+File+'-'+str(i)+'_'+str(akey),shell=True,executable='/bin/sh')
                subp.call('./gCMD_0.21.5 cmd_'+File+'-'+str(i)+'_'+str(akey),shell=True,executable='/bin/sh')
        except:
            print 'CMDs missing in the output'
    print "Calculations of the models complete!" 
    
def Main(File, Iterations):
#    params=ReadModelParameters(os.getcwd()+'/',File)
    MODELS=Model_IO(os.getcwd()+'/', File, Iterations)
    
#    MainRun(File, MODELS, Iterations)
    
    MODELS.ReadModelOutput()

    MainPlots(File, MODELS, Iterations)
    
if __name__=='__main__':
    Main(File=sys.argv[1], Iterations=int(sys.argv[2]))
    
    
    
