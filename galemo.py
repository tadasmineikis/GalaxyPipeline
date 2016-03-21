#! /usr/bin/env python
# Python script for 2-D galaxy disk evolution model output initial analysis

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import subprocess as subp
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
import os
import random

#global variables
CMD_COADDED_OUTPUT=True
CLEAN_UP=True

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def clean_up(Path, File, OutPrm, Prm, Iterations):
    ZIP_FILE_LIST=[]
    RM_FILE_LIST=[]
    
    DICT={}
    for key in OutPrm.keys():
        if key == '0d':
            DICT['_igal.dat']=True
        elif key == '1d':
            for akey in OutPrm[key]:
                DICT['_rings_'+str(akey)+'.dat']=True
        elif key == '2d':
            for akey in OutPrm[key]:
                DICT['_cells_'+str(akey)+'.dat']=True
        elif key == 'cmd':
            for akey in OutPrm[key]:
                DICT['_cmd_'+str(akey)+'.dat']=True
                DICT['_cmd_'+str(akey)+'.cmd']=True
    DICT['']=True
    DICT['.log']=True
    for itr in range(Iterations):
        for key in DICT.keys():
            ZIP_FILE_LIST.append(File+'-'+str(itr)+key)
            RM_FILE_LIST.append(File+'-'+str(itr)+key)
    for itr in range(Iterations):
        for key in ['cmd_']:
            ZIP_FILE_LIST.append(key+File+'-'+str(itr))
            RM_FILE_LIST.append(key+File+'-'+str(itr))

    ZIP_FILE_LIST.append(Prm['ACCRETION'][0])
    ZIP_FILE_LIST.append(Prm['ACCRETION'][1])
    
    if is_number(Prm['SFE'][0]):
        pass
    else:
        ZIP_FILE_LIST.append(Prm['SFE'][0])
    
    if is_number(Prm['SFE_POW'][0]):
        pass
    else:
        ZIP_FILE_LIST.append(Prm['SFE_POW'][0])
        
    if is_number(Prm['TRIGGERED'][0]):
        pass
    else:
        ZIP_FILE_LIST.append(Prm['TRIGGERED'][0])
    
    ZIP_FLIST=''
    RM_FLIST=''
    for line in ZIP_FILE_LIST:
        ZIP_FLIST+=line+' '
    
    for line in RM_FILE_LIST:
        RM_FLIST+=line+' '
    
    subp.call('zip '+Path+File+'.zip '+ZIP_FLIST,shell=True,executable='/bin/sh',cwd=os.getcwd())
    subp.call('mv '+File+' '+Path+File,shell=True,executable='/bin/sh',cwd=os.getcwd())
    subp.call('rm '+RM_FLIST,shell=True,executable='/bin/sh',cwd=os.getcwd())

def ReadModelParameters(Path,File):
    PARAMETERS={}
    for line in open(Path+File,'r'):
        if line[0]!='#':
            sline=line.split()
            PARAMETERS[sline[0]]=[]
            for i in range(len(sline)-1):
                PARAMETERS[sline[0]].append(sline[i+1])
    return PARAMETERS

def ModelOutputFiles(Parameters):
    OUTPUT_FILES={}
    numOfTypes=int(Parameters['OUTPUT'][0])
    for otype in Parameters['OUTPUT'][1:numOfTypes+1]:
        OUTPUT_FILES[otype]=[]
        if otype=='0d':
            OUTPUT_FILES[otype].append(int( Parameters['GALAXY_AGE'][0] ))
        else:
            # new type output times specification
            if Parameters['OUTPUT'][numOfTypes+1] =='var':
                t_time0=int(Parameters['OUTPUT'][numOfTypes+2])
                t_time1=int(Parameters['OUTPUT'][numOfTypes+3])
                t_step=int(Parameters['OUTPUT'][numOfTypes+4])
                for ind in range( (t_time1-t_time0)/t_step+1):
                    OUTPUT_FILES[otype].append(int(t_time0+t_step*ind))
            #old type output times specification
            elif is_number(Parameters['OUTPUT'][numOfTypes+1]):
                for t_time in Parameters['OUTPUT'][numOfTypes+2:]:
                    OUTPUT_FILES[otype].append(int(t_time))
    return OUTPUT_FILES
    
def MakePlot(ax,X,Y,PLOT_PRM):
    if 'LINE' in PLOT_PRM.keys():
        ax.plot(X,Y,c='k',lw=1)
    else:
        if 'POINTS_COLOR' in PLOT_PRM.keys():
            ax.scatter(X,Y,edgecolor='none',facecolor=PLOT_PRM['POINTS_COLOR'],s=5)
        else:
            ax.scatter(X,Y,edgecolor='none',facecolor='k',s=5)
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
           'YTICKS':YTICKS,'XTICKS':XTICKS}
         
def ReadModelOutput(Path,Pfile,Iterations,Parameters):
    MODEL={}
    for key in Parameters:
        if key == '2d':
            MODEL['2d']={}
            for idx in range(Iterations):
                MODEL['2d'][idx]={}
                for akey in Parameters['2d']:
                    MODEL['2d'][idx][akey]=np.genfromtxt(Path+Pfile+'-'+str(idx)+'_cells_'+str(akey)+'.dat',\
                                   names=True,dtype='float32')
        elif key == '1d':
            MODEL['1d']={}
            for idx in range(Iterations):
                MODEL['1d'][idx]={}
                for akey in Parameters['1d']:
                    MODEL['1d'][idx][akey]=np.genfromtxt(Path+Pfile+'-'+str(idx)+'_rings_'+str(akey)+'.dat',\
                                   names=True,dtype='float32')
        elif key == '0d':
            MODEL['0d']={}
            for idx in range(Iterations):
                MODEL['0d'][idx]={}
                for akey in Parameters['0d']:
                    MODEL['0d'][idx][akey]=np.genfromtxt(Path+Pfile+'-'+str(idx)+'_igal.dat',\
                                   names=True,dtype='float32')
        elif key == 'cmd':
            MODEL['cmd']={}
            for idx in range(Iterations):
                MODEL['cmd'][idx]={}
                for akey in Parameters['cmd']:
                    MODEL['cmd'][idx][akey]=Path+Pfile+'-'+str(idx)+'_cmd_'+str(akey)+'.cmd'
    return MODEL
    
    
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
            
            CMD[idx][akey]=np.genfromtxt(Model_cmd[idx][akey],names=True,dtype='f4')
            cmd=CMD[idx][akey]
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
            MakePlot(ax,m0d['t']*1e-3,m0d['TSFR'],PLOT_PRM)
#            ax.plot(m0d['t']*1e-3,m0d['ACC'], color='b')
            tmin=np.amin(m0d['t']*1e-3)
            tmax=np.amax(m0d['t']*1e-3)
            tstep=(m0d['t'][1]-m0d['t'][0])
            bins=np.arange(tmin,tmax+0.1,0.5)
            ax.hist(m0d['t']*1e-3,weights=m0d['TSFR']/(500./tstep),bins=bins,histtype='step',lw=2,zorder=1.,log=True, color='m')
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
#            ax.plot(m0d['t']*1e-3,m0d['Zgas'])
#            tmin=np.amin(m0d['t']*1e-3)
#            tmax=np.amax(m0d['t']*1e-3)
#            tstep=(m0d['t'][1]-m0d['t'][0])
#            bins=np.arange(tmin,tmax+0.1,0.5)
#            ax.hist(m0d['t']*1e-3,weights=m0d['TSFR']/(500./tstep),bins=bins,histtype='step',lw=2,zorder=1.,log=True)
#            ax.set_ylim(bottom=1e-5)
            
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

#def WriteVOFiles(Path, File, Model, Iterations, OutPrm):
#    for key in OutPrm.keys():
#        if key =='cmd':
#            pass
#        else:
#            Res={}
#            # Create a new VOTable file...
#            votable = VOTableFile()
#            
#            Res[key]={}
#            Res[key]['resource']=Resource(name=key)
#            votable.resources.append(Res[key]['resource'])
#            for midx in range(Iterations):
#                if key =='0d':
#                    Res[key]['table'] = Table(votable,name='Table_'+key+'-'+str(midx),\
#                    nrows=Model[key][midx][OutPrm[key][0]].size)
#                else:
#                    Res[key]['table'] = Table(votable,name='Table_'+key+'-'+str(midx),\
#                    nrows=Model[key][midx][OutPrm[key][0]].size)
#                Res[key]['resource'].tables.append(Res[key]['table'])
#                FIELDS=[]
#                for akey in OutPrm[key]:
#                    for ckey in Model[key][midx][akey].dtype.names:
#                        FIELDS.append( Field(votable, name=ckey,\
#                                         datatype="float", arraysize="1") )
#                   
#                Res[key]['table'].fields.extend(FIELDS)
#                Res[key]['table'].create_arrays(Model[key][midx][OutPrm[key][0]].size)
#                Names=Model[key][midx][OutPrm[key][0]].dtype.names
#                for akey in OutPrm[key]:
#                    Res[key]['table'].array[:]=Model[key][midx][akey][list(Names)][:]
#            votable.to_xml(Path+File+'_'+str(key)+ ".xml", tabledata_format="binary")
#            del votable, Res, FIELDS
def WriteVOFiles(Path, File, Model, Iterations, OutPrm):
    for key in OutPrm.keys():
        if key =='cmd':
            pass
        else:
            Res={}
            # Create a new VOTable file...
            votable = VOTableFile()
            
            Res[key]={}
            Res[key]['resource']=Resource(name=key)
            votable.resources.append(Res[key]['resource'])
            
            a_age=[]
            a_itr=[]
            a_model=[]
            for itr in Model[key].keys():
                for akey in Model[key][itr].keys():
                    a_model.append(Model[key][itr][akey])
                    a_age.append(np.repeat(akey,Model[key][itr][akey].size))
                    a_itr.append(np.repeat(itr,Model[key][itr][akey].size))
            st_model=np.hstack(a_model)
            st_age=np.hstack(a_age)
            st_itr=np.hstack(a_itr)

            a_stack=np.hstack( [st_model.view(np.float32).reshape(st_model.shape+(-1,)),\
                                st_age.reshape(st_age.shape+(-1,)),\
                              st_itr.reshape(st_itr.shape+(-1,))])
            
            FIELDS=[]
            for midx in range(Iterations):
                for akey in OutPrm[key]:
                    for ckey in Model[key][midx][akey].dtype.names:
                        FIELDS.append( Field(votable, name=ckey,\
                                         datatype="float", arraysize="1") )
                    size=Model[key][midx][akey].size*Iterations*len(OutPrm[key])
                    names=Model[key][midx][akey].dtype.names
                    dt=Model[key][midx][akey].dtype
                    break
                break

            names+=('age','iteration',)
            dt=[]
            for name in names:
                if name=='age':
                    dt.append((name,'i4'))
                elif name=='iteration' or name =='age':
                    dt.append((name,'i4'))
                else:
                    dt.append((name,'f4'))
            
            stack=np.empty(a_stack.shape[0],dtype=dt)
            
            for i, nm in zip(range(a_stack.shape[1]),stack.dtype.names):
                stack[nm]=a_stack[:,i]
    
            FIELDS.append( Field(votable, name='age',\
                                         datatype="int", arraysize="1") )
            FIELDS.append( Field(votable, name='iteration',\
                                         datatype="int", arraysize="1") )
            if key =='0d':
                Res[key]['table'] = Table(votable,name='Table_'+key,\
                    nrows=size)
            else:
                Res[key]['table'] = Table(votable,name='Table_'+key,\
                    nrows=size)
            Res[key]['resource'].tables.append(Res[key]['table'])
                   
            Res[key]['table'].fields.extend(FIELDS)
            Res[key]['table'].create_arrays(size)

            Res[key]['table'].array[:]=stack[list(stack.dtype.names)][:]
            
            votable.to_xml(Path+File+'_'+str(key)+ ".xml", tabledata_format="binary", compressed=True)
            del votable, Res, FIELDS
            
def WriteVOFiles_CMD_Only(Path, File, Cmd,Iterations, OutPrm):
    itr=[]
    ages=[]
    cmd=[]
    for iterat in Cmd.keys():
        for akey in Cmd[iterat].keys():
            cmd.append(Cmd[iterat][akey])
            itr.append(np.repeat(iterat,Cmd[iterat][akey].size))
            ages.append(np.repeat(akey,Cmd[iterat][akey].size))
        
    a_cmd=np.hstack(cmd)
    a_iter=np.hstack(itr)
    a_ages=np.hstack(ages)
    
    dt=[]
    for name in a_cmd.dtype.names:
        dt.append((name,'f4'))
    dt.append(('ModelAge','i4'))
    dt.append(('iteration','i4'))
    
    cmd_stack=np.empty(a_cmd.size,dtype=dt)
    
    for name in a_cmd.dtype.names:
        cmd_stack[name]=a_cmd[name]
    cmd_stack['ModelAge']=a_ages
    cmd_stack['iteration']=a_iter
    
    Res={}
    # Create a new VOTable file...
    votable = VOTableFile()

    Res={}
    Res['resource']=Resource(name='cmd')
    votable.resources.append(Res['resource'])

    FIELDS=[]
    for ckey in a_cmd.dtype.names:
        if ckey=='ssp_index':
            FIELDS.append( Field(votable, name=ckey,datatype="int", arraysize="1") )
        else:
            FIELDS.append( Field(votable, name=ckey,datatype="float", arraysize="1") )
    size=a_cmd.size

    FIELDS.append( Field(votable, name='ModelAge', datatype="int", arraysize="1") )
    FIELDS.append( Field(votable, name='iteration',datatype="int", arraysize="1") )

    Res['table'] = Table(votable,name='Table_cmd',  nrows=size)
    Res['resource'].tables.append(Res['table'])

    Res['table'].fields.extend(FIELDS)
    Res['table'].create_arrays(size)
    

    Res['table'].array[:]=cmd_stack[list(cmd_stack.dtype.names)][:]

    votable.to_xml(Path+File+"_cmd.xml", tabledata_format="binary", compressed=True)
    del votable, Res, FIELDS, Cmd, a_cmd, cmd_stack
    
def MainPlots(File, Params, Iterations=15):
    subp.call('mkdir '+'Dir_'+File,shell=True,executable='/bin/sh', cwd=os.getcwd())
    PATH='Dir_'+File+"/"
    OutPrm=ModelOutputFiles(Params)
    MODEL=ReadModelOutput(os.getcwd()+'/',File,Iterations,OutPrm)
    
    try:
        PAIRS={0:['r','mgas','YLOG'],1:['r','mstr','YLOG'],2:['r','zgas','YLOG'],3:['r','SFR','YLOG'],\
           4:['r','SFR100','YLOG'],5:['r','Ogas_tot','YLOG'],6:['r','Ogas_cur','YLOG'],\
           7:['r','Ometals_tot','YLOG'], 8:['r','Ometals_cur','YLOG'],\
           9:['r','SF_events','YLOG'],10:['r','SP_events','YLOG'],11:['r','Tgas','YLOG']}
        
        PlotGenericType(OutPrm['1d'], MODEL['1d'], PAIRS, Iterations, File+'_1d_', PATH, PlotColumns=3)
    except:
        print '1d missing in the ouput'
        
    try:
        PAIRS={0:['r','mgas','YLOG'],   1:['r','mstr','YLOG'], 2:['r','zgas','YLOG'],\
                   3:['r','sfe','YLOG'],4:['r','last_mstr','YLOG'],5:['r','ref_t','YLOG'],\
                   6:['r','sfr_t','YLOG'],  7:['r','TVEL']}
        
        PlotGenericType(OutPrm['2d'], MODEL['2d'], PAIRS, Iterations, File+'_2d_', PATH, PlotColumns=2)
    except:
        print '2d missing in the ouput'

    try:
        PAIRS={0:['t','TSFR','YLOG'],   1:['t','ACC','YLOG'], 2:['t','SP_E','YLOG'],\
                   3:['t','TR_E','YLOG'],4:['t','OTFL','YLOG'],5:['t','ST_GAS_ACC','YLOG']}
        
        PlotGenericType(OutPrm['0d'], MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except:
        print '0d missing in the ouput'
    Cmd=0
    try:
        Cmd=PlotCmds(MODEL['cmd'], MODEL['0d'], OutPrm['cmd'], ['o_B', 'o_I'], Iterations, File, PATH)
        
        PlotGenericType(OutPrm['0d'], MODEL['0d'], PAIRS, Iterations, File+'_0d_', PATH, PlotColumns=2)
    except:
        print 'cmd missing in the ouput'
    
    WriteVOFiles(PATH, File, MODEL, Iterations, OutPrm)
    
    if Cmd !=0 :
        WriteVOFiles_CMD_Only(PATH, File, Cmd,Iterations, OutPrm)
    
    global CLEAN_UP
    if CLEAN_UP==True:
        clean_up(PATH, File, OutPrm, Params, Iterations)

def MainRun(File, Params, Iterations):
    for i in range(Iterations):
        subp.call('cp '+File+' '+File+'-'+str(i),shell=True,executable='/bin/sh')
        SEED=int(random.SystemRandom().random()*1e8)
        subp.call('sed -i \"s/SEED.*/SEED '+str(SEED)+'/\" '+File+'-'+str(i)+'',shell=True,executable='/bin/sh')
        subp.call('./galaxy_2.0 '+File+'-'+str(i),shell=True,executable='/bin/sh')
        
        OutPrm=ModelOutputFiles(Params)
        
        try:
            for akey in OutPrm['cmd']:
                cmd_out=File+'-'+str(i)+'_cmd_'+str(akey)+'.dat'
                
                lines=subp.Popen('cat '+cmd_out +' | wc -l', executable='/bin/sh', shell=True, stdout=subp.PIPE)
                lines=lines.communicate()[0].rstrip('\n')
                lines=int(lines)-1 # to account for the header
                template=open('template').read()
                SEED=int(random.SystemRandom().random()*1e8)
                tmp=template.replace('SEED', 'seed '+str(SEED))
                tmp=tmp.replace('GALEMO_RESULTS', 'galemo_results '+str(lines)+' '+cmd_out)
                tmp=tmp.replace('OUT', 'out '+cmd_out.replace('dat', 'cmd'))
                subp.call('echo \"'+tmp+'\" >'+'cmd_'+File+'-'+str(i),shell=True,executable='/bin/sh')
                subp.call('./gCMD_0.21.5 cmd_'+File+'-'+str(i),shell=True,executable='/bin/sh')
        except:
            print 'CMDs missing in the output'
    print "Calculations of the models complete!" 
    
def Main(File, Iterations):
    params=ReadModelParameters(os.getcwd()+'/',File)
    MainRun(File, params, Iterations)
    MainPlots(File, params, Iterations)
    
if __name__=='__main__':
    Main(File=sys.argv[1], Iterations=int(sys.argv[2]))
    
    
    
