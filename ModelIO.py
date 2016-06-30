# fileIO class
import numpy as np
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
import pandas as pd
import subprocess as subp
import os

class Model_IO:
    def __init__(self, Path, Pfile, Iterations):
        self.FULL_PATH=Path+Pfile
        self.FILE=Pfile
        self.ITERATIONS=Iterations
        
        
        self.ReadModelParameters()
        self.ModelOutputFiles()
        self.FileEndings()
#        self.ReadModelOutput()
        
    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False
            
    def FileEndings(self):
        self.FILE_NAMES={}
        for key in self.OUTPUT_FILES.keys():
            if key == '0d':
                self.FILE_NAMES[key]={}
                for akey in self.OUTPUT_FILES[key]:
                    self.FILE_NAMES[key][akey]='_igal.dat'
            elif key == '1d':
                self.FILE_NAMES[key]={}
                for akey in self.OUTPUT_FILES[key]:
                    self.FILE_NAMES[key][akey]='_rings_'+str(akey)+'.dat'
            elif key == '2d':
                self.FILE_NAMES[key]={}
                for akey in self.OUTPUT_FILES[key]:
                    self.FILE_NAMES[key][akey]='_cells_'+str(akey)+'.dat'
            elif key == 'cmd':
                self.FILE_NAMES[key]={}
                self.FILE_NAMES[key+'_raw']={}
                for akey in self.OUTPUT_FILES[key]:
                    self.FILE_NAMES[key][akey]='_cmd_'+str(akey)+'.cmd'
                    self.FILE_NAMES[key+'_raw'][akey]='_cmd_'+str(akey)+'.dat'
            
    def ReadModelParameters(self):
        self.PARAMETERS={}
        for line in open(self.FULL_PATH,'r'):
            if line[0]!='#':
                sline=line.split()
                self.PARAMETERS[sline[0]]=[]
                for i in range(len(sline)-1):
                    self.PARAMETERS[sline[0]].append(sline[i+1])

    def ModelOutputFiles(self):
        self.OUTPUT_FILES={}
        numOfTypes=int(self.PARAMETERS['OUTPUT'][0])
        for otype in self.PARAMETERS['OUTPUT'][1:numOfTypes+1]:
            self.OUTPUT_FILES[otype]=[]
            if otype=='0d':
                self.OUTPUT_FILES[otype].append(int( self.PARAMETERS['GALAXY_AGE'][0] ))
            else:
                # new type output times specification
                if self.PARAMETERS['OUTPUT'][numOfTypes+1] =='var':
                    t_time0=int(self.PARAMETERS['OUTPUT'][numOfTypes+2])
                    t_time1=int(self.PARAMETERS['OUTPUT'][numOfTypes+3])
                    t_step=int(self.PARAMETERS['OUTPUT'][numOfTypes+4])
                    for ind in range( (t_time1-t_time0)/t_step+1):
                        self.OUTPUT_FILES[otype].append(int(t_time0+t_step*ind))
                #old type output times specification
                elif self.is_number(self.PARAMETERS['OUTPUT'][numOfTypes+1]):
                    for t_time in self.PARAMETERS['OUTPUT'][numOfTypes+2:]:
                        self.OUTPUT_FILES[otype].append(int(t_time))

    def read_file(self, name):
        f=open(name,'r')
        head=f.readline().split()
        f.close()
        head[0]=head[0][1:]
        DTYPE={}
        for item in head:
            DTYPE[item]=np.float32
        df=pd.read_table(name,na_values='\"\"',names=head,header=0,\
        delim_whitespace=True,index_col=False,dtype=DTYPE,as_recarray=True)
        return df
    
    def ReadModelOutput(self):
        self.MODEL={}
        for key in self.OUTPUT_FILES.keys():
            self.MODEL[key]={}
            for idx in range(self.ITERATIONS):
                self.MODEL[key][idx]={}
                for akey in self.OUTPUT_FILES[key]:
                    self.MODEL[key][idx][akey]=self.read_file(self.FULL_PATH+'-'+str(idx)+ self.FILE_NAMES[key][akey])
                    
    def CoadFilesIntoRecArray(self, key):
        a_age=[]
        a_itr=[]
        a_model=[]
        for itr in self.MODEL[key].keys():
            for akey in self.MODEL[key][itr].keys():
                a_model.append(self.MODEL[key][itr][akey])
                a_age.append(np.repeat(akey,self.MODEL[key][itr][akey].size))
                a_itr.append(np.repeat(itr,self.MODEL[key][itr][akey].size))
        st_model=np.hstack(a_model)
        st_age=np.hstack(a_age)
        st_itr=np.hstack(a_itr)

        a_stack=np.hstack( [st_model.view(np.float32).reshape(st_model.shape+(-1,)),\
                            st_age.reshape(st_age.shape+(-1,)),\
                          st_itr.reshape(st_itr.shape+(-1,))])
        
        names=st_model.dtype.names+('ModelAge','iteration',)
        dt=[]
        for name in names:
            if name=='iteration' or name =='ModelAge' or name=='ssp_index':
                dt.append((name,'i4'))
            else:
                dt.append((name,'f4'))

        stack=np.empty(a_stack.shape[0],dtype=dt)
        for i, nm in zip(range(a_stack.shape[1]),stack.dtype.names):
            stack[nm]=a_stack[:,i]
        return stack
        
    def WriteCsvFiles(self, key, Output_Path, Compression=None):
        RecArray=self. CoadFilesIntoRecArray(key)
        df = pd.DataFrame(RecArray, columns=RecArray.dtype.names)
        if Compression == 'gzip':
            import gzip
            f=gzip.GzipFile(Output_Path+self.FILE+'_'+key+'.csv.gz', 'wb')
            df.to_csv(f,index=False,float_format='%.3f')
            f.close()
            return True
        elif Compression==None:
            df.to_csv(Output_Path+self.FILE+'_'+key+'.csv',index=False,float_format='%.3f')
            return True
        else:
            print "Compresion parameter", Compression, "not supported."
            return False
        
    def WriteVOFiles(self, key, Output_Path):
        Res={}
        # Create a new VOTable file...
        votable = VOTableFile()
        Res[key]={}
        Res[key]['resource']=Resource(name=key)
        votable.resources.append(Res[key]['resource'])
        
        RecArray=self. CoadFilesIntoRecArray(key)
        
        size=RecArray.size
        FIELDS=[]
        for ckey in RecArray.dtype.names:
            if ckey=='ssp_index' or ckey == 'iteration' or ckey =='ModelAge':
                FIELDS.append( Field(votable, name=ckey,datatype="int", arraysize="1") )
            else:
                FIELDS.append( Field(votable, name=ckey,datatype="float", arraysize="1") )
        Res[key]['table'] = Table(votable,name='Table_'+key,\
                nrows=size)
        Res[key]['resource'].tables.append(Res[key]['table'])
        Res[key]['table'].fields.extend(FIELDS)
        Res[key]['table'].create_arrays(size)
        Res[key]['table'].array[:]=RecArray[list(RecArray.dtype.names)][:]
        votable.to_xml(Output_Path+self.FILE+'_'+str(key)+ ".xml", tabledata_format="binary", compressed=True)
        del votable, Res, FIELDS

    def clean_up(self, Output_Path):
        ZIP_FILE_LIST=[]
        RM_FILE_LIST=[]

        for itr in range(self.ITERATIONS):
            ZIP_FILE_LIST.append(self.FILE+'-'+str(itr)+'.log')
            RM_FILE_LIST.append(self.FILE+'-'+str(itr)+'.log')
            ZIP_FILE_LIST.append(self.FILE+'-'+str(itr))
            RM_FILE_LIST.append(self.FILE+'-'+str(itr))
            for key in self.FILE_NAMES.keys():
                for akey in self.FILE_NAMES[key].keys():
                    ZIP_FILE_LIST.append(self.FILE+'-'+str(itr)+self.FILE_NAMES[key][akey])
                    RM_FILE_LIST.append(self.FILE+'-'+str(itr)+self.FILE_NAMES[key][akey])
            for key in ['cmd_']:
                try:
                    for akey in self.OUTPUT_FILES['cmd']:
                        ZIP_FILE_LIST.append(key+self.FILE+'-'+str(itr)+'_'+str(akey))
                        RM_FILE_LIST.append(key+self.FILE+'-'+str(itr)+'_'+str(akey))
                except:
                    pass
                
        ZIP_FILE_LIST.append(self.PARAMETERS['ACCRETION'][0])
        ZIP_FILE_LIST.append(self.PARAMETERS['ACCRETION'][1])
        
        for key in ['SFE', 'SFE_POW', 'TRIGGERED']:
            if not self.is_number(self.PARAMETERS[key][0]):
                ZIP_FILE_LIST.append(self.PARAMETERS[key][0])
        
        ZIP_FLIST=self.FILE+' '
        for line in ZIP_FILE_LIST:
            ZIP_FLIST+=line+' '
        RM_FLIST=''
        for line in RM_FILE_LIST:
            RM_FLIST+=line+' '
        
        subp.call('zip '+Output_Path+self.FILE+'.zip '+ZIP_FLIST,shell=True,executable='/bin/sh',cwd=os.getcwd())
        subp.call('mv '+self.FILE+' '+Output_Path,shell=True,executable='/bin/sh',cwd=os.getcwd())
        subp.call('rm '+RM_FLIST,shell=True,executable='/bin/sh',cwd=os.getcwd())
