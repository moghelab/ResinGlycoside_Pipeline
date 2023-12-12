from __future__ import division
import sys
'''
This is Version 3
more general in scope than the HortRes paper version, can handle different formats,
better annotation of script
Written for the Methods paper
GM, Oct 2023
'''

#Reorganizes all files

class step0:
    def __init__(self):
        pass

    def writeTabs(self, tab, mtabs):
        #Extract values from the file
        scan=tab[0]; rt=tab[1]; mz=tab[2]; adduct=tab[4]
        samples=mtabs[0].split(','); blanks=mtabs[1].split(',')
        
        pcindex=int(mtabs[2])
        fillper=int(mtabs[3]); sn=int(mtabs[4]); msms=int(mtabs[5])        
            
        #Remake line
        nlist=[scan, rt, mz, adduct] #0,1,2,3
        nlist.append(tab[pcindex]) #4
        nlist.append(tab[fillper]) #5
        nlist.append(tab[sn])      #6
        nlist.append(tab[msms])    #7    

        #Get max for blank
        blist=[]
        for item in blanks:
            blist.append(tab[int(item)])
        bmax=max(blist)            #8
        nlist.append(str(bmax))  

        #Get sample values
        slist=[]
        for item in samples:
            slist.append(tab[int(item)])
        str1=';'.join(slist)
        nlist.append(str1)         #9

        #Make line
        mynline='\t'.join(nlist)
        return (mynline)
        
    def full_script(self, syslist):
        #Read inputs
        file1=open(syslist[0],'r')
        tabslist=syslist[1]

        #Create outout file
        out1=open(syslist[0]+".reorder",'w')
        line1=file1.readline()
        startProcessing=0
        while line1:
            #Header line
            if line1.startswith('Alignment'):
                startProcessing=1
                tab1=line1.strip().split('\t')
                nline=self.writeTabs(tab1, tabslist)
                out1.write(f'{nline}\n')

            #Non-header line
            elif startProcessing==1:
                tab1=line1.strip().split('\t')
                tab1=line1.strip().split('\t')
                nline=self.writeTabs(tab1, tabslist)
                out1.write(f'{nline}\n')

            #Blank lines at the beginning
            else:
                pass
            line1=file1.readline()
        file1.close(); out1.close()

if __name__ == '__main__':
    step0 = step0()    
    print ("INP1: Peak area file")
    print ("INP2: List of indexes with specific data")
    print ("#PeakAreaFile-name      Species-name    Samples-index   Blank-index" \
           "Post-curation-index     Fillper-index   S/N-index       MSMS-index  GnpsMGF-file")
    areafile=sys.argv[2]
    datatabs=eval(sys.argv[3])
    step0.full_script([areafile, datatabs])
    print ("Done!")
    

