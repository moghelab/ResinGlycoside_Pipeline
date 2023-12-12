from __future__ import division
import sys

'''
This is Version 3
more general in scope than the HortRes paper version, can handle different formats,
better annotation of script
Written for the Methods paper
GM, Oct 2023
'''

class step1:
    def __init__(self):
        pass
    def full_script(self,syslist):
        '''
        The first blocks read in information from INP1 and INP2
        '''
        #Read fragment file
        file1=open(syslist[0],'r')
        line1=file1.readline()
        fdict={}
        while line1:
            if line1.startswith('#'):
                pass
            else:        
                frags=line1.strip()#.split(',')
                #[float(x) for x in frags]
                fdict[frags]=1
            line1=file1.readline()
        file1.close()

        #Read all lines in INP2 (output of step0) first
        file1=open(syslist[1],'r').readlines()
        dictx={}
        for line in file1:
            sp=line.strip().split('\t')
            #Associate each line with the ScanID
            id1=sp[0]
            dictx[id1]=line.strip()
        file1=''        

        #Re-read INP2 line by line now
        file1=open(syslist[1],'r')
        out1=open(syslist[1]+".parsed",'w')
        line1=file1.readline()
        startProcessing=0; m=0
        while line1:
            #Header line
            if line1.startswith('Alignment'):
                startProcessing=1
                out1.write(line1)

            #Data lines
            elif startProcessing==1:
                '''
                The next few blocks filter the peaks based on area, fold-change-over-blank
                and other metrics
                '''
                tab1=line1.strip().split('\t')               

                #Specify values
                #Filter by peak area vs. blank
                id1=tab1[0]; postCuration=tab1[4]; flag=0; finStatus='fail'; hflag=0
                sn=float(tab1[6])
                
                #Remove in-source fragments defined by 'found in higher' tag
                if 'found in higher' in postCuration:
                    hflag=1               
                
                if hflag==0:                    
                    #Remove peaks without MS/MS 
                    if tab1[7]!='' and sn>=10:  #Signal to noise >10          
                        #Estimate if peak area of at least 1 sample is 5X max_blank
                        max_blank=int(tab1[8])
                        pa_samples=tab1[9].split(';')                        
                        pa_samples2=[int(x) for x in pa_samples]                        
                        
                        for area in pa_samples2:                
                            if area>1000:
                                if max_blank==0:
                                    flag=1
                                else:
                                    norm=area/max_blank
                                    if norm>5:
                                        flag=1                    
                                        
                #Only if the peak passes blank filter and area threshold in at least 1 of the samples
                if flag==1:
                    printstat=0
                    if syslist[1]=='34-Cent.txt' and id1=='7119':
                        printstat=1
                        
                    #Get reference peak intensity from the isotope file
                    plist=[]
                    #refpeaks=tab1[20].split(' ')                    
                    #for px in refpeaks:
                    #    plist.append(int(px.split(':')[1]))                    
                    for px in pa_samples:
                        plist.append(int(px))
                    maxpeak=max(plist)
                    
                    #Get MS/MS fragments and their intensities
                    #Get maximum intensity
                    msms=tab1[7].split(' ')                        
                    
                    #Initialize gdict and get fragments from fdict into gdict
                    #fdict is the fragment pairs we are looking for
                    gdict={}
                    for linex in fdict:                                        
                        frags=[float(x) for x in linex.split(',')]                        
                        for frag in frags:                        
                            #Initialize a dict to capture all intensities
                            if frag not in gdict:
                                gdict[frag]=[]
                    '''
                    The next few blocks gather the MSMS peak intensities of each fragment
                    and ask if the peak intensity is greater than 5% of the max peak intensity
                    of all fragments of the given peak
                    '''
                    #Gather all MSMS fragment intensities for fragment combinations listed in file 1
                    idict={}
                    for item in msms:
                        if item!='':
                            sp=item.split(':')
                            #print (id1, sp)
                            #print (msms)
                            ion=float(sp[0]); intensity=float(sp[1])
                            if intensity>0:
                                if intensity not in idict:
                                    idict[intensity]=[ion]
                                else:
                                    idict[intensity].append(ion)

                                #Now see if this MSMS fragment is present in the pairs we are looking for
                                for frag in gdict:                            
                                    if abs(frag-ion)<=0.01:    #Mass deviation                            
                                        gdict[frag].append(intensity)
                                    
                    #Get maximum ion intensity for each fragment in the desired MSMS fragment pair
                    hdict={}
                    for linex in fdict:
                        frags=[float(x) for x in linex.split(',')]                        
                        for frag in frags:  #These are all MSMS pairs we are looking for            
                            if gdict[frag]!=[]: #If the desired MSMS peak is present in the current peak
                                maxion=max(gdict[frag]) #v2
                                #maxion=sum(gdict[frag])
                                if frag not in hdict:
                                    if maxion>0.0:
                                        hdict[frag]=maxion                                
                        
                    #Get maximum intensity of ALL fragments
                    max_all=max(idict.keys())      #v2      
                    
                    #Filter using MS/MS fragment combinations
                    for linex in fdict:
                        frags=[float(x) for x in linex.split(',')]                        
                        count=0
                        for frag in frags:
                            if frag in hdict:
                                maxion=hdict[frag]
                                ratio=max_all/maxion #v2                                                                
                                if ratio<=20: #At least 5% of the max peak of this row                                
                                    count+=1                                    
                                    
                        #BOTH MS/MS peaks in the combination need to have
                        #intensities above 5% of max ion intensity
                        if count>=len(frags):
                            out1.write(line1)
                            m+=1
                            finStatus='pass'
            else:
                pass
            line1=file1.readline()
        file1.close(); out1.close()
        print ("  # of RG-like peaks in {}: {}".format(syslist[2],m))        

if __name__ == '__main__':
    step1 = step1()
    print ("INP1: Fragment combinations to filter using this script")
    print ("INP2: Peak area file")
    print ("INP3: Species name as Genus_species")
    fragfile=sys.argv[1]; areafile=sys.argv[2]; species=sys.argv[3]
    step1.full_script([fragfile,areafile, species])
    print ("Done!")
    


        
        
