from __future__ import division
import sys, operator, random
import pandas as pd

class rgid_step1:
    def __init__(self):
        pass

    def full_script (self, tab):
        #Parse inputs
        filelist=tab[0] #list of files to analyze
        file1=open(tab[1],'r')  #motiffile        
        desfrags=int(tab[2])
        sufs=tab[3]
        
        print ("  -->Only top {} MSMS fragments by intensity per peak will be analyzed<--".format(desfrags))        
        line1=file1.readline()
        motifdict={}; longd={}; smalld={}; sugard={}; otherd={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                grp=tab1[0]; mass=float(tab1[1]); annot=tab1[2]
                if grp.startswith('comb'):
                    #mass=mass
                    pass
                else:
                    if grp=='longacyl':
                        nmass=mass-1.007
                        longd[nmass]=annot
                    elif grp=='sugar':
                        sugard[mass]=annot
                    elif grp=='smallacyls':
                        nmass=mass-1.007
                        smalld[nmass]=annot
                    elif grp=='other':
                        nmass=mass-1.007
                        otherd[nmass]=annot
                        

                mass2=mass-1.007
                if mass2 not in motifdict:
                    motifdict[mass2]=annot
                else:
                    print ("Mass repeat: ", mass)
            line1=file1.readline()
        file1.close()
        motiflist=list(motifdict.keys())
        mlongs=list(longd.keys())
        msmalls=list(smalld.keys())
        mothers=list(otherd.keys())

        #Make all combinations of masses for each longd
        longsugard={}
        for long in longd:
            longsugard[long]={}
            for sugar in sugard:
                longsugar=long+sugar-18.01
                str1='{}+{}'.format(longd[long],sugard[sugar])
                if longsugar not in longsugard[long]:
                    longsugard[long][longsugar]=str1
                else:
                    print ("Repeat: ", long, longsugar, str1)
                    sys.exit()
        mlongsugars=list(longsugard.keys())

        #Make all combinations of masses for each longsugar
        longsugarsugard={}
        for long in longsugard:
            longsugarsugard[long]={}
            
            for longsugar in longsugard[long]:
                longsugarsugard[long][longsugar]={}
                
                for sugar in sugard:
                    m1=longsugar+sugar-18.01
                    m2=longsugar+sugar-18.01-18.01 #circular

                    #Assign name
                    longsugarname=longsugard[long][longsugar]
                    str1='{}+{}'.format(longsugarname,sugard[sugar])
                    str2='{}+{}-Water'.format(longsugarname,sugard[sugar])
                    #print (type(long), type(longsugar), type(m1), type(m2), type(longsugard))
                    
                    if m1 not in longsugarsugard[long][longsugar]:
                        longsugarsugard[long][longsugar][m1]=str1
                    else:
                        print ("Repeat2: ", long, longsugar, m1, str1)
                        sys.exit()

                    if m2 not in longsugarsugard[long][longsugar]:
                        longsugarsugard[long][longsugar][m2]=str2
                    else:
                        print ("Repeat3: ", long, longsugar, m2, str2)
                        sys.exit()
        mlongsugarsugars=list(longsugarsugard.keys())           

        file1=open(filelist,'r') #List of files to analyze
        out1=open(filelist+".nid1",'w')        
        out1.write('#FileID\tSeriesName\tScanID\tMaxPeaks\tIdentifiedPeaks\tUnexplainedPeaks\n')
        line1=file1.readline()

        ############
        def getpeaki(pdictx,scanx,peakx):
            try:
                intensex=pdict[scan][peak]
            except:
                #First MS1 peak doesn't have any intensity value associated with it
                intensex=0.0
            return intensex
        ############
            
        while line1:
            if line1.startswith('#'):
                pass
            else:
                sp=line1.strip().split('\t')

                #Rename the file to account for different fragpairs files
                fcount=sp[0].zfill(2)
                areafile=sp[0]; mgffile=sp[-1]; specname=sp[1]
                dict1={}; dict2={}; mdict={}; scandict={}; pdict={}
                for suf in sufs:
                    fn=f'2_MGF/{mgffile}.{suf}.fil.mgf'
                    #os.system('cp {} {}.{}\n'.format(fnx,fnx,suf))
                    print ("  Now reading: ", fn)
                    file2=open(fn,'r')
                    line2=file2.readline()
                    start=0; count=0
                    dict1[suf]={}
                    while line2:
                        if line2.startswith('BEGIN IONS'):
                            start=1; count+=1; peak271=0
                        elif line2.startswith('END IONS'):
                            start=0
                            #print ("######")
                            if count==10:
                                #sys.exit()
                                pass
                            
                            #Get top N ions in the dict sorted by their intensities
                            sorted_d = sorted(msdict.items(), key=operator.itemgetter(1))
                            sorted_d.reverse()

                            if len(sorted_d)>=desfrags:
                                ftuples=sorted_d[0:desfrags]
                            else:
                                ftuples=sorted_d

                            for str1 in ftuples:
                                #print ("file-{}".format(fcount), scan, pepmass, str1)
                                mass=list(str1)[0]; intensity=int(list(str1)[1])
                                pdict[scan][mass]=intensity
                                if intensity>1000:
                                    dict1[suf][scan].append(mass)
                                    if mass not in mdict:
                                        mdict[mass]=1
                                
                            msdict={}
                            #sys.exit()
                            
                        elif line2.strip()=='':
                            pass
                        else:
                            if start==1:
                                if '=' in line2:
                                    #Get MS1 peak information
                                    if 'PEPMASS' in line2:
                                        ms1=float(line2.strip().split('=')[1])
                                        pepmass=ms1
                                        dict1[suf][scan].append(ms1)
                                        if ms1 not in mdict:
                                            mdict[ms1]=1
                                            
                                    #Get Scan information
                                    elif 'SCANS' in line2:
                                        scan=int(line2.strip().split('=')[1])
                                        dict1[suf][scan]=[]
                                        scandict[scan]=1
                                        if scan not in pdict:
                                            pdict[scan]={}
                                        else:
                                            print ("    Scan repeat. Second instance will be excluded: ", scan, suf)
                                            
                                        #Initialize msdict for collecting msms to get top 10
                                        msdict={} 
                                    else:
                                        pass
                                else:
                                    #Get MSMS peak information
                                    sp=line2.strip().split()
                                    v1=float(sp[0]); v2=int(sp[1])
                                    if int(v1)==271:
                                        peak271=1
                                    if v1 not in msdict:
                                        msdict[v1]=v2
                                    else:
                                        if v2 > msdict[v1]:
                                            msdict[v1]=v2
                        line2=file2.readline()
                    file2.close()
                    #print ("Finished: ", fn)

                #Get all peaks list
                mlist=sorted(list(mdict.keys()))
                scanlist=sorted(list(scandict.keys()))
                sline='\t'.join([str(f) for f in scanlist])

                #Combine all masses from all 3 scans into one file
                dict2={}; m=0
                for suf in dict1:
                    for scan in dict1[suf]:
                        if scan not in dict2:
                            dict2[scan]=dict1[suf][scan] #peaks
                        else:
                            print ("    Scan repeat. Second instance will be excluded: ", fn, scan)
                            
                #For each scan classify the peaks based on motifs
                #print (pdict)
                #sys.exit()
                for scan in dict2:
                    peaks=dict2[scan]; useful=[]; unexplained=[]; donepeak=[]
                    maxpeak=max(peaks)

                    #Identify longacyls and sugars hierarchically
                    flonglist=[]; flonglist2=[]
                    flongsugarlist=[]; flongsugarlist2=[]
                    flongsugarsugarlist=[]; flongsugarsugarlist2=[]

                    #Get the longacyl peaks first
                    for peak in peaks:
                        flag=0
                        for motif in mlongs:
                            if peak not in donepeak:
                                diffpeak=abs(peak-motif)
                                if diffpeak<0.04:
                                    flonglist.append(motif)
                                    intense=getpeaki(pdict,scan,peak)
                                    str1='{}|{}|{}|{}'.format(peak, motif, intense, longd[motif])
                                    flonglist2.append(str1)
                                    donepeak.append(peak)

                    #Get longsugar peaks
                    if flonglist!=[]:
                        for long in flonglist:
                            longsugars=list(longsugard[long].keys())
                            
                            for peak in peaks:
                                for motif in longsugars:
                                    if peak not in donepeak:
                                        diffpeak=abs(peak-motif)
                                        if diffpeak<0.04:
                                            flongsugarlist.append(motif)
                                            intense=getpeaki(pdict,scan,peak)
                                            annot=longsugard[long][motif]
                                            str1='{}|{}|{}|{}'.format(peak, motif, intense, annot)
                                            flongsugarlist2.append(str1)
                                            donepeak.append(peak)
                                    
                    #Get longsugarsugar peaks
                    if flongsugarlist!=[]:
                        for long in flonglist:
                            nflongsugarlist=list(longsugard[long].keys())
                            
                            for longsugar in nflongsugarlist:
                                #print (">>>>", long, longsugar)
                                longsugarsugars=list(longsugarsugard[long][longsugar].keys())
                                
                                for peak in peaks:
                                    for motif in longsugarsugars:
                                        if peak not in donepeak:
                                            diffpeak=abs(peak-motif)
                                            if diffpeak<0.04:
                                                flongsugarsugarlist.append(motif)
                                                intense=getpeaki(pdict,scan,peak)
                                                annot=longsugarsugard[long][longsugar][motif]
                                                str1='{}|{}|{}|{}'.format(peak, motif, intense, annot)
                                                flongsugarsugarlist2.append(str1)
                                                donepeak.append(peak)

                    #Combine all long peaks
                    longcombined=flonglist+flongsugarlist+flongsugarsugarlist
                    longcombined2=flonglist2+flongsugarlist2+flongsugarsugarlist2

                    #Get other peaks
                    otherpeaks=[]
                    for peak in peaks:
                        for motif in mothers:
                            diffpeak=abs(peak-motif)
                            if diffpeak<0.04:
                                intense=getpeaki(pdict,scan,peak)
                                str1='{}|{}|{}|{}'.format(peak, motif, intense, motifdict[motif])
                                otherpeaks.append(str1); donepeak.append(peak)
                                

                    #Get acyl peaks
                    useful=[]; unexplained=[]
                    for peak in peaks:
                        flag=0
                        for motif in msmalls:
                            diffpeak=abs(peak-motif)
                            if diffpeak<0.04:
                                intense=getpeaki(pdict,scan,peak)
                                
                                #First differentiate between pentanoic and acetoacetic masses
                                #acetoacetic is likely hexose degradation
                                if round(motif,2)-101.06==0 or round(motif,2)-101.02==0:
                                    if diffpeak<0.02:
                                        str1='{}|{}|{}|{}'.format(peak, motif, intense, motifdict[motif])
                                        useful.append(str1); donepeak.append(peak)
                                        flag=1
                                        
                                #then between caffeic acid and Hexose
                                elif round(motif,2)-179.05==0 or round(motif,2)-179.03==0:
                                    if diffpeak<0.01:
                                        str1='{}|{}|{}|{}'.format(peak, motif, intense, motifdict[motif])
                                        useful.append(str1); donepeak.append(peak)
                                        flag=1

                                #then between coumaric and deoxyhexose
                                elif round(motif,2)-163.06==0 or round(motif,2)-163.04==0:
                                    if diffpeak<0.01:
                                        str1='{}|{}|{}|{}'.format(peak, motif, intense, motifdict[motif])
                                        useful.append(str1); donepeak.append(peak)
                                        flag=1
                                        
                                else:
                                    str1='{}|{}|{}|{}'.format(peak, motif, intense, motifdict[motif])
                                    useful.append(str1); donepeak.append(peak)
                                    flag=1

                        if flag==0:
                            if peak not in donepeak:
                                intense=getpeaki(pdict,scan,peak)
                                str2='{}|{}'.format(peak, intense)
                                unexplained.append(str2)
                            

                    #Write to output
                    nuseful=longcombined2+useful+otherpeaks
                    out1.write('file-{}\t{}\t{}\t{}\t{}\t{}\n'. \
                               format(fcount,specname,scan,maxpeak,nuseful,unexplained))
                print ("  Finished with: ", areafile, mgffile)
                print ("###############")
            line1=file1.readline()
        file1.close(); out1.close()        
        
if __name__ == '__main__':
    id_step1=id_step1()    
    print ("INP1: filelist.tab")    
    flist=sys.arv[1]
    suffixez=['mz132','mz146','mz162']        
    mfile='3_predictRGs/motifs.groups'
    dfrags=10
    rgid_step1.full_script([flist,mfile,dfrags,suffixez])      



        
            
