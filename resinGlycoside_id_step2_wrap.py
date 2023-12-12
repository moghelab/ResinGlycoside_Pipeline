from __future__ import division
import sys, operator, random
import pandas as pd

class rgid_step2:
    def __init__(self):
        pass

    def full_script (self, tab):
        #Parse inputs
        file1=open(tab[1],'r') #motiffile
        s1fn=tab[0] #output of step1
        
        line1=file1.readline()
        motifdict={}; motifdict2={}
        while line1:
            if line1.startswith('#') or line1.startswith('comb'):
                pass
            else:
                tab1=line1.strip().split('\t')
                grp=tab1[0]; mass=float(tab1[1]); annot=tab1[2]; flag=0
                if grp.startswith('smallacyls'):
                    name=annot.split('$')[1]; flag=1
                elif grp.startswith('sugar'):
                    name=annot; flag=1

                if flag==1:
                    neutLoss=mass-18.01; nx=int(round(neutLoss))            
                    if nx not in motifdict:                
                        motifdict[nx]='{}-${}-loss'.format(neutLoss, name)
                        
                    neutLoss2=mass; nx2=int(round(neutLoss2))
                    str1='{}-${}-$Water-loss'.format(neutLoss2, name)            
                    if nx2 not in motifdict2:                
                        motifdict2[nx2]='{}-${}-loss'.format(neutLoss2, name)
                        
                else:
                    #print ("Mass repeat: ", mass)
                    print ("  Skipped: ", line1.strip())
            line1=file1.readline()
        file1.close()
        motiflist=list(motifdict.keys())        

        file1=open(s1fn,'r')
        out1=open(s1fn+".id2",'w')        
        out1.write('#FileID\tSeriesName\tScanID\tMaxPeaks\tIdentifiedPeaks\tNeutralLosses\tUnexplainedPeaks\n')
        out2=open(s1fn+".id2.tab",'w')        
        out2.write('#FileID\tSeriesName\tScanID\tMaxPeaks\tInfo1\tInfo2\tInfo3\n')
        line1=file1.readline()
        pdict={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                fn=tab1[0]; specs=tab1[1]
                scan=tab1[2]; mp=float(tab1[3])
                idd=eval(tab1[4])
                un=eval(tab1[5])

                #Make a peaklist
                plist1=[mp]; plist2=[]; xdict={}
                for peak in idd:
                    sp=peak.split('|')
                    p=float(sp[0])
                    plist1.append(p)
                    xdict[p]=peak
                    
                for peak in un:
                    sp=peak.split('|')
                    p=float(sp[0])
                    plist2.append(p)
                    xdict[p]=peak
                plist=plist1+plist2

                #Obtain info about neutral losses
                nlist=[]; slist1=[]; nlistx=[]
                for i in range(0,len(plist)):
                    peak1=plist[i]
                    for j in range(i+1,len(plist)):
                        peak2=plist[j]
                        if peak1 in plist2 or peak2 in plist2:
                            if peak1>300.0 and peak2>300.0: #Losses calculated only if peaks>300 m/z                        
                                diff=(abs(peak1-peak2)); flag=0
                                #print (peak1, peak2, diff)

                                #First block of ifs                    
                                if abs(46-diff)<0.05:
                                    annot='loss46|{:.4f}-$Formate-loss'.format(diff); flag=1                            
                                elif abs(43.99-diff)<=0.03:
                                    annot='loss44|{:.4f}-$CO2-loss'.format(diff); flag=1
                                elif abs(18-diff)<=0.03:
                                    annot='loss18|{:.4f}-$Water-loss'.format(diff); flag=1
                                elif abs(56-diff)<=0.05:
                                    annot='loss56|{:.4f}-$Something?-loss'.format(diff); flag=1
                                elif abs(84.05-diff)<=0.02:
                                    annot='loss84|{:.4f}-$Pentanoyl-loss'.format(diff); flag=1
                                elif abs(102.055-diff)<=0.02: #actual diff=102.08
                                    annot='loss102|{:.4f}-$Pentanoyl$water-loss'.format(diff); flag=1
                                #elif abs(42-diff)<=0.03:
                                #    annot='loss42|{:.4f}-$Acetyl-loss'.format(diff); flag=1
                                elif int(diff)==0:
                                    annot='loss0|{:.4f}-$ParentIon-loss'.format(diff); flag=1
                                else:
                                    pass

                                #if flag==1:
                                #    print (peak1, peak2, diff, annot)
                                

                                rdiff=int(round(diff))                        
                                if rdiff in motifdict:
                                    val=motifdict[rdiff]
                                    expDiff=float(val.split('-')[0])
                                    if abs(diff-expDiff)<0.02:
                                        annot='loss{}|Actual:{:.4f}-Pred:{}'.format(rdiff,diff,val)
                                        flag=1
                                else:                            
                                    if rdiff in motifdict2:
                                        val=motifdict2[rdiff] #Actual mass
                                        expDiff=float(val.split('-')[0])
                                        if abs(diff-expDiff)<0.02:
                                            annot='loss{}|Actual:{:.4f}-Pred:{}'.format(rdiff,diff,val)
                                            flag=1
                                    else:
                                        pass
          
                                if flag==1:
                                    str1='{}-{}|{}'.format(max([peak1,peak2]),min([peak1,peak2]),annot)
                                    nlist.append(str1); slist1.append(peak1); slist1.append(peak2)                
                        else:
                            pass

                #Annotate the 125.0241 fragment
                for i in range(0,len(plist)):
                    peak1=plist[i]            
                    if abs(peak1-125.0239)<0.04:                
                        annot='sugardeg|Actual:{:.4f}-Pred:125.0239'.format(peak1)
                        str1='{}\t{}'.format(peak1,annot)
                        nlistx.append(str1); slist1.append(peak1)            

                #Write to output
                slist2=[]
                for peak in plist2:
                    if peak not in slist1:                
                        str1=xdict[peak]
                        slist2.append(str1)

                #Write to output
                out1.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(fn,specs,scan,mp,idd,nlist,slist2))
                if fn not in pdict:
                    pdict[fn]=[scan]
                else:
                    if scan not in pdict[fn]:
                        pdict[fn].append(scan)
                    
                for item in idd:
                    sp=item.split('|')
                    n1=sp[0]; n2='|'.join([sp[0],sp[1],sp[2]]); n3=sp[3]
                    out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(fn,specs,scan,mp,n1,n2,n3))
                for item in nlist:
                    sp=item.split('|')
                    m1=sp[0]; m2='|'.join([sp[1],sp[2]])
                    out2.write('{}\t{}\t{}\t{}\t{}\t{}\tNeutralLoss\n'.format(fn,specs,scan,mp,m1,m2))
                for item in nlistx:
                    sp=item.split('\t')
                    m1=sp[0]; m2=sp[1]
                    out2.write('{}\t{}\t{}\t{}\t{}\t{}\tSugarDegradation\n'.format(fn,specs,scan,mp,m1,m2))
                for item in slist2:
                    sp=item.split('|')
                    d1=sp[0]; d2=sp[1]
                    out2.write('{}\t{}\t{}\t{}\t{}|{}\tNA\tUnmatchedFrags\n'.format(fn,specs,scan,mp,d1,d2))
                out2.write('######\n')
            line1=file1.readline()
        file1.close(); out1.close(); out2.close()

        for fname in pdict:
            print (f"# of peaks annotated for {fname}: {len(pdict[fname])}")            
        
if __name__ == '__main__':
    rgid_step2=rgid_step2()    
    print ("INP1: Output of step1")    
    flist=sys.arv[1]    
    mfile='3_predictRGs/motifs.groups'    
    rgid_step2.full_script([flist,mfile])
    




        
            
