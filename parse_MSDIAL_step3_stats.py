from __future__ import division
import sys, operator

class step3:
    def __init__(self):
        pass
    
    def stats(self,syslist):
        file1=open(syslist[0],'r')
        line1=file1.readline()
        flist=[]; fdict={}
        while line1:
            if line1.startswith('#'):
                pass
            else:        
                frag=float(line1.strip().split(',')[0])
                flist.append(frag)
                fdict[frag]=0
            line1=file1.readline()
        file1.close()

        file1=open(syslist[1],'r')
        out1=open(syslist[1]+".sigfrags",'w')
        line1=file1.readline()
        startProcessing=0; m=0; rcount=0; scount=0; oroot=0; oshoot=0; both=0; twofrags=0; threefrags=0; mulfrags=0; ambig=0
        sampledict={}; totallist=[]
        while line1:
            if line1.startswith('Alignment'):
                startProcessing=1
                samplenames=line1.strip().split('\t')[9].split(';')
                
            elif startProcessing==1:
                tab1=line1.strip().split('\t')
                #Specify values
                id1=tab1[0]; avgrt=tab1[1]; avgmz=tab1[2]; postCuration=tab1[4]; flag=0
                avgstr=f'{avgrt}|{avgmz}'

                #Get root vs shoot areas
                #pa_samples=[float(x) for x in tab1[28:30]]
                #for v in range(0,len(tab1)):
                #    print (v, tab1[v])
                pa_samples=[float(x) for x in tab1[9].split(';')] #for GNPS input
                for i in range(0,len(pa_samples)):
                    pa=pa_samples[i]; samplen=samplenames[i]                
                    if pa>0:
                        if id1 not in totallist:
                            totallist.append(id1)

                        #Get RG counts per sample
                        if samplen not in sampledict:
                            sampledict[samplen]=[id1]
                        else:
                            if id1 not in sampledict[samplen]:
                                sampledict[samplen].append(id1)
                
                '''
                shoot=pa_samples[0]; root=pa_samples[1]
                if root>10000:
                    rcount+=1
                if shoot>10000:
                    scount+=1

                #Get only root vs shoot
                if root>10000 and shoot==0:
                    oroot+=1
                elif root==0 and shoot>10000:
                    oshoot+=1
                else:
                    if root>10000 and shoot>10000:
                        both+=1
                    else:
                        ambig+=1
                '''

                #Parse MSMS fragments
                msms=tab1[7].split(' ')
                idict={}
                for item in msms:                    
                    if item!='':
                        sp=item.split(':')                
                        ion=float(sp[0]); intensity=float(sp[1])
                        #Combine intensities of all fragments of interest
                        #assuming a mass error of 0.05
                        for frag in flist:
                            if abs(frag-ion)<=0.05:                        
                                if intensity>0:
                                    if frag not in idict:
                                        idict[frag]=[intensity]
                                    else:
                                        idict[frag].append(intensity)                

                #Sum all intensities of a given fragment
                for frag in idict:
                    totint=sum(idict[frag])
                    idict[frag]=totint
                

                #Get fragment with the maximum total intensity     
                #tup=max(idict.iteritems(), key=operator.itemgetter(1)) #Python2
                tup=max(idict.items(), key=operator.itemgetter(1)) #Python3                
                maxfrag=tup[0]; maxint=tup[1]        
                select=[maxfrag]

                #Get other signature fragments that are atleast 5% maxfrag
                for frag in idict:
                    if frag!=maxfrag:
                        intensity=idict[frag]
                        if (intensity/maxint)<=0.05:
                            select.append(frag)

                #Write this information to output
                v1=','.join([str(x) for x in select])
                out1.write(f'{avgstr}\t{v1}\n')                
                if len(select)==2:
                    twofrags+=1
                elif len(select)==3:
                    threefrags+=1
                elif len(select)>3:
                    mulfrags+=1                
                
                
                #Count peaks of each type. Are any signature -OH-FA fragments in there?
                #flist doesnt include FA+sugars
                for frag in flist:            
                    if frag in select:                        
                        fdict[frag]+=1                    
            else:
                pass
                    
            line1=file1.readline()
        file1.close(); out1.close()

        #print (flist)
        #sys.exit()

        out1=open(syslist[1]+"."+syslist[3]+".stats",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))

        #Make headers for different types of fragments and different samples
        flist2=[str(x) for x in flist]
        fline='\t'.join(flist2)
        samplelist=list(sampledict.keys())
        sampleline='\t'.join(samplelist)
        
        #out1.write('#FileName\tSpecies\tTotal\tTotal2\tShoot\tRoot\tOnlyShoot\tOnlyRoot\tBoth\tAmbig\t{}\tTwoFrags\tThreeFrags\tMulFrags\n'.format('\t'.join(flist2)))
        out1.write('#FileName\tSpecies\tTotal\tTotal2\t{}\t{}\tTwoFrags\tThreeFrags\tMulFrags\n'.format(sampleline,fline))
        
        #Count # of fragments of each long acyl chain
        fx=[]; nfcount=0
        for frag in flist:
            fcount=fdict[frag]
            nfcount+=fcount
            fx.append(str(fcount))
        fline='\t'.join(fx)
        #total=both+oroot+oshoot

        #Get total RGs identified using two different means
        total=len(totallist)
        total2=nfcount-(twofrags+(threefrags*2))

        #Get RG counts on a per sample basis
        rglist=[]
        for samplex in samplelist:
            rgcount=len(sampledict[samplex])
            rglist.append(str(rgcount))
        rgline='\t'.join(rglist)
        
        #out1.write(f'{syslist[1]}\t{syslist[2]}\t{total}\t{total2}\t{scount}\t{rcount}\t{oshoot}\t{oroot}\t{both}\t{ambig}\t{fline}\t{twofrags}\t{threefrags}\t{mulfrags}\n')
        out1.write(f'{syslist[1]}\t{syslist[2]}\t{total}\t{total2}\t{rgline}\t{fline}\t{twofrags}\t{threefrags}\t{mulfrags}\n')
        out1.close()
        print ("  Final count: ", total)

if __name__ == '__main__':
    step3 = step3()
    print ("INP1: Frag. combinations file, where the first column is the long-acyl-fragment")
    print ("INP2: Output of Step1 (*.parsed)")
    print ("INP3: Species name as Genus_species")
    print ("INP4: Suffix for the output file of stats")
    fragfile=sys.argv[1]; areafile=sys.argv[2]; species=sys.argv[3]
    sufx=sys.argv[4]
    step3.stats([fragfile,areafile,species,sufx])
    print ("Done!")


    
                    
                    
                        

