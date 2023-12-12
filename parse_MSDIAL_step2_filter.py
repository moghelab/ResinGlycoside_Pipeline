from __future__ import division
import sys, operator
import statistics as st

'''
This is Version 3
more general in scope than the HortRes paper version, can handle different formats,
better annotation of script
Written for the Methods paper
GM, Oct 2023

This script removes lines containing "found in higher" tags as well as selects better of
the two lines containing adducts
'''

class step2:
    def __init__(self):
        pass

    def mergeLines(self,fn1):
        #Read all lines in the file first
        file1=open(fn1,'r').readlines()
        dict1={}
        for line in file1:
            sp=line.split('\t')
            id1=sp[0]; pc=sp[5]            
            dict1[id1]=line.strip()
        file1=''

        #Now read file line by line
        file1=open(fn1,'r')        
        out1=open(fn1+".fil",'w')
        line1=file1.readline()
        removedLine=[]
        while line1:
            if line1.startswith('Alignment ID'):
                out1.write(line1)                
            else:
                tab1=line1.strip().split('\t')
                id1=tab1[0]; pc=tab1[4]
                if id1 not in removedLine:
                    #Get peak areas ORIGINAL            
                    org1=int(tab1[8]) #Blank
                    org2list=[float(x) for x in tab1[9].split(';')] #Samples
                    org2=max(org2list) #Maximum peak area of the samples

                    #Adduct linked lines
                    sp=pc.split('; ')
                    if 'adduct linked' in pc: 
                        hflag=0                        
                        if 'found in higher' in pc: #mostly taking care of in-source frag                           
                            for item in sp:
                                if item.startswith('found in higher'):
                                    idx=item.split('_')[1]
                                    if idx in dict1:
                                        hflag=1                        

                        #hflag=0 means found in higher ID is not present in the fil file
                        #so now we can proceed to assess the adducts
                        if hflag==0:                         
                            ilist=[]; flag=0
                            for item in sp:
                                if item.startswith('adduct'):
                                    sp2=item.split('_')
                                    adduct=sp2[1]; idx=sp2[0].split()[-1]
                                    if idx in dict1: #Dictionary of all IDs in this file
                                        flag=1
                                        
                                        #Choose line where both organs have reads
                                        iline=dict1[idx]
                                        tab2=iline.strip().split('\t')
                                        id2=tab2[0]; pc2=tab2[4]; org1x=int(tab2[8])
                                        org2xlist=[float(x) for x in tab2[9].split(';')]
                                        org2x=max(org2xlist)

                                        #org1x is blank; org2x is max(samples)
                                        if org1>org1x and org2>org2x:
                                            out1.write(line1)  #select new line                                           
                                            removedLine.append(id2)
                                        elif org1<org1x and org2<org2x:
                                            out1.write(iline)  #select old line                                          
                                            removedLine.append(id1)
                                        elif org1>org1x and org2>0:
                                            out1.write(line1)  #select new line                                            
                                            removedLine.append(id2)
                                        elif org1>0 and org2>org2x:
                                            out1.write(line1)  #select new line                                            
                                            removedLine.append(id2)
                                        elif org1==0 or org2==0:
                                            out1.write(line1)  #select new line                                            
                                            removedLine.append(id2)
                                        else:
                                            print (line1)
                                            print (org1,org1x)
                                            print (org2,org2x)

                            #If none of the linked adduct IDs are there in the input file
                            #write this line to output
                            if flag==0:                                
                                out1.write(line1)
                        else:
                            pass
                    else:
                        hflag=0                        
                        if 'found in higher' in pc:                            
                            for item in sp:
                                if item.startswith('found in higher'):
                                    idx=item.split('_')[1]
                                    if idx in dict1:
                                        hflag=1
                        
                        if hflag==0:                            
                            out1.write(line1)
            line1=file1.readline()
        file1.close()       
        print ("  IDs removed due to adducts/in-source: ", len(list(set(removedLine))))

    
if __name__ == '__main__':
    step2=step2()
    print ("INP1: Output of Step1 (*.parsed)")    
    areafile=sys.argv[1]
    step2.mergeLines(areafile)
    #parseOut.filter([fragfile,areafile,species])
    print ("Done!")
