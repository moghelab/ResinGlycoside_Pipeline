from __future__ import division
import sys, operator

class step4:
    def __init__(self):
        pass
    
    def filmgf(self,syslist):        
        file1=open(syslist[0],'r')
        line1=file1.readline()
        flist=[]; fdict={}
        while line1:
            if line1.startswith('#') or line1.startswith('Alignment'):
                pass
            else:        
                alnid=(line1.strip().split('\t')[0])
                flist.append(alnid)
                fdict[alnid]=0
            line1=file1.readline()
        file1.close()

        #print ("Peaks in fdict: ", len(fdict.keys()))
        #print (fdict)
        #sys.exit()

        file1=open(syslist[1],'r')
        out1=open(syslist[1]+f".{syslist[2]}.fil.mgf",'w')
        line1=file1.readline()
        startProcessing=0; m=0; rcount=0; goahead=0
        while line1:
            if line1.startswith('BEGIN IONS'):
                startProcessing=1
                
            elif line1.startswith('END IONS'):
                if goahead==1:
                    out1.write('{}\n\n'.format(line1.strip()))
                startProcessing=0; goahead=0
                
            elif startProcessing==1:
                if goahead==0:
                    element=line1.strip().split('=')
                    if element[0]=='SCANS':
                        eid=element[1]
                        if eid in fdict:
                            goahead=1
                            out1.write('BEGIN IONS\n')
                            out1.write(line1)
                            rcount+=1
                        else:
                            goahead=0
                elif goahead==1:
                    out1.write(line1)
            else:
                pass
            line1=file1.readline()
        file1.close(); out1.close()
        print ("  # of SCANS in parsed.fil.uniq: ", len(fdict.keys()))
        print ("  # of SCANS in OUT: ", rcount)
        #sys.exit()

if __name__ == '__main__':
    step3 = step3()
    print ("INP1: parsed.fil.uniq output of step2")
    print ("INP2: MGF file")
    print ("INP3: Suffix to add to output file")
    fragfile=sys.argv[1]; mgffile=sys.argv[2]; sufx=sys.argv[3]
    step3.stats([fragfile,mgffile,sufx])
    print ("Done!")


    
                    
                    
                        
