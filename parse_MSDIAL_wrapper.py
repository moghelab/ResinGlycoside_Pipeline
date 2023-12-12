import sys, os
import parse_MSDIAL_step0_reorder
import parse_MSDIAL_step1_gnps
import parse_MSDIAL_step2_filter
import parse_MSDIAL_step3_stats
import parse_MSDIAL_step4_MGFfilter
import resinGlycoside_id_step1_wrap
import resinGlycoside_id_step2_wrap

step0=parse_MSDIAL_step0_reorder.step0()
step1=parse_MSDIAL_step1_gnps.step1()
step2=parse_MSDIAL_step2_filter.step2()
step3=parse_MSDIAL_step3_stats.step3()
step4=parse_MSDIAL_step4_MGFfilter.step4()
rgid_step1=resinGlycoside_id_step1_wrap.rgid_step1()
rgid_step2=resinGlycoside_id_step2_wrap.rgid_step2()

def help1():    
    print ("---------- PRE-SCRIPT INSTRUCTIONS -----------")
    print ("* Make four sub-folders called 0_RgFrags, 1_PeakArea, 2_MGF, 3_predictMotifs")
    print ("       mkdir 0_RgFrags 1_PeakArea 2_MGF 3_predictMotifs")
    print ("* Put mz132, mz146, mz162 frag pairs in 0_RgFrags")
    print ("* Put all Area files in 1_PeakArea")
    print ("* Put all GnpsMgf files in 2_MGF")
    print ("* Put motifs.groups in 3_predictMotifs")
    print ("* Ensure all scripts are in this working directory")
    print ("---------- REQUIRED ARGUMENTS -----------")    
    print ("-filelist: File containing a list of files to process")
    
    print ("---------- INSTRUCTIONS FOR HOW TO FORMAT -filelist -----------")
    print ("      Tab-delimited file should have following columns\n" \
       "#PeakAreaFile-name      Species-name    Samples-index   Blank-index" \
       "Post-curation-index     Fillper-index   S/N-index       MSMS-index  GnpsMGF-file")
    print ("      Columns should start with index0 and should be csv")
    print ("---------------------------------------------------------------")
    sys.exit()

filelist=""
for i in range(1,len(sys.argv)):    
    if sys.argv[i]=='-filelist':
        filelist=sys.argv[i+1]
        
if filelist=="":
    help1()
else:
    suffixes=['mz132','mz146','mz162']

    #MODULE 1: PARSE THE MS-DIAL OUTPUT FILES
    for suf in suffixes:
        print ("###############")
        print (f"**Extracting all peaks with {suf} fragment combinations...")
        file1=open(filelist,'r')
        fragfile=f'0_RgFrags/frags.list.{suf}.pairs'        
        out1=open(sys.argv[2]+f".{suf}.combined.stats",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        
        line1=file1.readline()
        added=0
        while line1:
            if line1.startswith('#'):
                pass
            else:
                sp=line1.strip().split('\t')                
                datatabs=sp[2:]
                mgffile=sp[-1]
                print ("PROCESSING ===================> ", sp[0])
                
                #Rename the file to account for different fragpairs files        
                fnx=f'1_PeakArea/{sp[0]}'
                os.system('cp {} {}.{}\n'.format(fnx,fnx,suf))
                fn=f'{fnx}.{suf}'

                #Reorganize table
                species=sp[1]
                print (f"Processing {species} through file {fn}")
                print ("Step 0...Reformat...")                
                step0.full_script([fn, datatabs])                        
                
                #Run pipeline
                print ("Step 1...Extract RG-like peaks...")
                step1.full_script([fragfile, fn+".reorder", species])
                
                print ("Step 2...Filter...")
                step2.mergeLines(fn+".reorder.parsed")        

                #Remove duplicated lines
                lines=open(fn+".reorder.parsed.fil",'r')
                outx=open(fn+".reorder.parsed.fil.uniq",'w')
                doneline=[]
                for line in lines:
                    if line not in doneline:
                        outx.write(line)
                        doneline.append(line)
                outx.close()

                #Make a readable output
                lines=open(fn+".reorder.parsed.fil.uniq",'r').readlines()
                outx=open(fn+".reorder.parsed.fil.uniq.tab",'w')        
                for line in lines:
                    tab=line.strip().split('\t')
                    fsp=tab[-1].split(';')
                    oline='\t'.join(tab[0:-1])
                    fline='\t'.join(fsp)
                    outx.write('{}\n'.format('\t'.join([oline, fline])))
                outx.close()

                print ("Step 3...Make individual stats files...")
                step3.stats([fragfile, fn+".reorder.parsed.fil.uniq", species,suf])

                print ("Step 4...Process mgf files")
                step4.filmgf([fn+".reorder.parsed.fil.uniq", f'2_MGF/{mgffile}', suf])        

                #Get stats based on filtering of GNPS table file
                print ("Step 5...Combine stats files")
                file2=open(fn+".reorder.parsed.fil.uniq.{}.stats".format(suf),'r')
                line2=file2.readline()        
                while line2:
                    if line2.startswith('#python'):
                        pass
                    else:                        
                        out1.write(line2)                    
                    line2=file2.readline()
                out1.write('======================\n')
                file2.close()
                print (f"  See output: {fn}.reorder.parsed.fil.uniq.{suf}.stats")
                
            line1=file1.readline()
        file1.close(); out1.close()
        
    #MODULE 2: ANNOTATE THE RESIN GLYCOSIDE LIKE PEAKS
    print ("----------------------------------------------------------------------------------------")
    print ("Annotating peaks Step 1...")
    suffixes=['mz132','mz146','mz162']
    motiffile='3_predictMotifs/motifs.groups'
    dfrags=10
    rgid_step1.full_script([filelist,motiffile,dfrags,suffixes])

    print ("Annotating peaks Step 2...")
    nfile=filelist+".nid1"
    rgid_step2.full_script([nfile,motiffile])
    

print ("All done!")
    
    
