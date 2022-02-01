
import os
import subprocess
import pandas as pd
import argparse
import sys



print("\n\nif conda environment /data/anaconda3/envs/porecov not activated, plase activate......\n\n")

parser=argparse.ArgumentParser(description="It is for QC and transfering covid data(Passed fastq) to LMO Server")
parser.add_argument('-FastqPassFolder','--FastqPassFolder', help="Full path of fastq_pass should end with /", required=True)
parser.add_argument('-SampleSheet','--SampleSheet', help="Sample sheet; should have two columns namely Sample_ID and ONT_Barcode_ID \
                    file should be separated by coma", required=True)

parser.add_argument('-OutputName','--OutputName', help="Provide Set name and number of Samples, like SET46_93Samples", required=True)

parser.add_argument('-ServerName','--ServerName', help="ServerName with path; like host@10.12.13.14:/Data/Data/", required=True)


args=parser.parse_args()
FastqPassFolder=args.FastqPassFolder
SampleSheet=args.SampleSheet
OutputName=args.OutputName
ServerName=args.ServerName #Server name and location to transfer



#FastqPassFolder="/data/NGC_Data/Nanopore_Output/2022/Covid/COVID_SET_56_89_samples_190122/no_sample/20220119_1209_X1_FAR39728_4584f146/fastq_pass/"
#SampleSheet="SampleSheet.csv"
#OutputName="SET_56_89_samples"
#ServerName="ngc@10.117.173.91:/Data/NGC_Data/" #Server name and location to transfer



CovidNanopore_1200Bed=" /data/Software_Resources/covid_refseq/Porecov/CovidNanopore_1200.bed" #Bed file

CombinedFastqFolder="CombinedFastq"
NTC_PTC_fOLDER="NTC_PTC"
Date = os.popen("date +%d-%b-%Y").read().split("\n")[0]
FolderName="Nanopore_Data_"+Date+"_"+OutputName

os.system("mkdir -p "+ FolderName)
os.system("mkdir -p "+ NTC_PTC_fOLDER)
os.system("mkdir -p "+ CombinedFastqFolder)
Sampledf=pd.read_csv(SampleSheet)

Sampledf['Sample_ID']=Sampledf['Sample_ID'].str.strip()
Sampledf['ONT_Barcode_ID']=Sampledf['ONT_Barcode_ID'].str.strip()




for Sample in Sampledf['Sample_ID']:
    if Sample =="NTC":
        SampleName=str(Sampledf[Sampledf['Sample_ID']==Sample]['ONT_Barcode_ID']).split()[1]
        print(SampleName)
        os.system("zcat "+ FastqPassFolder + SampleName +"/*fastq.gz > " + CombinedFastqFolder +"/"+ SampleName +".fastq" )
        os.system("cp -r "+ FastqPassFolder + SampleName + " "+ NTC_PTC_fOLDER  )
    elif Sample =="PTC":
        SampleName=str(Sampledf[Sampledf['Sample_ID']==Sample]['ONT_Barcode_ID']).split()[1]
        print(SampleName)
        os.system("zcat "+ FastqPassFolder + SampleName +"/*fastq.gz > " + CombinedFastqFolder  +"/" + SampleName +".fastq" )
        os.system("cp -r "+ FastqPassFolder + SampleName + " "+ NTC_PTC_fOLDER  )
    else:
        SampleName=str(Sampledf[Sampledf['Sample_ID']==Sample]['ONT_Barcode_ID']).split()[1]
        print(SampleName)
        os.system("zcat "+ FastqPassFolder + SampleName +"/*fastq.gz > " + CombinedFastqFolder  +"/" + SampleName +".fastq" )
        os.system("cp -r "+ FastqPassFolder + SampleName + " "+ FolderName  )



#Total Base, Reads and Mb data
def BaseCount():
        cmd= "ls " + CombinedFastqFolder + "/*.fastq >fasqfile.txt && rm Sample_Names.txt temp.txt"
        os.system(cmd)
        os.system("cat fasqfile.txt")
        os.system("while read Fastq; do echo $Fastq  | cut -d '.' -f 1 >> temp.txt ;done < fasqfile.txt")
        os.system("cat temp.txt |sort| uniq |sed '/^$/d' > Sample_Names.txt")
        os.system("cat -n Sample_Names.txt")
        cmd='''
        while read Samplename ; do
        Name=${Samplename}
        Totalreads=$(cat ${Samplename}.fastq | egrep "@" -A 1| egrep -v "@|--" | wc -l)
        TotalBases=$(cat ${Samplename}.fastq | egrep "@" -A 1| egrep -v "@|--" | wc -m)
        TotalDataInMb=$((TotalBases / 1000000))
        echo  "${Samplename} ${Totalreads} ${TotalBases} ${TotalDataInMb}" >> PerSample_Data_Details.txt; done < Sample_Names.txt   '''
        os.system(cmd)

try:
    os.system("rm PerSample_Data_Details.txt")
except:
    pass

BaseCount()

os.system("sed  -i '1i Samplename Totalreads TotalBases TotalDataInMb' PerSample_Data_Details.txt")
os.system("sed -i 's|CombinedFastq/||g' PerSample_Data_Details.txt")
os.system("cat PerSample_Data_Details.txt")

SamleCoveDf=pd.read_csv("PerSample_Data_Details.txt",sep="\s",engine='python')

##----------------------------------------------------------------------------------------------------------
try:
    Sampledf.drop("Totalreads",axis=1,inplace=True)
except:
    pass

try:
    Sampledf.drop("TotalBases",axis=1,inplace=True)
except:
    pass


try:
    Sampledf.drop("TotalDataInMb",axis=1,inplace=True)
except:
    pass

try:
    Sampledf.drop("MeanDepth",axis=1,inplace=True)
except:
    pass

#-------------------------------------------------------------------------------------------------------------



SamleCoveDf=pd.merge(Sampledf,SamleCoveDf,left_on="ONT_Barcode_ID",right_on="Samplename",how="outer").drop("Samplename",axis=1)
SamleCoveDf['MeanDepth']=SamleCoveDf['TotalBases']/29903
SamleCoveDf.to_csv(Date+"_"+OutputName+"_PerSample_Data_Details.csv",index=None)
os.system("rm PerSample_Data_Details.txt")
os.system("cat " + Date+"_"+OutputName+"_PerSample_Data_Details.csv")

#Data Transfer to LMO
print("Amount of data in each samples to be send to LMO\n")
os.system("du -hs " +FolderName+"/*" )
#os.system("tree -d --du -hs " +FolderName+"/" )


LMO=input("Needs to transfere data to LMO(Print Yes or No: ")
LMO=LMO.upper()

##Checking samples present in the sample sheet and Passed fastq folders are same or not

Sampledf['Sample_ID']=Sampledf['Sample_ID'].str.strip()
Sampledf['ONT_Barcode_ID']=Sampledf['ONT_Barcode_ID'].str.strip()



LMO_Barcodes=os.popen('ls '+FolderName).read().split()
LMO_Barcodes = [x.strip(' ') for x in LMO_Barcodes]
LMO_Barcodes.sort()

SampleBarcode=Sampledf[((Sampledf['Sample_ID']!='NTC') & (Sampledf['Sample_ID']!='PTC'))]['ONT_Barcode_ID'].to_list()
SampleBarcode = [x.strip(' ') for x in SampleBarcode]
SampleBarcode.sort()



if LMO_Barcodes != SampleBarcode:
    print ("Number of samples in the sample sheet not mathing with Barcode")
    print(LMO_Barcodes)
    print(SampleBarcode)
    sys.exit("Please correct")
else :
    pass



if LMO=="YES":
    os.system("scp -r "+ FolderName + " " +ServerName)
else:
    print("Data Not transfering to LMO server\n\n\n")

#FastQC
def fastqc():
    cmd="mkdir -p QC/FastQC/Raw_Fastq_QC"
    os.system(cmd)
    
    cmd="for Samplename in $(cat Sample_Names.txt); do echo \"fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/ ${Samplename}.fastq\" ; done | parallel -j2"
    os.system(cmd)

fastqc()


#MultiQC
os.system("rm -rf QC/MultiQC/FastQC/RawFastq_QC/; mkdir -p QC/MultiQC/FastQC/RawFastq_QC/")
os.system("multiqc -o QC/MultiQC/FastQC/RawFastq_QC QC/FastQC/Raw_Fastq_QC/")
os.system("rm temp.txt fasqfile.txt Sample_Names.txt")



#Run PORECOV

#Run PORECOV

Analysis=input("Needs to Complete Covid analysis or not Print Yes or No: ")

Analysis=Analysis.upper()

print(Sampledf.head())


if Analysis=="YES" :
        samplesheetdf=Sampledf[['Sample_ID','ONT_Barcode_ID']]
        samplesheetdf=samplesheetdf.rename(columns={"Sample_ID":"_id","ONT_Barcode_ID":"Status"})
        samplesheetdf.loc[:,"Description"]="good"
        samplesheetdf=samplesheetdf[['_id','Status','Description']]
        print(samplesheetdf)
        samplesheetdf.to_csv("samplesheet.csv",index=None)
        cmd="cp -r NTC_PTC/* "+FolderName
        os.system(cmd)
        cmd="nextflow run replikation/poreCov  --fastq_pass " + FolderName + " \
            --cores 2 --samples samplesheet.csv --primerV V1200 --output poreCov_results \
            -profile local,docker -r 0.11.0  --rapid --update \
            --databases /data/Software_Resources/covid_refseq/Porecov/databases/ --cachedir /data/Software_Resources/covid_refseq/Porecov/singularity_images/ "
        os.system(cmd)
else:
        print("Analysis not performing")



if Analysis=="YES" :
        for sample in Sampledf['Sample_ID']:
                cmd='mkdir -p QC/Mosedepth/'
                os.system(cmd)
                cmd="mosdepth --by "+ CovidNanopore_1200Bed + " QC/Mosedepth/"+sample +" poreCov_results/2.Genomes/"+sample+"/"+sample+'''_mapped_MN908947.3.primertrimmed.sorted.bam \
        --thresholds 1,10,20,30,50,100,200,300,400,500,1000,1500,2000,5000 --quantize 0:1:5:10:50:100:150: --fast-mode --flag 772'''
                os.system(cmd)
        os.system("python /data/Software_Resources/NGC_scripts/MosedepthTargetComparison.py --MosedepthOutPutFolder QC/Mosedepth/ --Outputfolder QC/Mosedepth/")



os.system("rm -rf work")
os.system("rm -rf NTC_PTC CombinedFastq " +FolderName)