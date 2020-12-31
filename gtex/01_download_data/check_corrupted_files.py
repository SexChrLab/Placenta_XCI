# README
# The files on line 14 are now loaded on shared google team drive at: wilsonlab/gtex/version8/chrX_WASP_raw_counts
# After downloading the data with the script `download_asereadcoutner_counts.sh`, there are some files with content like this: No such object: fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_WASP_chrX_raw_counts_by_subject/GTEX-1J8EW.v8.readcounts.chrX.txt.gz. We want to remove these individuals from further analyses

# In this script, I want to subset the ASEReadCounter counts downloaded from GTEX for each rna_id
import pandas as pd
import os

if os.path.exists('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/failed_files.txt'):
    outfile = open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/failed_files.txt', 'a')
else:
    outfile = open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/failed_files.txt', 'w')

for file in os.listdir('/home/tphung3/mwilsons/public_data/controlled_access/gtex/version8/chrX_WASP_raw_counts/'):
    try:
        data = pd.read_csv(os.path.join('/home/tphung3/mwilsons/public_data/controlled_access/gtex/version8/chrX_WASP_raw_counts/', file), sep='\t',
                           compression='gzip')
    except:
        j = file.split('.')
        id = j[0]
        print (id, file=outfile)
        continue




