# README:
# In this script, I want to return a list of gtex individuals from the file `participant.tsv` because I need this list in order to download the count files from google cloud
# If you need to rerun this script, replace the path on line 7.

individuals_list = []

with open("/scratch/tphung3/Placenta_XCI/gtex/01_download_data/participant.tsv", 'r') as f:
    for line in f:
        if not line.startswith('entity'):
            individuals_list.append(line.rstrip('\n').split('\t')[0])

print (len(individuals_list))
print (individuals_list)
