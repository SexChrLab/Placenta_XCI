# In this script, we are looking for samples with rnaseq data for 2 hearts (Heart_Left_Ventricle and Heart_Atrial_Appendage)
import json

heart_left_ventricle_sampleid = set()
heart_left_ventricle = {}

heart_atrial_appendage_sampleid = set()
heart_atrial_appendage = {}

with open("/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json") as json_file:
    data = json.load(json_file)
    for i in data["Heart_Left_Ventricle"]:
        items = i.split("-")
        subjectid = items[0] + "-" + items[1]
        heart_left_ventricle_sampleid.add(subjectid)
        heart_left_ventricle[subjectid] = i
    for i in data["Heart_Atrial_Appendage"]:
        items = i.split("-")
        subjectid = items[0] + "-" + items[1]
        heart_atrial_appendage_sampleid.add(subjectid)
        heart_atrial_appendage[subjectid] = i

shared_sampleid = heart_left_ventricle_sampleid.intersection(heart_atrial_appendage_sampleid)
# print(len(heart_left_ventricle_sampleid))
# print(len(heart_atrial_appendage_sampleid))
# print(len(shared_sampleid))

# Save as a json file
data = {}
data["subjectids"] = []
data["shared_heart_left_ventrible"] = []
data["shared_heart_atrial_appendage"] = []
for i in shared_sampleid:
    data["subjectids"].append(i)
    data["shared_heart_left_ventrible"].append(heart_left_ventricle[i])
    data["shared_heart_atrial_appendage"].append(heart_atrial_appendage[i])

with open("/scratch/tphung3/Placenta_XCI/heart/subjects_with_two_hearts_config.json", "w") as outfile:
    json.dump(data, outfile)

outfile = open("/scratch/tphung3/Placenta_XCI/heart/subjects_with_two_hearts.txt", "w")
header = ["subjectid", "heart_left_ventricle", "heart_atrial_appendage"]
print("\t".join(header), file=outfile)
for i in shared_sampleid:
    out = [i, heart_left_ventricle[i], heart_atrial_appendage[i]]
    print("\t".join(out), file=outfile)

