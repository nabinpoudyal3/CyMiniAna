## Merge json files
import json

print " Load JSON data "
rsg400 = json.load(open("data/rsg_400_77_27Jan.json"))
rsg500 = json.load(open("data/rsg_500_77_27Jan.json"))
rsg750 = json.load(open("data/rsg_750_77_27Jan.json"))
rsg1000 = json.load(open("data/rsg_1000_77_27Jan.json"))
rsg2000 = json.load(open("data/rsg_2000_77_27Jan.json"))
rsg3000 = json.load(open("data/rsg_3000_77_27Jan.json"))

print " Merge files "
data = {}
for key in rsg400.keys():
    if key=='metadata':
        data[key] = rsg400[key]
    else:
        data[key] = rsg400[key]+rsg500[key]+rsg750[key]+rsg1000[key]+rsg2000[key]+rsg3000[key]

print " Dump JSON data "
with open('data/rsg_inclu_77_27Jan.json', 'w') as outfile:
    json.dump(data, outfile)
