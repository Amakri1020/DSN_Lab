import sys
import csv

out_name = sys.argv[1]
f = open("results.csv")
csv_f = csv.reader(f)
output = open("%s.csv"%out_name, 'w')
writer = csv.writer(output)
writer.writerow(["Heterodimer", "Tm", "Sense", "Tm", "Antisense", "Tm"])
flag = 0
heteros = []
senses = []
antisenses = []

for row in csv_f:
    res = ""
    for element in row[0]:
        if (element == 'C' or element == 'T' or element == 'G' or element == 'A'):
            res = res + str(element)
    #writer.writerow([res, row[1]])
    if flag == 0:
        if row[1] == 'k':
            flag = 1
            continue
        else:
            heteros.append([res, row[1]])
    if flag == 1:
        if row[1] == 'k':
            flag = 2
            continue
        else:
            senses.append([res, row[1]])
    if flag == 2:
        antisenses.append([res, row[1]])
        
        
for i in range(len(heteros)):
    writer.writerow([heteros[i][0], heteros[i][1], senses[i][0], senses[i][1], antisenses[i][0], antisenses[i][1]])
    
    
    
    
    
    
    
    
    
    
    
    
    
