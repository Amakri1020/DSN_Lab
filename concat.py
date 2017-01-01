import sys
import csv

f = open("results.csv")
csv_f = csv.reader(f)
output = open("results_c.csv", 'wb')
writer = csv.writer(output)
writer.writerow(["Heterodimer", "Tm"])
flag = 0

for row in csv_f:
    res = ""
    for element in row[0]:
        if (element == 'C' or element == 'T' or element == 'G' or element == 'A'):
            res = res + str(element)
    print res
    if row[1] == 'k':
        if flag == 0:
            row = ["Sense", "Tm"]
            flag = 1
            writer.writerow(row)
            continue
        else:
            row = ["Antisense", "Tm"]
            writer.writerow(row)
            continue
    writer.writerow([res, row[1]])
