import sys
import os.path
import os
import time
import csv
import melting
from tkinter import filedialog
from tkinter import *
MIN_PROBE_TM = 64
MAX_PROBE_TM = 66
MIN_OVERLAP_TM = 44
MAX_OVERLAP_TM = 60
DEFAULT_PROBE_LENGTH = 15
DEFAULT_OVERLAP_LENGTH = 10
FIRST_OFFSET = 10

class tiler(object):
    
    def __init__(self):
        root = Tk()
        root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
        with open(root.filename, 'r') as f:
            reader = csv.reader(f)
            data = list(reader)
        
        self.results = [[0]*14 for x in range(len(data))]
        self.results[0] = ['Target Name', 'Probe Overlap Dimer','Tm', 'Sense_name', 'Sense', 'Tm', 
        'Length','Mutation Position', 'Antisense_name', 'Antisense', 'Tm', 'Length', 
        'Mutation Position', 'Reverse Complement']
        
        for i in range(1, len(data)):
            try:
                self.results[i][0] = data[i][0]
                self.results[i][3] = data[i][0] + "_sense"
                self.results[i][8] = data[i][0] + "_antisense"
                hd, Tm = self.heterodimer(data[i][2],data[i][1])
                self.results[i][1] = hd
                self.results[i][2] = int(Tm)
                
                sense, Tm, pos = self.extend_sense(hd, data[i][2],data[i][1])
                self.results[i][4] = sense
                self.results[i][5] = int(Tm)
                self.results[i][6] = len(sense)
                self.results[i][7] = pos  
                #self.results[i][14] = sense[:-3]
                
                antisense, Tm = self.extend_antisense(hd, data[i][2],data[i][1])
                copy = antisense
                self.results[i][13] = copy
                antisense = list(antisense)
                antisense.reverse()
                antisense = ''.join(antisense)
                self.results[i][9] = self.convert(antisense)
                self.results[i][10] = int(Tm)
                self.results[i][11] = len(antisense)
                self.results[i][12] = len(antisense) - 10 
                #self.results[i][15] = self.convert(antisense[:-3])
                
            except IndexError:
                i = i-1
                continue
                
                   
            print ("Heterodimer: %s " % hd)
            print ("Sense: %s " % sense)
            print ("Antisense: %s " % antisense)

        
        #Results: [0 Target_name, 1 Probe Overlap Dimer, 2 Tm, 3 Target_name_sense, 4 Sense_tail, 5 Tm, 6 Length, 
        #			7 Mut_position, 8 Target_name_antisense, 9 Antisense_tail, 10 Tm, 11 Length, 12 Mutation Position
        #           13 Reverse Complement, 14 Sense_no_tail, 15 Antisense_no_tail]
        name = root.filename[:-4] + "_results_tail.csv"
        f = open(name, "w")
        writer = csv.writer(f)

        for i in range(len(self.results)):
            writer.writerow(self.results[i])

        self.results2 = self.results.copy()
        for i in range(1, len(self.results)):
        	self.results2[i][4] = self.results2[i][4][:-3]
        	self.results2[i][6] = int(self.results2[i][6]) - 3
        	self.results2[i][9] = self.results2[i][9][:-3]
        	self.results2[i][11] = int(self.results2[i][11]) - 3


        name = root.filename[:-4] + "_results_notail.csv"
        f = open(name, "w")
        writer = csv.writer(f)

        for i in range(len(self.results2)):
        	writer.writerow(self.results2[i])

        
    ##Attaches bases to antisense probe until TM is satisfactory
    def extend_antisense(self, hd, seq, start):
        antisense = []
        start = int(start)
        end = start - 8
        for element in hd:
            antisense.append(element)
        probe_lengthener = start-8 + len(hd)
        antisense = ''.join(antisense)
        Tm = self.analyze(antisense)
            
        while (Tm < MIN_PROBE_TM):
            antisense = list(antisense)
            antisense.append(seq[probe_lengthener])
            antisense = ''.join(antisense)
            probe_lengthener = probe_lengthener + 1
            
            Tm = self.analyze(antisense)
            #print Tm
        
        end = end - 1
        for i in range(3):
        	antisense = list(antisense)
        	if seq[end] != 'T':
        		antisense.insert(0, 'T')
        	else:
        		antisense.insert(0, 'G')
        	end = end - 1
        	antisense = ''.join(antisense)

        return (antisense, Tm)
    
    #Attaches bases to sense probe until TM is satisfactory
    def extend_sense(self, hd, seq, start):
        count = 8
        sense = []
        start = int(start)
        for element in hd:
            sense.append(element)
        probe_lengthener = start-9
        end  = probe_lengthener + len(hd)
        sense = ''.join(sense)
        Tm = self.analyze(sense)
            
        while (Tm < MIN_PROBE_TM):
            sense = list(sense)
            sense.insert(0, seq[probe_lengthener])
            sense = ''.join(sense)
            probe_lengthener = probe_lengthener - 1
            count = count + 1
            
            Tm = self.analyze(sense)
            
        end = end + 1
        for i in range(3):
        	sense = list(sense)
        	if seq[end] != 'A':
        		sense.append('A')
        	else:
        		sense.append('C')
        	end = end + 1
        	sense = ''.join(sense)
        ## Tm should never be > MAX_PROBE_TM because it starts with hd's tm, which is significantly lower, no second while required
        return (sense, Tm, count)
    
    
    def heterodimer(self, seq, start):
        probe_lengthener = 0
        start = int(start)
        hd = seq[start-8:start+7]
            
        Tm = self.analyze(hd)
        
        ## Use probe_lengthener to update index when extending/shortening probes
        ## Extend to the left until Tm is sufficiently high
        while (Tm < MIN_OVERLAP_TM):
            probe_lengthener = probe_lengthener + 1
            hd = seq[start-8:(start+7+probe_lengthener)]
            
            Tm = self.analyze(hd)
            
                
        while (Tm > MAX_OVERLAP_TM):
            probe_lengthener = probe_lengthener - 1
            hd = hd[:-1]
            
            
            Tm = self.analyze(hd)
            
        #print("Hetero Done")
        return (hd, Tm)
    
    ##Simple method for rewriting sense probes into antisense
    def convert(self, probe):
        probe = list(probe)
        for i in range(len(probe)):
            if (probe[i] == 'A' or probe[i] == 'a'):
                probe[i] = 'T'
            elif (probe[i] == 'T' or probe[i] == 't'):
                probe[i] = 'A'
            elif (probe[i] == 'G' or probe[i] == 'g'):
                probe[i] = 'C'
            elif (probe[i] == 'C' or probe[i] == 'c'):
                probe[i] = 'G'
        probe = ''.join(probe)
        return probe
        
    ##Returns TM of probe by using analyzer in IDTDNA.com
    def analyze(self, probe):
        #print(probe)
        #print(str(probe))
        tm = melting.temp(str(probe), DNA_c=250,Na_c=50, Mg_c=5, dNTPs_c=0)
        #print(tm)
        return tm
        
        
## Initializer with "python tiler.py [filename]"
if __name__=='__main__':
    my_tiler = tiler()
