import sys
import os.path
import os
import time
import csv
import melting
MIN_PROBE_TM = 64
MAX_PROBE_TM = 66
MIN_OVERLAP_TM = 44
MAX_OVERLAP_TM = 49
DEFAULT_PROBE_LENGTH = 15
DEFAULT_OVERLAP_LENGTH = 10
FIRST_OFFSET = 10

class tiler(object):
    
    def __init__(self, arg1):
        with open(arg1, 'r') as f:
            reader = csv.reader(f)
            data = list(reader)
        
        self.results = [[0]*6 for x in range(len(data))]
        self.results[0] = ['Heterodimer','Tm', 'Sense', 'Tm', 'Antisense', 'Tm']
        
        for i in range(1, len(data)):
            try:
                hd, Tm = self.heterodimer(data[i][2],data[i][1])
                self.results[i][0] = hd
                self.results[i][1] = int(Tm)
                
                sense, Tm = self.extend_sense(hd, data[i][2],data[i][1])
                self.results[i][2] = sense
                self.results[i][3] = int(Tm)   
                
                antisense, Tm = self.extend_antisense(hd, data[i][2],data[i][1])
                antisense = list(antisense)
                antisense.reverse()
                antisense = ''.join(antisense)
                self.results[i][4] = self.convert(antisense)
                self.results[i][5] = int(Tm)
                
            except IndexError:
                i = i-1
                continue
                
                   
            print ("Heterodimer: %s " % hd)
            print ("Sense: %s " % sense)
            print ("Antisense: %s " % antisense)

        
                
        f = open("results.csv", "w")
        writer = csv.writer(f)

        for i in range(len(self.results)):
            writer.writerow(self.results[i])
        
    ##Attaches bases to antisense probe until TM is satisfactory
    def extend_antisense(self, hd, seq, start):
        antisense = []
        start = int(start)
        for element in hd:
            antisense.append(element)
        probe_lengthener = start-3 + len(hd)
        antisense = ''.join(antisense)
        Tm = self.analyze(antisense)
            
        while (Tm < MIN_PROBE_TM):
            antisense = list(antisense)
            antisense.append(seq[probe_lengthener])
            antisense = ''.join(antisense)
            probe_lengthener = probe_lengthener + 1
            
            Tm = self.analyze(antisense)
            #print Tm
                
        return (antisense, Tm)
    
    #Attaches bases to sense probe until TM is satisfactory
    def extend_sense(self, hd, seq, start):
        sense = []
        start = int(start)
        for element in hd:
            sense.append(element)
        probe_lengthener = start-4
        sense = ''.join(sense)
        Tm = self.analyze(sense)
            
        while (Tm < MIN_PROBE_TM):
            sense = list(sense)
            sense.insert(0, seq[probe_lengthener])
            sense = ''.join(sense)
            probe_lengthener = probe_lengthener - 1
            
            Tm = self.analyze(sense)
            
        ## Tm should never be > MAX_PROBE_TM because it starts with hd's tm, which is significantly lower, no second while required
        return (sense, Tm)
    
    
    def heterodimer(self, seq, start):
        probe_lengthener = 0
        start = int(start)
        hd = seq[start-3:start+3]
            
        Tm = self.analyze(hd)
        
        ## Use probe_lengthener to update index when extending/shortening probes
        ## Extend to the left until Tm is sufficiently high
        while (Tm < MIN_OVERLAP_TM):
            probe_lengthener = probe_lengthener + 1
            hd = seq[start-3:(start+3+probe_lengthener)]
            
            Tm = self.analyze(hd)
            
                
        while (Tm > MAX_OVERLAP_TM):
            probe_lengthener = probe_lengthener - 1
            hd = hd[:-1]
            
            
            Tm = self.analyze(hd)
            
        print("Hetero Done")
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
        print(probe)
        print(str(probe))
        tm = melting.temp(str(probe), DNA_c=250,Na_c=50, Mg_c=5, dNTPs_c=0)
        print(tm)
        return tm
        
        
## Initializer with "python tiler.py [filename]"
if __name__=='__main__':
    my_tiler = tiler(sys.argv[1])
