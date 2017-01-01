from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
import sys
import os.path
import os
import time
import csv
MIN_PROBE_TM = 64
MAX_PROBE_TM = 66
MIN_OVERLAP_TM = 44
MAX_OVERLAP_TM = 49
DEFAULT_PROBE_LENGTH = 15
DEFAULT_OVERLAP_LENGTH = 10
FIRST_OFFSET = 10

class tiler(object):
    
    def __init__(self, arg1):
        with open(arg1, 'rb') as f:
            reader = csv.reader(f)
            data = list(reader)
        ## str contains a list of ALL bases irrespective of start/stop entries
        
        ## open site, make settings for TM calculation 
        self.driver = webdriver.Firefox()
        self.driver.get("http://www.idtdna.com/calc/analyzer")
        Mg_conc = self.driver.find_element_by_xpath("//div[@id='OligoAnalyzer']/div[2]/div[1]/div/div[2]/table/tbody/tr[4]/td[2]/input")
        Mg_conc.send_keys("5")
        
        self.heterodimers = []
        self.heteroTM = []
        self.senses = []
        self.senseTM = []
        self.antisenses = []
        self.antiTM = []
        
        for seq in data:
            hd = "rip"
            sense = "ripper"
            antisense = "ripperino"
            try:
                hd, Tm = self.heterodimer(seq[1],151)
                self.heterodimers.append(hd)
                self.heteroTM.append(Tm)
            except IndexError:
                self.heterodimers.append(["fail"])
                self.heteroTM.append(["fail"])
                self.senses.append(["fail"])
                self.senseTM.append(["fail"])
                self.antisenses.append(["fail"])
                self.antiTM.append(["fail"])
                pass
            try:
                sense, Tm = self.extend_sense(hd, seq[1])
                self.senses.append(sense)
                self.senseTM.append(Tm)
            except IndexError:
                self.senses.append(["fail"])
                self.antisenses.append(["fail"])
                self.senseTM.append(["fail"])
                self.antiTM.append(["fail"])
                pass
            try:    
                antisense, Tm = self.extend_antisense(hd, seq[1])
                self.antisenses.append(antisense)
                self.antiTM.append(Tm)
            except IndexError:
                self.antisenses.append(["fail"])
                self.antiTM.append(["fail"])
                pass    
            print ("Heterodimer: %s " % hd)
            print ("Sense: %s " % sense)
            print ("Antisense: %s " % antisense)

        
        
        self.driver.close()
        
        f = open("results.csv", "wb")
        writer = csv.writer(f)
        for i in range(len(self.heterodimers)):
            writer.writerow([self.heterodimers[i],self.heteroTM[i]])
        writer.writerow("skip")
        for i in range(len(self.senses)):
            writer.writerow([self.senses[i], self.senseTM[i]])
        writer.writerow("skip")
        for i in range(len(self.antisenses)):
            writer.writerow([self.antisenses[i], self.antiTM[i]])
        
        
    ##Attaches bases to antisense probe until TM is satisfactory
    def extend_antisense(self, hd, seq):
        antisense = []
        for element in hd:
            antisense.append(element)
        probe_lengthener = 144 + len(hd)
        
        Tm = self.analyze(antisense)
        if (Tm[2] == " "):
            Tm = float(Tm[0] + Tm[1])
        else:
            Tm = float(Tm[0] + Tm[1])  
            
        while (Tm < MIN_PROBE_TM):
            antisense.append(seq[probe_lengthener])
            probe_lengthener = probe_lengthener + 1
            
            Tm = self.analyze(antisense)
            #print Tm
            if (Tm[2] == " "):
                Tm = float(Tm[0] + Tm[1])
            else:
                Tm = float(Tm[0] + Tm[1])
                
        return (antisense, Tm)
    
    #Attaches bases to sense probe until TM is satisfactory
    def extend_sense(self, hd, seq):
        sense = []
        for element in hd:
            sense.append(element)
        probe_lengthener = 143
        
        Tm = self.analyze(sense)
        if (Tm[2] == " "):
            Tm = float(Tm[0] + Tm[1])
        else:
            Tm = float(Tm[0] + Tm[1])
            
        while (Tm < MIN_PROBE_TM):
            sense.insert(0, seq[probe_lengthener])
            probe_lengthener = probe_lengthener - 1
            
            Tm = self.analyze(sense)
            #print Tm
            if (Tm[2] == " "):
                Tm = float(Tm[0] + Tm[1])
            else:
                Tm = float(Tm[0] + Tm[1])
            
        ## Tm should never be > MAX_PROBE_TM because it starts with hd's tm, which is significantly lower, no second while required
        return (sense, Tm)
    
    
    def heterodimer(self, seq, start):
        probe_lengthener = 0
        hd = seq[146:152]
            
        Tm = self.analyze(hd)
        print Tm
        if (Tm[2] == " "):
            Tm = float(Tm[0] + Tm[1])
        elif (Tm[0] == "0"): 
                Tm = 0
        else:
            Tm = float(Tm[0] + Tm[1])
        
        ## Use probe_lengthener to update index when extending/shortening probes
        ## Extend to the left until Tm is sufficiently high
        while (Tm < MIN_OVERLAP_TM):
            probe_lengthener = probe_lengthener + 2
            hd = seq[144:(152+probe_lengthener)]
            
            Tm = self.analyze(hd)
            print Tm
            if (Tm[2] == " "):
                Tm = float(Tm[0] + Tm[1])
            elif (Tm[0] == "0"): 
                Tm = 0
            else:
                print Tm
                Tm = float(Tm[0] + Tm[1])
            
                
        while (Tm > MAX_OVERLAP_TM):
            probe_lengthener = probe_lengthener - 1
            hd = hd[:-1]
            
            
            Tm = self.analyze(hd)
            #In case analysis returns int vs float
            if (Tm[2] == " "):
                Tm = float(Tm[0] + Tm[1])
            else:
                Tm = float(Tm[0] + Tm[1])
            
        
        self.start = start + DEFAULT_OVERLAP_LENGTH + probe_lengthener + 1
        return (hd, Tm)
    
    ##Simple method for rewriting sense probes into antisense
    def convert(self, probe):
        for i in range(len(probe)):
            if (probe[i] == 'A' or probe[i] == 'a'):
                probe[i] = 'T'
            elif (probe[i] == 'T' or probe[i] == 't'):
                probe[i] = 'A'
            elif (probe[i] == 'G' or probe[i] == 'g'):
                probe[i] = 'C'
            elif (probe[i] == 'C' or probe[i] == 'c'):
                probe[i] = 'G'
        return probe
        
    ##Returns TM of probe by using analyzer in IDTDNA.com
    def analyze(self, probe):
        
        self.driver.implicitly_wait(1)
        time.sleep(1)
        probe_string = ''.join(probe)
        self.driver.find_element_by_css_selector('textarea.form-control').clear()
        elem = self.driver.find_element_by_css_selector('textarea.form-control')
        elem.send_keys(probe_string)
        
        #time.sleep(3)
        self.driver.implicitly_wait(3)
        
        self.driver.find_element_by_css_selector("button.btn.btn-primary.btn-md.btn-block").click()
        
        time.sleep(3)
        self.driver.implicitly_wait(3)
        
        result = self.driver.find_element_by_xpath("//div[@id='OAResults']/div/div[1]/div[3]/div/div/table/tbody/tr[5]/td[2]/span")
        result = result.text
        while result == None:
            result = self.driver.find_element_by_xpath("//div[@id='OAResults']/div/div[1]/div[3]/div/div/table/tbody/tr[5]/td[2]/span")
            result = result.text
        return result
        
        
## Initializer with "python tiler.py [filename]"
if __name__=='__main__':
    my_tiler = tiler(sys.argv[1])
