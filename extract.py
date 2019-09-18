from Bio import Entrez, SeqIO
import re

      

Entrez.email = 'vandoorslaer@email.arizona.edu'


def countoverlappingdistinct(pattern, thestring):
  total = 0
  start = 0
  there = re.compile(pattern)
  while True:
    mo = there.search(thestring, start)
    if mo is None: return total
    total += 1
    start = 1 + mo.start()


accessions=[]

with open("accessions.txt","rU") as f:
    for l in f:
        accessions.append(l.strip())


with open("results.csv","w") as out:
    print >>out, ",".join(["accession","protein","CxC count","100 nt upstream"])
    for accession in accessions:
        with open("temp.gb","w") as GB:
            handle=Entrez.efetch(db='nucleotide',id=accession,rettype='gb', retmod="text")
            GB.write(handle.read())
            handle.close()
    
        with open("temp.gb","rU") as f:
            for line in f:
                if "DBSOURCE" in line:
                    x = line.strip().split(" ")[-1]
        
    
        
        with open("temp.gb","w") as GB:
            handle=Entrez.efetch(db='nucleotide',id=x,rettype='gb', retmod="text")
            GB.write(handle.read())
            handle.close()
        for rec in SeqIO.parse("temp.gb","genbank"):
            for feature in rec.features:
                if feature.type == "CDS":
                    if accession in feature.qualifiers["protein_id"]:
                        m = re.search("(\d*):(\d*)",str(feature.location))
                        if m:
                            start, end = int(m.groups()[0]), int(m.groups()[1])
                            upstream = rec.seq[start-100:start]
                        protein = feature.location.extract(rec).seq.translate()
                        print >>out, ",".join(map(str,[accession, protein, countoverlappingdistinct("C.C",str(protein[-10:])), upstream]))