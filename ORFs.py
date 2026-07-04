
genetic_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "Stop",
    "UAG": "Stop",
    "UGU": "C",
    "UGC": "C",
    "UGA": "Stop",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

def read_fasta(filename):
    names = []
    seqs = []
    current = []

    with open(filename) as file:
        for line in file:
            line = line.strip()
            if line[0] == ">":
                names.append(line[1:])
                if current:
                    seqs.append("".join(current))
                    current = []
            else:
                current.append(line)

    seqs.append("".join(current))
    return names, seqs
ORFs=[]
ORF=''
names, seqs = read_fasta("rosalind_orf.txt")
code=seqs[0]
code=code.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
code1=code.replace('T', 'U')
code2=seqs[0].replace('T', 'U')
codes = [code1, code2]
for m in codes:
    for i in range(len(m)-2):
        if m[i:i+3] == "AUG":
            j=i
            while j+3 <= len(m) and genetic_code[m[j:j+3]] != "Stop":
                ORF+=genetic_code[m[j:j+3]]
                j+=3
            if j+3 <= len(m) and genetic_code[m[j:j + 3]] == "Stop":
                ORF += "Stop"
            ORFs.append(ORF)
            ORF=''
for k in range(len(ORFs)):
    if "Stop" in ORFs[k] and ORFs.index(ORFs[k])==k:
        print(ORFs[k][0:-4])

