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


def overlap(a, b):
    min_len = min(len(a), len(b))
    limit = min_len // 2

    for k in range(min_len, limit, -1):
        if a[-k:] == b[:k]:
            return k

    return 0


names, reads = read_fasta("rosalind_long.txt")

next_read = {}
prev_read = {}
overlap_len = {}

for i in range(len(reads)):
    for j in range(len(reads)):
        if i != j:
            k = overlap(reads[i], reads[j])
            if k > 0:
                next_read[i] = j
                prev_read[j] = i
                overlap_len[i] = k
start = None
for i in range(len(reads)):
    if i not in prev_read:
        start = i
        break

result = reads[start]
current = start

while current in next_read:
    nxt = next_read[current]
    k = overlap_len[current]
    result += reads[nxt][k:]
    current = nxt

print(result)






