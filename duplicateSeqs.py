names = []
with open('/home/labs/fleishman/jonathaw/t_coffee_aln/all_docs_23.12/all_dockerins_Sept2014b.one', 'r') as f:
    for line in f:
        if line[0] == '>':
            names.append(line)

for i in range(0, len(names)):
    for j in range(i+1, len(names)):
        if str(names[i]) == str(names[j]):
            print names[i], 'equals', names[j]