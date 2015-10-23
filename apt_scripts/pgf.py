dat = {}
with open('HuGene-1_0-st-v1.r4.pgf') as f:
  for line in f:
    if "#" in line: continue
    line = line.strip().split('\t')
    if len(line) == 1: continue
    elif len(line) == 2:
      probeset = line[0]
    else:
      probe = line[0]
      if probe in dat: 
        dat[probe].append(probeset)
      else: 
        dat[probe] = [ probeset ]

with open('out.txt','w') as out_f:
  out_f.write('probe_id\tprobeset_id\n')
  with open('old_snps.txt') as f:
    for line in f:
      line=line.strip()
      probesets = dat[line]
      for i in probesets:
        out_f.write('%s\t%s\n' % (line, i))

