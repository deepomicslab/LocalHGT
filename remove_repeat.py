"""
some breakpoint pairs are too close, we only keep a 
representive one from the nearby bkps. 
"""

import sys

infile = sys.argv[1]
outfile = sys.argv[2]
cutoff = 50

def main():
    inf = open(infile)
    outf = open(outfile, 'w')
    record = []
    for line in inf:
        flag = True
        array = line.split(",")
        if array[0] != "from_ref":
            array[1] = int(array[1])
            array[3] = int(array[3])       
            for rec in record:
                if array[0] == rec[0] and abs(array[1] - rec[1])<cutoff \
                and array[2] == rec[2] and abs(array[3] - rec[3])<cutoff:
                    flag = False
                    break
                elif array[2] == rec[0] and abs(array[3] - rec[1])<cutoff \
                and array[0] == rec[2] and abs(array[1] - rec[3])<cutoff:
                    flag = False
                    break
        if flag:
            record.append(array[:5])
            print (line, end='', file = outf)

    inf.close()
    outf.close()

main()