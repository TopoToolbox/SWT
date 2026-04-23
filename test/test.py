import sys
import csv

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    rval = 0
    
    with open(file1, newline='') as f1:
        with open(file2, newline='') as f2:
            r1 = csv.reader(f1, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
            r2 = csv.reader(f2, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
            try:
                for (l1, l2) in zip(r1, r2, strict=True):
                    for (x, y) in zip(l1, l2, strict=True):
                        if abs(x - y) > abs(x) * 0.0:
                            print("<", l1)
                            print('---')
                            print(">", l2)
                            rval = 1
            except:
                pass
        
    return rval
    
if __name__ == '__main__':
    sys.exit(main())
