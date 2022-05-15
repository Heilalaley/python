import sys

def cheak(x0, y0, x1, y1, R):
    if (R>( (x1-x0)**2 + (y1-y0)**2 )**(1/2)):
        return 1
    if (R==( (x1-x0)**2 + (y1-y0)**2 )**(1/2)):
        return 0
    return 2

def read_files(inputfile1, inputfile2):
    result = []
    f1 = open(inputfile1)
    Cy = f1.readline().split()
    R = float(f1.readline())
    f2 = open(inputfile2)
    for line in f2:
        Po = line.split()
        result.append(cheak(float(Cy[0]), float(Cy[1]), float(Po[0]), float(Po[1]), R))
    f1.close()
    f2.close()
    print(*result, sep='\n')
    return None
        
def main():
    try:
        f1 = open(sys.argv[1])
        f2 = open(sys.argv[2])
        f1.close()
        f2.close()
    except FileNotFoundError:
        print("Files not found Error")
        return None
    read_files(sys.argv[1], sys.argv[2])
    return None

main()