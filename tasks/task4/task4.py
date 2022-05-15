import sys
import math

def function(input):
    f1 = open(input)
    arr = []
    for line in f1:
        arr.append(int(line))
    f1.close()
    arr.sort()
    if (len(arr) % 2 == 1):
        med = arr[(math.floor(len(arr)/2) + 1)]
    else:
        med = (arr[len(arr)//2] + arr[len(arr)//2 + 1])/2
        med = math.floor(med)
    resault = 0
    for i in arr:
        resault = abs(med - i) + resault
    print(resault)
    return 0


def main():
    try:
        f1 = open(sys.argv[1])
        f1.close()
    except FileNotFoundError:
        print("Files not found Error")
        return 1
    function(sys.argv[1])
    return 0

main()