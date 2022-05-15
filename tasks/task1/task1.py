import sys

def path(n, m):
    if (m>n):
        m = m % n
    a = 1
    c = [1]
    while ((a != 1) or (len(c) ==1)):
        a = a + m-1
        if (a>n):
            a = a%n
        c.append(a)
    c.pop()
    print(*c, sep="")
    return 0

def main():
    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
    except ValueError:
            print("argument type error")
            return None
    path(n, m)
    return 0

main()
