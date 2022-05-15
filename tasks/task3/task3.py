import json
import sys

def perebor(tests, values):
    for x in tests['tests']:
        if_dict(x, values)
    return tests

def if_dict(x, values):
    x['value'] = find(x['id'], values)
    for y in x.values():
        if (type(y) == list):
            if_list(y, values)
    return None

def if_list(x, values):
    for y in x:
        if (type(y) == dict):
            if_dict(y, values)
    return None

def find(x,  values):
    for y in values['values']:
        if (x == y['id']):
            return y['value']
    print ("Id ",x," not found in values")
    return ""


def main():
    with open(sys.argv[1]) as file1:
        tests = json.load(file1)
    with open(sys.argv[2]) as file2:    
        values = json.load(file2)
    report = perebor(tests, values)
    with open("report.json","w") as write_file: 
        json.dump(report ,write_file, indent=2) 
    return 0

main()