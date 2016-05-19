import re
import shutil

def main(file_name,pattern):
    f = open(file_name,'r')
    fix_file_name = file_name + ".fix"
    new = open(fix_file_name,"w")
    lines = f.readlines()
    for l in lines:
        if re.search(pattern,l):
            r = re.match('([^"]*")(.*)("[^"]*)',l)
            inner = r.group(2)
            inner = re.sub(r'([^\\])"',r'\1\"',inner)
            l = r.group(1) + inner + r.group(3)
        new.write(l)
    new.close()
    f.close()
    shutil.move(fix_file_name,file_name)
            
if __name__ == "__main__":
    import sys
    main(sys.argv[1],"INTERFACE_COMPILE_DEFINITIONS")
