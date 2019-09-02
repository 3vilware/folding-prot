pdb = open("1vbk.pdb", "r")
for line in pdb:
    if line[:4] == 'ATOM': #and line[12:16] == " CA ":
        print (line)
        x = input('...')