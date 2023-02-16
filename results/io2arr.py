
import numpy as np
def file2twoarr(iname, nol=None):
    num_of_lines = 0
    time = []
    dens = []
    with open(iname, 'r') as f:
        for line in f:
            if not(line[0].isdigit()):
                continue
            num_of_lines += 1
            ll = line.split()
            time.append(float(ll[0]))
            dens.append(float(ll[1]))
            if nol is not None and num_of_lines >= nol:
                break
    return time, dens

# Copy a two-columned array to file
def save_to_file(name, x, y):
    with open(name, 'w') as f:
        for i in range(len(x)):
            f.write('{:.8f}\t{:.8f}\n'.format(x[i], y[i]))
            
    return;

# Divide catalogue <iname> into subcatalogues of time length <time_lng>
# The names are obtained from the input name by adding the number of the catalogue to the end of the filename
def divide(iname, time_lng, oname=None, datacolumn=2, num_file=0):
    #eliminate the extension
    idot = -1
    for i in range(len(iname)):
        if iname[i] == '.':
            idot = i 
    if oname != None:
        name = oname
    elif idot >= 0:
        name = iname[:idot]
    else:
        name = iname
    isfrst_input_line = True
    files = []
    with open(iname, 'r') as f:
        isoutnew = True
        for line in f:
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            time = int(np.round(float(ll[0])))
            if isfrst_input_line:
                time_break = time
                isfrst_input_line = False
            if isoutnew:
                num_file += 1
                oname = f'{name}{num_file}part.txt'
                fo = open(oname, 'w')
                files.append(fo)
                time_break += time_lng
                fo.write('Time\tHeight\n\n')
                isoutnew = False
            if time < time_break:
                ln_split = line.split()
                if len(ln_split) >= datacolumn+1:
                    fo.write(f'{ln_split[0]}\t{ln_split[datacolumn]}\n')
            else:
                fo.close()
                isoutnew = True
    for fo in files:
        try:
            fo.close()
        except:
            pass
            
#Concatenate files from filenames[] to fileout
def concatenate(filenames, fileout):
    with open(fileout, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


#arr[] is in the reverse ordered
#code finds the indices of the array that correspond to the pair of values <val>
def fnd_pair_of_ind_in_ordered_arr(arr, val):
    i = len(arr)-1
    ind = []
    for jj in range(2):
        j = 1 - jj
        while i > 0 and arr[i] <= val[j]:
            i -= 1
        ind.append(i)
        ind.reverse()
    return ind

            
