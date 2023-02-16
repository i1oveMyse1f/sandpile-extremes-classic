import numpy as np

#dens2smooth, version 25.05.2022
#smoothed values are in the 3d column
#time (1st columnt) is assumed to be integers
def dens2smooth(iname, wnd=1, gap=1):
    def incr(n, wnd): #return n + 1 (mod wnd)
        n += 1
        if n >= wnd:
            n = 0
        return n

    def adjust_wnd(time_fst_in_wnd, time_cur, wnd):
        while time_fst_in_wnd <= time_cur - wnd:
            time_fst_in_wnd += wnd
        time2save = time_fst_in_wnd + 0.5*(wnd-1) #time of the window center
        time_lst_in_wnd = time_fst_in_wnd + wnd - 1
        return time_fst_in_wnd, time_lst_in_wnd, time2save

    with open(iname, 'r') as f:
        time = []
        dens = []
        time_fst_in_wnd = 0
        cur_dens = 0
        n_in_wnd = 0 #number of current window elements
        wnd_time = np.zeros(wnd+gap-1)
        wnd_dens = np.zeros(wnd+gap-1)
        i_in = 0 #index of the first window element (first in first out)
        i_out = 0
        isfst = 1
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            time_cur = float(ll[0])
            size = int(float(ll[1]))
            dens_val = float(ll[2])
            if isfst: #window is empty
                time_fst_in_wnd = time_cur
                time_lst_in_wnd = time_fst_in_wnd + wnd - 1
                time2save = time_fst_in_wnd + 0.5*(wnd-1)
                isfst = 0
            if n_in_wnd == 0 and time_cur < time_fst_in_wnd:
                #possible if gap > wnd
                continue
            elif n_in_wnd == 0 and time_cur > time_lst_in_wnd:
                time_fst_in_wnd, time_lst_in_wnd, time2save = \
                    adjust_wnd(time_fst_in_wnd, time_cur, wnd)
            if time_cur > time_lst_in_wnd:
                #Save window data
                dens.append(cur_dens)
                time.append(time2save)
                #Definition of the next window
                time_fst_in_wnd += gap
                time_fst_in_wnd, time_lst_in_wnd, time2save = \
                    adjust_wnd(time_fst_in_wnd, time_cur, wnd)
                #--The window is defined
                #two cases: 
                #(i) time_fst_in_wnd < time_cur; 
                #(ii) time_fst_in_wnd >= time_cur (possible if gap > wnd); 
            #Delete the earlierst data from the window
            while n_in_wnd > 0 and wnd_time[i_out] < time_fst_in_wnd:
                n_in_wnd -= 1
                if n_in_wnd == 0:
                    cur_dens = 0
                elif n_in_wnd < 0:
                    print('Negative number of elements in the window')
                    print('Current time\n', time_cur)
                    return time, dens
                else:
                    cur_dens -= (wnd_dens[i_out] - cur_dens) / n_in_wnd
                i_out = incr(i_out, wnd)
            #Add to window
            wnd_dens[i_in] = dens_val
            wnd_time[i_in] = time_cur
            i_in = incr(i_in, wnd)
            n_in_wnd += 1
            #Find the average
            cur_dens += (dens_val - cur_dens)/n_in_wnd
            #print(time_cur, dens_val, cur_dens)
   
    return time, dens


#smooth the data from the second column over wnd time moment
#the first column is time; gap is the time difference between the windows
#wnd_max is the maximal number of elements gathered within a single window
#improvment of dens2smooth(): the maximal number of smoothed values in sliding windows is controlled with wnd_max
def arr2smooth(iname, wnd=1, gap=1, wnd_max=1000):
    def incr(n, mod): #return n + 1 (mod)
        n += 1
        if n >= mod:
            n = 0
        return n

    def adjust_wnd(time_fst_in_wnd, time_cur, wnd):
        while time_fst_in_wnd <= time_cur - wnd:
            time_fst_in_wnd += wnd
        time2save = time_fst_in_wnd + 0.5*(wnd-1) #time of the window center
        time_lst_in_wnd = time_fst_in_wnd + wnd - 1
        return time_fst_in_wnd, time_lst_in_wnd, time2save

    with open(iname, 'r') as f:
        time = []
        arr = []
        time_fst_in_wnd = 0
        cur_arr = 0
        n_in_wnd = 0 #number of current window elements
        wnd_time = np.zeros(wnd_max)
        wnd_arr = np.zeros(wnd_max)
        i_in = 0 #index of the first window element (first in first out)
        i_out = 0
        isfst = 1
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            time_cur = float(ll[0])
            arr_val = float(ll[1])
            if isfst: #window is empty
                time_fst_in_wnd = time_cur
                time_lst_in_wnd = time_fst_in_wnd + wnd - 1
                time2save = time_fst_in_wnd + 0.5*(wnd-1)
                isfst = 0
            if n_in_wnd == 0 and time_cur < time_fst_in_wnd:
                #possible if gap > wnd
                continue
            elif n_in_wnd == 0 and time_cur > time_lst_in_wnd:
                time_fst_in_wnd, time_lst_in_wnd, time2save = \
                    adjust_wnd(time_fst_in_wnd, time_cur, wnd)
            if time_cur > time_lst_in_wnd:
                #Save window data
                arr.append(cur_arr)
                time.append(time2save)
                #Definition of the next window
                time_fst_in_wnd += gap
                time_fst_in_wnd, time_lst_in_wnd, time2save = \
                    adjust_wnd(time_fst_in_wnd, time_cur, wnd)
                #--The window is defined
                #two cases: 
                #(i) time_fst_in_wnd < time_cur; 
                #(ii) time_fst_in_wnd >= time_cur (possible if gap > wnd); 
            #Delete the earlierst data from the window
            #print(n_in_wnd, wnd_time[i_out], time_fst_in_wnd)
            while n_in_wnd > 0 and wnd_time[i_out] < time_fst_in_wnd:
                n_in_wnd -= 1
                if n_in_wnd == 0:
                    cur_arr = 0
                elif n_in_wnd < 0:
                    print('Current time\n', time_cur)
                    raise ValueError('Negative number of elements in the window')
                else:
                    cur_arr -= (wnd_arr[i_out] - cur_arr) / n_in_wnd
                i_out = incr(i_out, wnd_max)
            #Add to window
            wnd_arr[i_in] = arr_val
            wnd_time[i_in] = time_cur
            i_in = incr(i_in, wnd_max)
            n_in_wnd += 1
            if n_in_wnd > wnd_max: #error: two many elements in the window
                raise ValueError('arr2smooth: error: two many elements in the window')
            #Find the average
            cur_arr += (arr_val - cur_arr)/n_in_wnd
   
    return time, arr

#collect sizes that are larger than <size_border> from <fname>
def tail_from_file(fname, size_border):
    is_frst_line = True
    size = []
    with open(fname, 'r') as f:
        for line in f:
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            s = int(round(float(ll[1])))
            if is_frst_line:
                t_frst = float(ll[0])
                is_frst_line = False
            if s > size_border:
                size.append(s)
    t_last = float(ll[0])
    t_lng = t_last-t_frst+1
    size.sort(reverse=True)
    return size, t_lng

#tail[] is the list of sizes placed in the reversed order
#t_lng is the normalization factor
#fname is the data file name
#Code writes sizes one-by-one
#  as a new written size coincides with the previous one, their number is increased by 1
#  otherwise the previous size and its quantity is saved 
def tail_cum_dstr_to_file(tail, t_lng, fname):
    is_new_val = False
    s4dst = []
    quant = []
    cum = []
    cum_compl = []
    s_prv = tail[0]
    q_cur = 0
    cum_prv = 0
    for s in tail:
        if s < s_prv:
            s4dst.append(s_prv)
            quant.append(q_cur)
            cum_prv += q_cur
            cum.append(1-cum_prv/t_lng)
            cum_compl.append(cum_prv/t_lng)
            q_cur = 1
        else:
            q_cur += 1
        s_prv = s
        
    with open(fname, 'w') as f:
        f.write(f'#size\tQuant\tQua_cum\tTime_length: {t_lng}\n\n')
        for i in range(len(s4dst)):
            f.write(f'{s4dst[i]}\t{quant[i]}\t{cum_compl[i]}\n')
    return s4dst, quant, cum, cum_compl

#Reads the tail cummulative distribution from file
#This is inverse to tail_cum_dstr_to_file()
def read_cum_dstr(fname):
    s4dst = []
    quant = []
    cum_compl = []
    with open(fname, 'r') as f:
        for line in f:
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            if len(ll) >= 3:
                s4dst.append(int(round(float(ll[0]))))
                quant.append(int(ll[1]))
                cum_compl.append(float(ll[2]))
    return s4dst, quant, cum_compl

#Reads the time length of the catalogue from the first line of the input file
#Assigns True to <is_err> and np.NaN to t_lng, if not found
import numpy as np
def read_time_lng(fname):
    with open(fname, 'r') as f:
        for line in f:
            fragment = 'Time_length: '
            ind = line.find(fragment)
            if line[0] != '#' or ind < 0:
                is_err = True; t_lng = np.NaN
            else:
                t_lng = int(float(line[ind+len(fragment):]))
                is_err = False
            break
    return is_err, t_lng

#Input file: 
#Lines: Time<tab>Size Any_symbols or nothing<EOL>
#Lines starting with # are treated as comments and ignored
import re 
def sfr(iname, smax, smin=1, gap=1.2, issave=1, istimenorm=1, pattern='sfr', oname=None):
    lgap = np.log10(gap)
    lgsmax = np.log(smax)/np.log(10)
    lgsmin = np.log(smin)/np.log(10)
    num_of_size_arr = int(np.floor((lgsmax-lgsmin)/lgap))+1 #array length
    
    # Insert <pattern> after the last slash / in <s>.
    # If the slashes are absent, <pattern> is added to the beginning of <s>
    def ins_sfr_to_fname(s, pattern='sfr'):
        slash = '/'
        inds = [_.start() for _ in re.finditer(slash, s)]
        if len(inds) > 0:
            slash_last_ind = inds[-1]
            snew = f'{s[:slash_last_ind+1]}{pattern}{s[slash_last_ind+1:]}'
        else:
            snew = f'{pattern}{s}'
        return snew
    
    if oname is None:
        oname = ins_sfr_to_fname(iname, pattern)
    # --End of definitions
    
    is_frst = True
    t_frst = 0
    with open(iname, 'r') as f:
        x = [10 ** (lgsmin + lgap*(i+0.5)) for i in range(num_of_size_arr)]
        q = [0 for _ in range(num_of_size_arr)]
        s_out_of_arr = 0
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            time = int(float(ll[0]))
            if is_frst:
                t_frst = time #Used to find the length of the time interval: t_last-t_frst+1
                is_frst = False
            size = int(np.round(float(ll[1])))
            if size <= 0:
                s_out_of_arr += 1
            else:
                ind = int(np.floor((np.log10(size) - lgsmin) / lgap)) #index of the bin
                if (ind >= 0 and ind < num_of_size_arr):
                    q[ind] += 1 #quantity of this bin increases
                else:
                    s_out_of_arr += 1

    total_event_num = sum(q) + s_out_of_arr
    t_lng = time - t_frst + 1
    if t_lng == 0:
        print('Unexpectedly, the length of the interval is 0')
        t_lng = 1
    #--Normalization and elimination of empty bins
    freq = []
    s4pos = []
    for ind in range(num_of_size_arr):
        if q[ind] > 0:
            s4pos.append(smin*np.exp(np.log(10)*ind*lgap))
            if istimenorm:
                freq.append(q[ind] / t_lng)
            else:
                freq.append(q[ind] / total_event_num)

    if issave:
        with open(oname, 'w') as f:
            f.write('Size\tFrequency\n')
            for i in range(len(freq)):
                f.write('{:.1f}\t{:.16f}\n'.format(s4pos[i], freq[i]))

    return s4pos, np.array(freq), t_lng, q, s_out_of_arr

#Exact empirical pdf
#input file <iname> (float)time (float)size anything
#Assumed that size is integer (but can be written in a float format);
#only size < smax are used in the construction
#oname = <oname> makes the function write the output to the specified file
import numpy as np
def sfr_exact(iname, smax=100, oname=None):
    arr = np.zeros(smax)
    s_max_real = 0
    is_first = True
    with open(iname, 'r') as f:
        sizes = [int(x) for x in f.readline().split()]
        for i, s in enumerate(sizes):
            if is_first:
                t_frst = i
                is_first = False
            if s < smax and s >= 0:
                arr[s] += 1
                if s > s_max_real:
                    s_max_real = s

    t_last = i
    t_lng = t_last - t_frst
    s_out = []
    freq = []
    s_arr_num = min(s_max_real, smax)
    for i in range(s_arr_num):
        if arr[i] > 0:
            s_out.append(i)
            freq.append(arr[i] / t_lng)
    if oname is not None:
        with open(oname, 'w') as f:
            f.write('Size\tFrequency\n')
            for i in range(len(freq)):
                f.write('{:.1f}\t{:.16f}\n'.format(s_out[i], freq[i]))
    return s_out, freq, t_lng

#average the density function over the bins with the length that is uniform in the logarithmic scale
#input file <iname> (float)time (float)size anything
#If issave=1 and oname is not specified then <pattern> is added to the begining of <iname>
#default normalization of frequencies, istimenorm=1, is by the number of avalanches computed as the difference between the last and first elements of the first column
import re 
def sfr_av(iname, smax, smin=1, gap=1.2, issave=1, istimenorm=1, isminone=1, pattern='pdf_av', oname=None):
    lgap = np.log10(gap)
    lgsmax = np.log(smax)/np.log(10)
    lgsmin = np.log(smin)/np.log(10)
    num_of_size_arr = int(np.floor((lgsmax-lgsmin)/lgap))+1 #array length
    
    # Insert <pattern> after the last slash / in <s>.
    # If the slashes are absent, <pattern> is added to the beginning of <s>
    def ins_sfr_to_fname(s, pattern='pdf_av'):
        slash = '/'
        inds = [_.start() for _ in re.finditer(slash, s)]
        if len(inds) > 0:
            slash_last_ind = inds[-1]
            snew = f'{s[:slash_last_ind+1]}{pattern}{s[slash_last_ind+1:]}'
        else:
            snew = f'{pattern}{s}'
        return snew
    
    if oname is None:
        oname = ins_sfr_to_fname(iname, pattern)
    # --End of definitions
    
    is_frst = True
    t_frst = 0
    with open(iname, 'r') as f:
        x = [10 ** (lgsmin + lgap*(i+0.5)) for i in range(num_of_size_arr)]
        q = [0 for _ in range(num_of_size_arr)]
        s_out_of_arr = 0
        num_of_lines = 0
        sizes = [int(x) for x in f.readline().split()]
        for i, size in enumerate(sizes):
            num_of_lines += 1
            time = i
            if is_frst:
                t_frst = time #Used to find the length of the time interval: t_last-t_frst+1
                is_frst = False
            if size <= 0:
                s_out_of_arr += 1
            else:
                ind = int(np.floor((np.log10(size) - lgsmin) / lgap)) #index of the bin
                if (ind >= 0 and ind < num_of_size_arr):
                    q[ind] += 1 #quantity of this bin increases
                else:
                    s_out_of_arr += 1

    total_event_num = sum(q) + s_out_of_arr
    t_lng = time - t_frst + 1
    if t_lng == 0:
        print('Unexpectedly, the length of the interval is 0')
        t_lng = 1
    #--Normalization and elimination of empty bins
    freq = []
    s4pos = []
    for ind in range(num_of_size_arr):
        bin_lng = smin*gap**ind*(gap-1)
        if isminone and bin_lng < 1:
            bin_lng = 1
        if q[ind] > 0:
            s4pos.append(x[ind])
            if istimenorm:
                freq.append(q[ind] / t_lng / bin_lng)
            else:
                freq.append(q[ind] / total_event_num / bin_lng)

    if issave:
        with open(oname, 'w') as f:
            f.write('Size\tFrequency\n')
            for i in range(len(freq)):
                f.write('{:.1f}\t{:.16f}\n'.format(s4pos[i], freq[i]))

    return s4pos, np.array(freq), t_lng, q, s_out_of_arr

import re
#Input file: Time <TAB> Value
#Values are smoothed over the windows with the constant logarithmic bins related to Time
def smooth_log(iname, tmin=1, gap=1.2, issave=0, pattern='smth_log_'):
    
    def ins_sfr_to_fname(s, pattern='smth_log_'):
        slash = '/'
        inds = [_.start() for _ in re.finditer(slash, s)]
        if len(inds) > 0:
            slash_last_ind = inds[-1]
            snew = f'{s[:slash_last_ind+1]}{pattern}{s[slash_last_ind+1:]}'
        else:
            snew = f'{pattern}{s}'
        return snew
    
    oname = ins_sfr_to_fname(iname, pattern)
    # --End of definitions

    is_frst = True
    with open(iname, 'r') as f:
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            ll = line.split()
            time = float(ll[0])
            if is_frst:
                tmax = time #Used to find the length of the bin
                lgap = np.log10(gap)
                lgtmax = np.log(tmax)/np.log(10)
                lgtmin = np.log(tmin)/np.log(10)
                lng_arr = int(np.floor((lgtmax-lgtmin)/lgap))+1 #array length
                x = [10 ** (lgtmin + lgap*(i+0.5)) for i in range(lng_arr)]
                val_av = [0 for _ in range(lng_arr)] #Average value in the bin
                n_cur = [0 for _ in range(lng_arr)] #Current number of values in the bin
                val_out_of_arr = 0
                is_frst = False
            val_cur = float(ll[1])
            if time < tmin or time > tmax:
                val_out_of_arr += 1
            else:
                ind = int(np.floor((np.log10(time) - lgtmin) / lgap)) #index of the bin
                if (ind >= 0 and ind < lng_arr):
                    n_cur[ind] += 1 #quantity of this bin increases
                    val_av[ind] += (val_cur - val_av[ind])/n_cur[ind]
    #----Array val_av[] with the smoothed values is created
    #Eliminate zeros from val_av[] and create the corresponding array with Time
    val_av_new = []
    time_new = []
    for ind in range(lng_arr):
        if n_cur[ind] > 0:
            #time_new.append(tmin*np.exp(np.log(10)*ind*lgap))
            time_new.append(x[ind])
            val_av_new.append(val_av[ind])
    #----Done
    
    if issave:
        with open(oname, 'w') as f:
            f.write('Time\tSmoothed_Value\n')
            for i in range(len(time_new)):
                f.write('{:f}\t{:f}\n'.format(time_new[i], val_av_new[i]))

    return time_new, val_av_new

import re
#Input file: Time <TAB> Value
#Values are summed up over the windows with the constant logarithmic bins related to Time
def sum_log(iname, tmin=1, gap=1.2, issave=0, pattern='sum_log_'):
    
    def ins_sfr_to_fname(s, pattern='smth_log_'):
        slash = '/'
        inds = [_.start() for _ in re.finditer(slash, s)]
        if len(inds) > 0:
            slash_last_ind = inds[-1]
            snew = f'{s[:slash_last_ind+1]}{pattern}{s[slash_last_ind+1:]}'
        else:
            snew = f'{pattern}{s}'
        return snew
    
    oname = ins_sfr_to_fname(iname, pattern)
    # --End of definitions

    is_frst = True
    with open(iname, 'r') as f:
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            ll = line.split()
            time = float(ll[0])
            if is_frst:
                tmax = time #Used to find the length of the time interval: t_last-t_frst+1
                lgap = np.log10(gap)
                lgtmax = np.log(tmax)/np.log(10)
                lgtmin = np.log(tmin)/np.log(10)
                lng_arr = int(np.floor((lgtmax-lgtmin)/lgap))+1 #array length
                x = [10 ** (lgtmin + lgap*(i+0.5)) for i in range(lng_arr)]
                val_sum = [0 for _ in range(lng_arr)] #Sum of values in the bin
                n_cur = [0 for _ in range(lng_arr)] #Current number of values in the bin
                val_out_of_arr = 0
                is_frst = False
            val_cur = float(ll[1])
            if time < tmin or time > tmax:
                val_out_of_arr += 1
            else:
                ind = int(np.floor((np.log10(time) - lgtmin) / lgap)) #index of the bin
                if (ind >= 0 and ind < lng_arr):
                    n_cur[ind] += 1 #quantity of this bin increases
                    val_sum[ind] += val_cur
    #----Array val_sum[] with the smoothed values is created
    #Eliminate zeros from val_sum[] and create the corresponding array with Time
    val_sum_new = []
    time_new = []
    for ind in range(lng_arr):
        if n_cur[ind] > 0:
            time_new.append(x[ind])
            val_sum_new.append(val_sum[ind])
    #----Done
    
    if issave:
        with open(oname, 'w') as f:
            f.write('Time\tSummed_Value\n')
            for i in range(len(time_new)):
                f.write('{:f}\t{:f}\n'.format(time_new[i], val_sum_new[i]))

    return time_new, val_sum_new

#Input file: Time <TAB> Value
#Values are summed up over the windows with the constant logarithmic bins related to Time and divided by the LENGTH of the bin
def av_log(iname, tmin=1, gap=1.2, issave=0, pattern='av_log_', oname=None):
    
    def ins_sfr_to_fname(s, pattern=pattern):
        slash = '/'
        inds = [_.start() for _ in re.finditer(slash, s)]
        if len(inds) > 0:
            slash_last_ind = inds[-1]
            snew = f'{s[:slash_last_ind+1]}{pattern}{s[slash_last_ind+1:]}'
        else:
            snew = f'{pattern}{s}'
        return snew
    # --End of definitions

    is_frst = True
    with open(iname, 'r') as f:
        num_of_lines = 0
        for line in f:
            num_of_lines += 1
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            ll = line.split()
            time = float(ll[0])
            if is_frst:
                tmax = time #Used to find the log length of the bins
                lgap = np.log10(gap)
                lgtmax = np.log(tmax)/np.log(10)
                lgtmin = np.log(tmin)/np.log(10)
                lng_arr = int(np.floor((lgtmax-lgtmin)/lgap))+1 #array length
                x = [10 ** (lgtmin + lgap*(i+0.5)) for i in range(lng_arr)]
                val_sum = [0 for _ in range(lng_arr)] #Sum of values in the bin
                n_cur = [0 for _ in range(lng_arr)] #Current number of values in the bin
                val_out_of_arr = 0
                is_frst = False
            val_cur = float(ll[1])
            if time < tmin or time > tmax:
                val_out_of_arr += 1
            else:
                ind = int(np.floor((np.log10(time) - lgtmin) / lgap)) #index of the bin
                if (ind >= 0 and ind < lng_arr):
                    n_cur[ind] += 1 #quantity of this bin increases; not used
                    val_sum[ind] += val_cur
    #----Array val_sum[] with the smoothed values is created
    #Eliminate zeros from val_sum[] and create a new array deviding each value by the bin length
    val_av = []
    time_new = []
    for ind in range(lng_arr):
        bin_lng = tmin*gap**ind*(gap-1)
        if bin_lng < 1:
            bin_lng = 1
        if n_cur[ind] > 0:
            time_new.append(x[ind])
            val_av.append(val_sum[ind] / bin_lng)
    #----Done
    
    if issave:
        if oname is None:
            oname = ins_sfr_to_fname(iname, pattern)
        with open(oname, 'w') as f:
            f.write('Time\tAv_Value\n')
            for i in range(len(time_new)):
                f.write('{:f}\t{:f}\n'.format(time_new[i], val_av[i]))

    return time_new, val_av

def split2two(iname, oname1, oname2):
    #Find time_break;
    isfrstline = True
    with open(iname, 'r') as f:
        for line in f:
            if not(line[0].isdigit()):
                continue
            if isfrstline:
                ll = line.split()
                try:
                    time_frst = int(ll[0])
                except:
                    print('Incorrect symbols in the first data line')
                isfrstline = False
        try:
            time_last = int(line.split()[0])
        except:
            print('Incorrect symbols in the last data line')
        time_break = int(round(0.5*(time_frst+time_last)))
    #--found
    with open(iname, 'r') as f, open(oname1, 'w') as fo1, open(oname2, 'w') as fo2:
        ishead = 1
        head = ''#head of the input file is collected and saved into both output files
        isheadtofrst = 1
        isheadtoscnd = 1
        for line in f:
            if ishead:
                if not(line[0].isdigit()):
                    head = head + line
                    continue
                ishead = 0
            if not(line[0].isdigit()):
                continue
            ll = line.split()
            time = int(ll[0])
            if time < time_break:
                if isheadtofrst == 1:
                    fo1.write('{:s}'.format(head))
                    isheadtofrst = 0
                ln_split = line.split()
                if len(ln_split) >= 3:
                    fo1.write(f'{ln_split[0]}\t{ln_split[1]}\t{ln_split[2]}\n')
            else:
                if isheadtoscnd == 1:
                    fo2.write('{:s}'.format(head))
                    isheadtoscnd = 0
                ln_split = line.split()
                if len(ln_split) >= 2:
                    fo2.write(f'{ln_split[0]}\t{ln_split[1]}\t{ln_split[2]}\n')
            

# Define function for string formatting of scientific notation
from math import floor, log10
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)

