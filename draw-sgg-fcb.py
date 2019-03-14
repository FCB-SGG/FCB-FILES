#! python3
# coding:utf-8
#you can download the FCB files of SGG on: https://github.com/FCB-SGG/FCB-FILES

#author:Pan Li(lipan.whu@gmail.com)


import os
import math
import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys, getopt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as mdates
from matplotlib.pylab import date2num

FRAC_VAL_GROSS = 0.2

class fcb_base_t:
    def __init__(self):
        self.id = 'X00'
        self.csys  = 'X'
        self.prn   = 0
        self.val   = 0.0


class fcb_1ep:
    def __init__(self):
        self.tm   = ''
        self.inf_fbs = []


def read_fcbfile(path):
    fp = open(path)
    f_Eps =[]

    ln = fp.readline()
    while ln:
        if '*' == ln[0]:

            tm = datetime.datetime.strptime(ln[2:-1 if len(ln)<30 else 30].strip(), '%Y %m %d %H %M %S.%f')

            fEp = fcb_1ep()
            fEp.tm = tm

            while 1:
                ln = fp.readline()
                if not ln:
                    break
                if '*' == ln[0]:
                    break
                if 'P' != ln[0]:
                    continue
                strs = ln.split()
                if 3 > len(strs):
                    raw_input('wrong line 3>len(strs)' + ln)

                ss = strs[0][1:]
                ch = ss[0]
                prn = int(strs[0][2:])

                fb = fcb_base_t()
                fb.id = ss
                fb.csys = ch
                fb.prn = prn
                fb.val = float(strs[1])

                fEp.inf_fbs.append(fb)

            if 3 < len(fEp.inf_fbs):
                f_Eps.append(fEp)
        else:
            ln = fp.readline()

    fp.close()
    return f_Eps


def cal_ave_std(dif):
    n = len(dif)
    if 0 == n:
        return 0.0, 999.9

    ave = sum(dif) / n

    std = 0.0
    for i in range(n):
        std += (dif[i] - ave) * (dif[i] - ave)
    std = math.sqrt(std / n)

    return ave, std


def find_cent_val(dif, ave):
    k = 0

    lst = []

    for v in dif:
        tmp = v - ave
        v = v - round(tmp)
        tmp1 = v - ave

        if abs(tmp1) < FRAC_VAL_GROSS:
            lst.append(v)
            k += 1

    if 0 == k:
        return ave, 999.9, k

    a0, s0 = cal_ave_std(lst)
    return a0, s0, k


def find_cent_val_0(dif):
    n = len(dif)
    nmax = -999
    ix = -1

    for i in range(0, n):
        ave, std, k = find_cent_val(dif, dif[i])
        if k > nmax:
            nmax = k
            ix = i
    return nmax, ix


def cal_ave_std_robust(dif):
    # print "src", dif
    n = len(dif)

    d0 = dif[0]
    for i in range(1, n):
        dt = dif[i] - d0
        dif[i] -= round(dt, 0)

    nm, ix0 = find_cent_val_0(dif)
    ave0 = dif[ix0]

    ave = ave0
    std = 999.9
    nmax = -1

    for i in range(0, 5):

        ab, sb, nb = find_cent_val(dif, ave)

        bBetter = False
        if nb > nmax:
            bBetter = True
        elif sb < std:
            bBetter = True

        if bBetter:
            nmax = nb
            std = sb
            ave = ab
        else:
            break

        if nmax == n:
            break
        if nmax > n * 0.9:
            if std < 0.005:
                break

    return ave, std


def adjust_l1_fcbs_(fEPs, chsys):
    nep = len(fEPs)

    for i in range(0, nep-2):
        fE0 = fEPs[i]
        fE1 = fEPs[i+1]
        nb0 = len(fE0.inf_fbs)
        nb1 = len(fE1.inf_fbs)

        k0 = 0
        k1 = 0

        dif = []
        ix0 = []
        ix1 = []
        difmax = 0.0
        while k0<nb0 and k1<nb1:
            fb0 = fE0.inf_fbs[k0]
            fb1 = fE1.inf_fbs[k1]
            bSkp = False
            if fb0.csys!=chsys:
                k0+=1
                bSkp=True
            if fb1.csys!=chsys:
                k1+=1
                bSkp=True
            if fb0.prn > fb1.prn:
                k1+=1
                bSkp=True
            if fb0.prn < fb1.prn:
                k0+=1
                bSkp=True
            if bSkp:
                continue

            ix0.append(k0)
            ix1.append(k1)
            d0 = fb1.val-fb0.val
            dif.append(d0)
            if abs(d0)>abs(difmax):
                difmax = d0

            k0+=1
            k1+=1

        if len(dif)<=1:
            continue

        if abs(difmax)<0.02:
            continue
        ave, std = cal_ave_std_robust(dif)

        for k1, d0 in zip(ix1, dif):
            d1 = d0 - ave
            fE1.inf_fbs[k1].val -= ave
            fE1.inf_fbs[k1].val -= round(d1)


def draw_l1_fcbs_(fEPs, sSatLst, plt):

    xtm = []
    for ss in sSatLst:
        xp=[]
        yp=[]

        for fE0 in fEPs:
            for fb in fE0.inf_fbs:
                if fb.id == ss:
                    xp.append(fE0.tm)
                    yp.append(fb.val)

        if len(yp)<3:
            continue

        plt.plot(xp, yp, '.', label=ss)
        xtm.append(xp[0])
        xtm.append(xp[-1])


    if len(xtm)>3:
        tmflt = date2num(xtm)
        plt.xlim(tmflt.min() - 0.000001, tmflt.max() + 0.000001)
    plt.grid(True)



def draw_l1_fcbs(fEPs, pth, ch_sys):
    plt.figure(figsize=(9, 5))

    ssLst = []
    for fE0 in fEPs:
        for fb in fE0.inf_fbs:
            if fb.csys == ch_sys or fb.csys in ch_sys:
                # set the satellite you need in the figure. eg.0!=fb.prn%3 means draw a pictuee contains the fcbs of G03,G06,G09.......
                # if you want to draw a pictur of all, please comment the two lines.
                if ch_sys=='G' and 0!=fb.prn%3:
                    continue
                #####################################################

                ssLst.append(fb.id)
    ssLst = list(set(ssLst))
    ssLst.sort()

    adjust_l1_fcbs_(fEPs, ch_sys)
    draw_l1_fcbs_(fEPs, ssLst, plt)

    plt.grid(linestyle='dashed', alpha=0.4)

    lab = ''
    if 'G'==ch_sys:
        lab = 'GPS NL'
    elif 'R'==ch_sys:
        lab = 'GLO NL'
    elif 'E'==ch_sys:
        lab = 'GAL NL'
    elif 'C'==ch_sys:
        lab = 'BDS NL'
    plt.ylabel(lab + ' Cycle')

    dateFmt = mdates.DateFormatter('%H:%M')
    ax = plt.subplot(111)
    ax.xaxis.set_major_formatter(dateFmt)

    ymajorLocator = MultipleLocator(0.1)
    ax.yaxis.set_major_locator(ymajorLocator)

    plt.ylim(-1.0, 1.0)

    plt.legend(fontsize=9, ncol=5)


    fpth_jpg = pth + "-" + ch_sys + ".jpg"
    print (fpth_jpg)
    plt.savefig(fpth_jpg, dpi=360, bbox_inches='tight')
    plt.close()


def prc_fcbfile(path, navsys, dirout):

    if len(dirout)<=1:
        dirout=os.path.split(path)[0]


    print ('inputfile: ', path)
    print ('navsys   : ', navsys)
    print ('outdir   : ', dirout)
    if 0==len(navsys):
        print ('no navsys assigned. exit' % path)
        return

    if not os.path.exists(path):
        print ('file [%s] not exist.' %path)
        return
    global out_dir
    out_dir = dirout
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    fdir, fname = os.path.split(path)
    fcbs = read_fcbfile(path)       #read your FCB files
    draw_l1_fcbs(fcbs, path, navsys)    #draw a picture of NL FCB

if __name__ == "__main__":


    dir  = 'G:/GNSS/product/FCB/'   #path of your FCB files
    navsys = 'G'                    #G/E/C

    for fl in os.listdir(dir):
        if fl.endswith('.fcb'):
            for ch in 'GEC':
                prc_fcbfile(os.path.join(dir, fl), ch, '')

    print ('\n*******  Draw upd time series. Finished.  *******')



