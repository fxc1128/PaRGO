# coding:utf-8
import os
import geopandas as gpd
import queue
import math
import time

workdir = r"G: & cd egc\\pargo-dev\\wkflow & "
work = r"G:\\egc\\pargo-dev\\wkflow\\"  # 文件存放位置
nrank = 8  # 进程数
mpi = r"mpiexec -n " + str(nrank)
middata = r"data\\process\\"
finaldata = r"data\\final\\"
nbr = r"data\\neighbor\\moore.nbr "
dem = r"data\\dem.tif "

demres = middata + r"res_dem.tif "
demresfel = middata + r"res_demfel.tif "
resdirfile = middata + r"res_dir.tif "
resscafile = middata + r"res_sca.tif "
weifile = middata + r"weight.tif "
ssafile = middata + r"w_sca.tif "
resnet = middata + r"res_net.tif "
watershed = middata + r"res_watershed.tif "
netshp = middata + r"res_rivernet.shp"

demfel = finaldata + r"demfel.tif "
dirfile = finaldata + r"dir.tif "
scafile = finaldata + r"sca.tif "
netfile = finaldata + r"net.tif "
streamfile = finaldata + r"streamOrder.tif "


def largearea(_g, buf):
    exe = r" apps\\Largearead8.exe "
    cmd = workdir + mpi + exe + watershed + dirfile + nbr + scafile + str(_g) + " " + str(buf) + " "
    net = gpd.read_file(netshp)
    q = queue.Queue()
    l = len(net)
    m = max(net['strmOrder'])

    while (q.qsize() < l):
        for it in net.iloc:
            if (it.USLINKNO1 == -1 and it.USLINKNO2 == -1 and it.LINKNO != -1):
                q.put(it.LINKNO)
                index1 = net.loc[net.USLINKNO1 == it.LINKNO].index.to_list()
                index2 = net.loc[net.USLINKNO2 == it.LINKNO].index.to_list()
                iin = net.loc[net.LINKNO == it.LINKNO].index.to_list()
                net.loc[iin[0], 'LINKNO'] = -1
                if (len(index1)):
                    net.loc[index1[0], 'USLINKNO1'] = -1
                # print('up1:',net.loc[index1[0],'USLINKNO1'])
                if (len(index2)):
                    net.loc[index2[0], 'USLINKNO2'] = -1

    c = cmd + str(1)  # init
    g = os.system(c)
    t = 0
    start = time.time()
    for i in range(0, l):
        id = str(q.get())
        c = cmd + id
        print("working process: large_aread8-", id)
        g = os.system(c)
        #print(g)
        t = t + g
    end = time.time()
    print("process time:", end - start)
    return t / 10000

def largePitremove(_g, buf):
    exe = r" apps\\Largepitremove.exe "
    cmd = workdir + mpi + exe + watershed + dem + nbr + demfel + str(_g) + " " + str(buf) + " "
    net = gpd.read_file(netshp)
    q = queue.Queue()
    l = len(net)
    m = max(net['strmOrder'])

    while (q.qsize() < l):
        for it in net.iloc:
            if (it.USLINKNO1 == -1 and it.USLINKNO2 == -1 and it.LINKNO != -1):
                q.put(it.LINKNO)
                index1 = net.loc[net.USLINKNO1 == it.LINKNO].index.to_list()
                index2 = net.loc[net.USLINKNO2 == it.LINKNO].index.to_list()
                iin = net.loc[net.LINKNO == it.LINKNO].index.to_list()
                net.loc[iin[0], 'LINKNO'] = -1
                if (len(index1)):
                    net.loc[index1[0], 'USLINKNO1'] = -1
                # print('up1:',net.loc[index1[0],'USLINKNO1'])
                if (len(index2)):
                    net.loc[index2[0], 'USLINKNO2'] = -1
    print("number of watersheds:", l)
    c = cmd + str(1)
    g = os.system(c)
    start = time.time()
    for i in range(0, l):
        id = str(q.get())
        c = cmd + id
        print("working process: large_aread8-", id)
        g = os.system(c)
        print(g)
    end = time.time()
    print("process time:", end - start)

def resample(dem, resdem, _g):
    exe = r" apps\\resample.exe "
    cmd = workdir + mpi + exe + "-input " + dem + "-output " + resdem + "-g " + str(_g)
    print("working process: resample(g=", _g)
    start = time.time()
    g = os.system(cmd)
    end = time.time()
    print("process time:", end - start)

def pitremove(dem, demfel):
    exe = r" apps\\pitremove.exe "
    cmd = workdir + mpi + exe + dem + nbr + demfel
    print("working process: PitRemove of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdird8(dem, dir):
    exe = r" apps\\flowdird8.exe "
    cmd = workdir + mpi + exe + dem + nbr + dir + str(8)
    print("working process: calculate d8 of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdirdinf(dem, dinf):
    exe = r" apps\\flowdirdinf.exe "
    cmd = workdir + mpi + exe + dem + nbr + dinf
    print("working process: calculate d-inf of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scadinf(dinf, scadinf):
    exe = r" apps\\scadinf.exe "
    cmd = workdir + mpi + exe + dinf + nbr + scadinf
    print("working process: calculate sca of ", dinf)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scad8(w, dir, sca, wei=weifile):
    exe = r" apps\\scad8.exe "
    cmd = workdir + mpi + exe + dir + nbr + sca
    wcmd = workdir + mpi + exe + dir + nbr + wei + sca
    start = time.time()
    if (w):
        print("working process: calculate weighted_sca of ", dir)
        start = time.time()
        g = os.system(wcmd)
    else:
        print("working process: calculate sca of ", dir)
        start = time.time()
        g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def PeukerDouglas(dem, weifile):
    exe = r" apps\\PeukerDouglas.exe "
    cmd = workdir + mpi + exe + dem + nbr + weifile
    print("working process: catch possible river net of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def dropanalysis(nthresh, threshmin, threshmax, dem, dir):
    exe = r" apps\\dropanalysis.exe "
    cmd = workdir + mpi + exe + dem + dir + resscafile + ssafile + nbr
    print("working process: drop analysis of ", dem)
    start = time.time()
    t = threshmin
    for th in range(0, nthresh):
        r = math.exp((math.log(threshmax) - math.log(threshmin)) / (nthresh - 1))
        thresh = threshmin * pow(r, th)
        c = cmd + str(thresh)
        g = os.system(c)
        print(g)
        if (g == 0):
            t = thresh
            break
    end = time.time()
    print("process time:", end - start)
    print("thresh:",t)
    return t

def streamnet(dem, sca, src):
    exe2 = r" apps\D8Flowdir.exe "
    exe = r" apps\Streamnet.exe "
    tdir = r"data\\process\\tau_d8.tif "
    tslope = r"data\\process\\tau_sd8.tif "
    tord = r"data\\process\\tau_ord.tif "
    ttree = r"data\\process\\tau_tree.dat "
    tcoord = r"data\\process\\tau_coord.dat "
    cmd2 = workdir + mpi + exe2 + "-p " + tdir + "-sd8 " + tslope + "-fel " + dem
    cmd = workdir + mpi + exe + "-fel " + dem + "-p " + tdir + "-ad8 " + sca + "-src " + src + "-ord " + tord + "-tree " + ttree + "-coord " + tcoord + "-net " + netshp + " -w " + watershed
    print("working process: calculate watersheds of ", dem)
    g = os.system(cmd2)
    g = os.system(cmd)
    print(g)

def threshold(thresh, sca, net):
    exe = r" apps\\threshold.exe "
    cmd = workdir + mpi + exe + sca + nbr + net + str(thresh)
    print("working process: extract river net with threshold:", thresh)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def stream(net, dir, sm):
    exe = r" apps\\stream.exe "
    cmd = workdir + mpi + exe + dir + net + nbr + sm
    print("working process: calculate stream order of ", net)
    g = os.system(cmd)
    print(g)

def main(_g, buf):
    resample(dem, demres, _g)
    pitremove(demres, demresfel)
    flowdird8(demresfel, resdirfile)
    scad8(0, resdirfile, resscafile)
    PeukerDouglas(demresfel, weifile)
    scad8(1, resdirfile, ssafile, weifile)
    thresh = dropanalysis(10, 2, 500, demresfel, resdirfile)
    #print(thresh)
    threshold(thresh * 1000, resscafile, resnet)
    streamnet(demresfel, resscafile, resnet)

    largePitremove(_g, buf)
    flowdird8(demfel, dirfile)
    t = largearea(_g, buf)
    print("process time:", t)
    threshold(5000, scafile, netfile)
    #stream(netfile, dirfile, streamfile)


if __name__ == "__main__":
    _g = 10  # 重采样倍数
    buf = 2  # 缓冲区大小
    nrank=8  # 进程数
    mpi = r"mpiexec -n " + str(nrank)
    print("rank number:", nrank)
    dem = r"data\\srtm_59_07.tif "  # 输入DEM

    main(_g, buf)



