import m_pycorr.live.c1_md as c1_md

h1 = c1_md.c1_md('C1_04_xcorr_glob_storm_test__coher/RTZ')

      
#h1.plot_station_map()


h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='N',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test_Z.mat')

h1.get_selection(idvs=['*'],idvr=['*'],icmp=1,filter_p=True,p1=3.,p2=12.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test_R.mat')

h1.get_selection(idvs=['*'],idvr=['*'],icmp=2,filter_p=True,p1=3.,p2=12.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test_T.mat')






h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test0.mat',iddate=0)

h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test2.mat',iddate=2)

h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test5.mat',iddate=5)

h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test8.mat',iddate=8)

h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=15.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test11.mat',iddate=11)

# test all zeros
#h5ls -d C1_04_xcorr_glob_storm__coher/xcorr_glob_06_glob_07_000.h5/ref/YT.SILY.00/XL.HV20.00/ZZ









import m_pycorr.live.c1_md as c1_md

h1 = c1_md.c1_md('C1_04_xcorr_glob_storm_test_param_overlap2__coher/')

h1.get_selection(idvs=['AU*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=12.,ctype='P',
    bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test_overlap2.mat')






import m_pycorr.live.c1_md as c1_md

h1 = c1_md.c1_md('C1_P-PKP-PWS/RTZ')

h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=3.,p2=12.,ctype='P',
    bin_distance=True,bin_width=0.05,dist_unit='deg',file_name='./test_P_pws.mat')







h1 = c1_md.c1_md('C1_PCP_PWS/RTZ')
h1.get_selection(idvs=['*'],idvr=['*'],icmp=0,filter_p=True,p1=4.,p2=10.,ctype='P',
    bin_distance=True,bin_width=0.05,dist_unit='deg',file_name='./test_PcP_pws.mat')







import m_pycorr.live.c1_md as c1_md

h1 = c1_md.c1_md('C1_PCP_PWS/RTZ')
for ii in ['4J','AU','AW','ER','NZ','XH','YH','YT','ZJ']:
    h1.get_selection(idvs=[ii + '*'],idvr=['*'],icmp=1,filter_p=True,p1=3.,p2=12.,ctype='P',
        bin_distance=True,bin_width=0.1,dist_unit='deg',file_name='./test_PcP_pws_'+ ii + '.mat')












