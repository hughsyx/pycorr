Bonjour Laurent,

Ici les modifications de pycorr/live/c1_md.py:

def build_db_file(file_list) :
     dd.dispc('  creating db file','y','b')
     nfile=len(file_list)
     db={}
     db['lat'] = []
     db['lon'] = []
     db['elev']= []
     db['id']  = []
     db['file_I']=np.zeros((nfile,2))
     db['ref_nstack'] = []
     npath_prev = 0
     k=0

     # adding metadata :
     ff = h5py.File(file_list[0])
     db['md_c']= h5.read_group_as_dict(ff['md_c'])
     db['in_'] = h5.read_group_as_dict(ff['in_'])
     db['nodata']={}


     for ifile in file_list :
         print(ifile)
         ff = h5py.File(ifile,'r')
         db['lat'].extend(ff['md']['lat'][:].tolist())
         db['lon'].extend(ff['md']['lon'][:].tolist())
         db['elev'].extend(ff['md']['elev'][:].tolist())
         db['id'].extend(ff['md']['id'][:].tolist())
         db['file_I'][k,:] = [npath_prev, len(db['lat'])-1]
         npath_prev = len(db['lat'])
         db['ref_nstack'].extend(np.array(ff['ref_nstack'])[0,:])
         for st1 in ff['ref']:
             #print(st1)
             if st1 not in db['nodata']:
                 db['nodata'][st1]={}
             for st2 in ff['ref/'+st1]:
                 #print(st2)
                 if st2 not in db['nodata'][st1]:
                     db['nodata'][st1][st2]={}
                 if st2 not in db['nodata']:
                     db['nodata'][st2]={}
                 if st1 not in db['nodata'][st2]:
                     db['nodata'][st2][st1]={}
                 for comp in ff['ref/'+st1+'/'+st2]:
                     #print(comp)
db['nodata'][st1][st2][comp]=np.count_nonzero(ff['ref/'+st1+'/'+st2+'/'+comp][:])!=0
db['nodata'][st2][st1][comp]=db['nodata'][st1][st2][comp]
                     #print(db['nodata'][st1][st2][comp])
                     if db['nodata'][st1][st2][comp]:
                         print(st1,' - ',st2,': ',comp)
         ff.close()
         k=k+1

     # converting list into array :
     db['id']  = np.array(db['id'])
     db['lon'] = np.array(db['lon'])
     db['lat'] = np.array(db['lat'])
     db['elev']= np.array(db['elev'])

     #computing distance and baz btw each station pair :
     npath = len(db['lat'])
     db['dist']= np.zeros((npath,1))
     db['az']  = np.zeros((npath,1))
     db['baz'] = np.zeros((npath,1))
     db['couples']={}
     db['id']=db['id'].astype('U13')
     for ipath in range(0, npath) :
         [dist,az,baz] =
gps2dist(db['lat'][ipath][0],db['lon'][ipath][0],db['lat'][ipath][1],db['lon'][ipath][1])
         db['dist'][ipath] = dist/1000.
         db['az'][ipath]   = az
         db['baz'][ipath]  = baz

         sta1=db['id'][ipath,0]
         sta2=db['id'][ipath,1]

         if sta1 in db['couples'].keys():
             db['couples'][sta1][sta2]=[ipath,False]
         else:
             db['couples'][sta1]={}
             db['couples'][sta1][sta2]=[ipath,False]
         if sta2 in db['couples'].keys():
             db['couples'][sta2][sta1]=[ipath,True]
         else:
             db['couples'][sta2]={}
             db['couples'][sta2][sta1]=[ipath,True]
     db['ref_nstack']  = np.array(db['ref_nstack'])
     db['md_c']['cmp']=db['md_c']['cmp'].astype('U13')
     db['in_']['cc_dtype']=db['in_']['cc_dtype'].decode('UTF-8')
     db['in_']['cc_func']=db['in_']['cc_func'].decode('UTF-8')
     db['in_']['tag']=db['in_']['tag'].decode('UTF-8')
     db['in_']['path']=db['in_']['path'].astype('U13')
     db['in_']['cc_cmp']=db['in_']['cc_cmp'].astype('U13')

     return db

