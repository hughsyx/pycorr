################################################
# dd.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################

import inspect 
try : import ipdb 
except : pass 


def main(input,lvl=1) :
    if lvl==0 :
        return 
    # if input is a class :
    if 'class' in str(type(input)) :
        dc(input,lvl)


# print list : 
def lst(input) :
    nel=len(input)
    dispc_help('this list has '+str(nel)+' elements','','',4)

    k=0
    for iel in input :
        dispc_help(str(k),'',input[k],4)
        k=k+1 
        if k> 10 :
            dispc('...','b','b')
            dispc_help(str(nel),'',input[-1],4)
            return 


#display class instance  :
def dc(db_sta,lvl=1) :
    dispc(str(db_sta.__class__),'y','b')
    dispc(str(db_sta.__doc__),'y','n')
    # list methods of the class 
    for func in dir(db_sta) : 
        if callable(getattr(db_sta,func)) == True  : # this a method 
            if func[0:2]!='__' :
                dispc('  '+func,'c','n')
                if type(getattr(db_sta,func).__doc__) == str  and lvl >=2 :
                    dispc('    '+getattr(db_sta,func).__doc__,'w','n')
    # list variables of the class instance :
    dd(db_sta.__dict__,tab=2, print_type=False,lvl=lvl)


# display dictionnary : 
def dd(input,tab=0, print_type=True,lvl=1) :
    """ print keys & data  of object input with fancy colors. Designed for dict """
    if print_type :
        print(type(input))
    if type(input)==dict :
        dict_=input

    #sert a aligner les : apres les nom des cles. 
    length= 0 
    for ikey in dict_ : 
        if len(ikey) > length : 
            length=len(ikey)
    length = length + tab 

    #tab/last tab sert a tabuler les cles des sous-object (qui sont des listes)
    last_tab=tab    
    for ikey in dict_ :
        type_str=str(type(dict_[ikey])).split("'")[1]
        if '.' in type_str :
            type_str = type_str.split('.')[-1]
        type_key=type(dict_[ikey])
        key_str=ikey.rjust(len(ikey)+last_tab)
        data = ' '.join(str(dict_[ikey]).split())
        dispc_help(key_str,type_str.ljust(15),'  '+data,length+last_tab)
        # if this is a list of dict 1D or 2D and lvl > 1 display the content of the list :
        if lvl< 3 : continue 
        if type(dict_[ikey])==list :
            if len(dict_[ikey])> 0 :
                if type(dict_[ikey][0]) == dict : # this a 1D list of dict
                    dd(dict_[ikey][0],tab=6,print_type=False)
                if type(dict_[ikey][0]) == list : # check is this a 2D list of dict 
                    if type(dict_[ikey][0][0]) == dict :
                        dd(dict_[ikey][0][0],tab=6,print_type=False)
        
        if type(dict_[ikey]) == dict :
            dd(dict_[ikey],tab=6,print_type=False)

def dispc_twice(text1,text2,length) :
    prefix='\x1b[1m \x1b[7m';
    str_=prefix+'\x1b['+str(36)+'m'+text1.ljust(length,' ')+':'+'\x1b[0m';
    str_=str_+'\x1b[1m \x1b['+str(36)+'m  '+text2.ljust(5,' ')+'\x1b[0m';
    print(str_)


def dispc_help(text1,text2,text3,length) :
    prefix='\x1b[1m';
    str_=prefix+'\x1b['+str(31)+'m'+text1.ljust(length,' ')+':'+'\x1b[0m';
    str_=str_+'\x1b['+str(33)+'m  '+text2.ljust(5,' ')+'\x1b[0m';
    str_=str_+'\x1b['+str(36)+'m'+text3[0:50]+'\x1b[0m';
    str_ = str_.replace('\n','')
    print(str_)

#------------------------------------------------
def dispc(text,color,attr) :
    # print text in color with attribute attr ! 
    if color =='r' :
        ncolor=31
    elif color=='g' :
        ncolor=32 
    elif color=='y' :
        ncolor=33
    elif color=='b' :
        ncolor=34
    elif color=='m' :
        ncolor=35
    elif color=='c' :
        ncolor=36
    elif color=='w' :
        ncolor=37
    elif color=='gray' :
        ncolor=38

    if attr=='b' :
        prefix='\x1b[1m';
    elif attr=='d' :
        prefix='\x1b[2m';
    elif attr=='u' :
        prefix ='\x1b[4m';
    elif attr=='blink' :
        prefix ='\x1b[5m'
    elif attr=='r' :
        prefix='\x1b[7m'
    else :
        prefix=''


    str_=prefix+'\x1b['+str(ncolor)+'m'+text+'\x1b[0m';
    print(str_)
    
    
    
    
    

