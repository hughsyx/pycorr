import inspect 


def dd(input,lvl=1) :
	if lvl==0 :
		return 
	# if input is a class :
	if type(input) == dict :
		_disp_dict(input,lvl)
	elif 'class' in str(type(input)) :
		_disp_class(input,lvl)
	elif type(input) == list :
		_disp_list(input,lvl)


#def disp(input,lvl=1) :
#	''' main function to disp any kind of variable'''
#	if lvl==0 :
#		return 
#	# if input is a class :
#	if type(input) ==dict :
#		_disp_dict(input,lvl)
#	elif 'class' in str(type(input)) :
#		_disp_class(input,lvl)
#	elif type(input) == list :
#		_disp_list(input,lvl)


def disp3(text1,text2,text3,c1=31,c2=33,c3=36,l3=50,l1=20,l2=20,bold=True) :
	prefix=''
	if bold :
		prefix='\x1b[1m';
	if type(c1)==str : c1 = int(c1) ; 
	if type(c2)==str : c2 = int(c2) ; 
	if type(c3)==str : c3 = int(c3) ; 


	str_=prefix+'\x1b['+str(c1)+'m'+text1.ljust(l1,' ')+':'+'\x1b[0m'
	str_=str_+'\x1b['+str(c2)+'m  '+text2.ljust(l2,' ')+'\x1b[0m'
	str_=str_+'\x1b['+str(c3)+'m'+text3[0:l3]+'\x1b[0m'
	str_ = str_.replace('\n','')
	print(str_)



def disp2(text1,text2,c1=31,c2=33,l1=20,l2=20,bold=True) :
	prefix=''
	if bold :
		prefix='\x1b[1m';
	if type(c1)==str : c1 = int(c1) ; 
	if type(c2)==str : c2 = int(c2) ; 

	str_=prefix+'\x1b['+str(c1)+'m'+text1.ljust(l1,' ')+':'+'\x1b[0m'
	str_=str_+'\x1b['+str(c2)+'m  '+text2.ljust(l2,' ')+'\x1b[0m'
	#str_=str_+'\x1b['+str(c3)+'m'+text3[0:l3]+'\x1b[0m'
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


#-----------------------------------------------------------------------------
#
#  SUB FUNCTION TO DISPLAY NICELY DIFFERENT KIND OF PYTHON VARIABLES
#
#-------------------------------------------------------------------------------


# print list : 
def _disp_list(input) :
	nel=len(input)
	disp3('this list has '+str(nel)+' elements','','')

	k=0
	for iel in input :
		disp3(str(k),'',input[k])
		k=k+1 
		if k> 10 :
			dispc('...','b','b')
			disp3(str(nel),'',input[-1])
			return 



#display class instance  :
def _disp_class(db_sta,lvl=1) :
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
	_disp_dict(db_sta.__dict__,tab=2, print_type=False,lvl=lvl)


# display dictionnary : 
def _disp_dict(input,tab=0, print_type=True,lvl=1) :
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
		#disp3(key_str,type_str.ljust(15),'  '+data,length+last_tab)
		disp3(key_str,type_str,'  '+data,l1=14,l2=10)
		# if this is a list of dict 1D or 2D and lvl > 1 display the content of the list :
		if lvl< 3 : continue 
		if type(dict_[ikey])==list :
			if len(dict_[ikey])> 0 :
				if type(dict_[ikey][0]) == dict : # this a 1D list of dict
					disp_dict(dict_[ikey][0],tab=6,print_type=False)
				if type(dict_[ikey][0]) == list : # check is this a 2D list of dict 
					if type(dict_[ikey][0][0]) == dict :
						disp_dict(dict_[ikey][0][0],tab=6,print_type=False)
        
		if type(dict_[ikey]) == dict :
			_disp_dict(dict_[ikey],tab=6,print_type=False)


