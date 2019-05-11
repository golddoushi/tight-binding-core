class quickSort:
    #####  WW   WWW   WW   ####
    #####   WW WW WW WW    ####
    #####    WW    WW      ####
    #=========================#
    # the input array will be changed forever!!!
    #=========================#
    sort_list=[] #the index enumerated with data list
    input_list=[] #the input initial data list
    
    def __init__(self):
        pass
    
    def sort(self,input_list=[]):
        #print '>>>>>>>>>>>>'
        #print 'qsort reporting...'
        #print 'length of input array is: ',len(input_list)
        sort_list=[]
        high=len(input_list)-1
        
        for i in range(0,high+1):
            sort_list.append(i)
            
        # print self.sort_list[high]
        self.firstSort(input_list,sort_list,
                       0,high)
            
        return sort_list
        
    def firstSort(self,input_list,sort_list,low,high):    
        
        l=low;h=high
        split=input_list[low]
        
        while(l<h):
            while(input_list[h]>split):
                h-=1
                
            if(l<h):
                input_list[l]=input_list[h]
                input_list[h]=split
                
                sort_tmp=sort_list[l]
                sort_list[l]=sort_list[h]
                sort_list[h]=sort_tmp
                
                l+=1
            
            while(input_list[l]<split):
                l+=1
            
            if l<h :
                input_list[h]=input_list[l]
                input_list[l]=split
                
                sort_tmp=sort_list[l]
                sort_list[l]=sort_list[h]
                sort_list[h]=sort_tmp
                
                h-=1
        
        self.loop(input_list,sort_list,l,h,low,high)
        
    def loop(self,input_list,sort_list,l,h,low,high):
        if l>low :
            self.firstSort(input_list,sort_list,low,h-1)
            
        if h<high:
            self.firstSort(input_list,sort_list,l+1,high)