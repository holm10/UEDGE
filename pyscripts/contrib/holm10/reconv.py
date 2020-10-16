# Holm10 Nov 6 2019, created from scratch

def natsort(l): 
    from re import split
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)



def reconv(path='.'):
    ''' Reconverges all cases in subfolders in path (default current folder)'''

    from os import chdir,getcwd,walk
    from uedge.rundt import rundt
    from uedge import bbb
    from importlib import reload


    path=getcwd()
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    f=open('reconv.log','w+')
    f.write('Runnin in '+path+':\n')
    f.close

    try:
        dirs.remove('grid')
    except:    
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass


    for d in dirs:
        # Move to dir with input file
        chdir(path+"/"+d+"/data")
        print('==================')
        print('CONVERGING DIRECOTRY '+d)
        print('==================')
        import input as i
        reload(i)
        i.restore_input()
        bbb.dtreal=1e-9;bbb.exmain()

        rundt()
        hdf5_save('../solutions/'+bbb.label[0].decode('UTF-8')+'.hdf5')
        f=open(path+'/reconv.log','a')
        if bbb.iterm==1:
            f.write('Case '+d+' reconverged successfully!\n')
        else:
            f.write('Case '+d+' NOT reconverged!\n')
        f.close
    chdir(path)


