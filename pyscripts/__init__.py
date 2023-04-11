from .uedge import *
try:
    from uedge.__version__ import __version__
    from uedge.__src__ import __src__
    import uedge.checkver
except:
    try:
        from __version__ import __version__
        from __src__ import __src__
        import checkver
    except:
        __version__ = 'unknown'
print('# UEDGE version:',__version__)        
from .UEDGESettings import *
#from .UEDGESimulation import *    
#from .UEDGEDoc import *
#from .UEDGEToolBox import *       
#from .UEDGEBas2Py import *
#__all__=["UEDGEToolBox","UEDGESettings"]


#
# Load the startup file .uedgerc.py from cwd or home.
#
_homepath = path.expanduser('~')
_homefile = Path('{}/.uedgerc.py'.format(_homepath))
_localpath = path.expanduser('.')
_localfile = Path('{}/.uedgerc.py'.format(_localpath))

if path.exists(_localfile):
   with open(_localfile) as f:
      exec(open(_localfile).read())
elif path.exists(_homefile):
   with open(_homefile) as f:
      exec(open(_homefile).read())
      
