from wx import * 
 
from wx.py.shell import Shell
from wx.py.crust import Crust
 
class ShellFrame(wx.Panel): 
    def __init__(self, parent=None, id=-1,locals=None ): 
        wx.Panel.__init__(self, parent, id) 
        self.shell = Shell(parent=self,locals=locals)
        #self.crust= Crust(parent=self,locals=locals)
        #self.shell = self.crust.shell
        #self.shell.SetSize(wxSize(800,100))
        EVT_SIZE(self,self.OnSize)

    def OnSize(self,event):
        self.shell.SetSize(self.GetClientSizeTuple())
        #self.crust.SetSize(self.GetClientSizeTuple())
class App(wx.App): 
    def OnInit(self): 
        self.frame = ShellFrame() 
        self.frame.Show(true) 
        self.SetTopWindow(self.frame) 
        return true 
 
def main(): 
    application = App(0) 
    application.MainLoop() 
 
if __name__ == '__main__': 
    main() 
