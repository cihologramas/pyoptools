def convert(d):
    try:
        return int(d)
    except ValueError:
        try:
            return float(d)
        except ValueError:
            return d

def zmx_read(fn):

    f=open(fn,"rU")
    data=f.read()
    return zmx_parse(data)

##Codigo tomado y modificado de https://github.com/jordens/rayopt/blob/master/rayopt/zemax.py
from struct import Struct
import numpy as np

def zmf2dict(fn):
    """Función que lee una librería de Zemax (archivo con terminación zmf), y genera un diccionario con las descripciones
    de cada componente. La llave es la referencia de cada componente
    """
    f=open(fn,"rb")
    rd={}
    head = Struct("<I")
    lens = Struct("<100sIIIIIIIdd")
    shapes = "?EBPM"
    version, = head.unpack(f.read(head.size))
    assert version in (1001, )
    while True:
        li = f.read(lens.size)
        if len(li) != lens.size:
            if len(li) > 0:
                print(f, "additional data", repr(li))
            break
        li = list(lens.unpack(li))
        li[0] = li[0].decode("latin1").strip("\0")
        li[3] = shapes[li[3]]
        description = f.read(li[7])
        assert len(description) == li[7]
        description = zmf_obfuscate(description, li[8], li[9])
        description = description.decode("latin1")
        assert description.startswith("VERS {:06d}\n".format(li[1]))
        rd[li[0]]=description
    return rd


def zmf_obfuscate(data, a, b):
    iv = np.cos(6*a + 3*b)
    iv = np.cos(655*(np.pi/180)*iv) + iv
    p = np.arange(len(data))
    k = 13.2*(iv + np.sin(17*(p + 3)))*(p + 1)
    k = (int(("{:.8e}".format(_))[4:7]) for _ in k)
    #data = np.fromstring(data, np.uint8)
    data = np.frombuffer(data, np.uint8).copy()
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tostring()

def find_key(key,lines):
    try:
        rv=""
        for s in list(filter(lambda x:x.startswith(key),lines)):
            rv = rv+s[len(key)+1:]
    except IndexError:
        # No se encontró la llave 
        rv=""
    return rv

def checktype(surflist,t):
    for s in surflist:
        if (t in s) or (s.get('TYPE',[None])[0])==t:
            return True
    return False
    
def checkglas(surflist,t):
    for s in surflist:
        if (s.get('GLAS',[None])[0])==t:
            return True
    return False
    


def zmx2pyoptoolsSP(libdata,key):
    """Take element from a zeemax library and return a dictionary with 
    its info in a compatible way to be used in pyoptools. If the element
    is not a component composed by spherical optical surfaces return None
    """
    
    data=libdata[key]
    surfaces=data.splitlines()
    #Separar encabezado de los datos
    header=[]
    # Interpretar el encabezado
    while True:
        if surfaces[0].startswith("SURF"):
            break
        line=surfaces.pop(0)
        header.append(line)
    
    pyot={}
    lens_data={}

    description=find_key("NOTE 0",header)
    name=find_key("NAME",header)
    
    #Thor uses the NAME key for description, while edmund puts there the 
    #lens reference, and de description in the NOTE 0
    if name==key:
        lens_data["description"]=description
    else:
        lens_data["description"]=name
    
    #lens_data["description"]=name+"\n"+description
    
    #Aun no se que hacer con el mode
    mode=find_key("MODE",header)
    
    if mode !="SEQ":
        raise ValueError("MODE '{}'' not recognized".format(mode))
    
    unit=find_key("UNIT",header).split()[0]
    
    if unit!="MM":
        raise ValueError("UNIT '{}' not recognized".format(unit))
        
    gcat=find_key("GCAT",header)
    

    # Separar las superficies en una lista de diccionarios
    surflist=[]
    for line in surfaces:
        #Ignore empty lines
        if line.split()==[]:
            continue
        if line.startswith("SURF"):
            surflist.append(dict())
            continue
        line=line.lstrip()
        code=line[:4]
        data=line[5:].split()
        data=[convert(d) for d in data]
        surflist[-1][code]=data

    #Eliminar el plano objeto y el plano imagen

    surflist.pop(0)
    surflist.pop()

    
    
    # Many lenses have a first surface with no glass, this mean air, so they will be
    # removed. It seems this is user to mark the mounting tube
    if "GLAS" not in surflist[0]:
        surflist.pop(0)
        
    # Many lenses have a last surface with no glass, this mean air, so they will be
    # removed. It seems this is user to mark the mounting tube
    if "GLAS" not in surflist[-1] and "GLAS" not in surflist[-2]:
        surflist.pop()
    
    
    #Remove mirrors
    if checkglas(surflist,"MIRROR"):
        return
    
    #Don't include aspheres
    if checktype(surflist,'EVENASPH'):
        return
    
    #Don't include difractive optics
    if checktype(surflist,'BINARY_2'):
        return
    
    if checktype(surflist,'CONI'):
        return

    #Don't include Toroidal lenses. It seems cylindrical lenses are modeled as toroidal
    if checktype(surflist,'TOROIDAL'):
        return
    
    #Don't include Fresnel lenses. It seems cylindrical lenses are modeled as toroidal
    if checktype(surflist,'FRESNELS'):
        return

    if checktype(surflist,'CFRESNEL'):
        return
    
    # Identificar el tipo de lentes a partir de el numero de superficies
    # validas

    ns=len(surflist)
        
    # Normal Singlets - Need to fix for asferical lenses
    if ns==2:
        c0=surflist[0]["CURV"][0]
        c1=surflist[1]["CURV"][0]
        d0=surflist[0]["DISZ"][0]
        #Not all lenses have DIAM in both surfaces
        if "DIAM" in surflist[0] and "DIAM" in surflist[1]:
            r0=surflist[0]["DIAM"][0]
            r1=surflist[1]["DIAM"][0]
        elif "DIAM" in surflist[0]:
            r0=surflist[0]["DIAM"][0]
            r1=r0
        elif "DIAM" in surflist[1]:
            r0=surflist[1]["DIAM"][0]
            r1=r0
        else:
            raise ValueError("'DIAM' not defined in the surfaces")
        
        g0=surflist[0]["GLAS"][0]

        ##c|Verificar que las superficies son iguales, si no emitir un error
        assert r0==r1
        lens_data["material"] = g0
        lens_data["glass_catalogs"] = gcat #This is not used for the moment
        lens_data["thickness"] = surflist[0]["DISZ"][0]
        lens_data["radius"] = r0
        lens_data["curvature_s2"] = surflist[1]["CURV"][0]
        lens_data["curvature_s1"] = surflist[0]["CURV"][0]
        lens_data["type"] = "SphericalLens"
        
        pyot[key]=lens_data
        return pyot
        
        ##return CL.SphericalLens(r0,d0,c0,c1,material=m0)
    
     
    elif ns==3:
        lens_data["curvature_s1"]=surflist[0]["CURV"][0]
        lens_data["curvature_s2"]=surflist[1]["CURV"][0]
        lens_data["curvature_s3"]=surflist[2]["CURV"][0]
        lens_data["thickness_l1"]=surflist[0]["DISZ"][0]
        lens_data["thickness_l2"]=surflist[1]["DISZ"][0]
        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        
        #Verificar que las superficies son iguales, si no emitir un error
        assert r0==r1 and r1== r2
        
        lens_data["radius"] = r0
        
        lens_data["material_l1"]=surflist[0]["GLAS"][0]
        lens_data["material_l2"]=surflist[1]["GLAS"][0]
        lens_data["glass_catalogs"] = gcat 
        #return CL.Doublet(r0,c0,c1,c2, d0,d1,m0,m1)
        lens_data["type"] = "Doublet"
        pyot[key]=lens_data
        return pyot
    
    elif ns==4:
        lens_data["curvature_s1"]=surflist[0]["CURV"][0]
        lens_data["curvature_s2"]=surflist[1]["CURV"][0]
        lens_data["curvature_s3"]=surflist[2]["CURV"][0]
        lens_data["curvature_s4"]=surflist[3]["CURV"][0]

        lens_data["thickness_l1"]=surflist[0]["DISZ"][0]
        lens_data["air_gap"]=surflist[1]["DISZ"][0]
        lens_data["thickness_l2"]=surflist[2]["DISZ"][0]

        #Verificar que las superficies son iguales, si no emitir un error
        #assert r0==r1 and r1== r2 and r2 ==r3 Esto no siempre se cumple

        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        r3=surflist[3]["DIAM"][0]
        
        #Verificar que las superficies son iguales, si no emitir un error
        #assert r0==r1 and r1== r2 and r2 ==r3 #Esto no siempre se cumple
        
        lens_data["radius"] = r0

        lens_data["material_l1"]=surflist[0]["GLAS"][0]
        # g1=surflist[1]["GLAS"][0] Este parametro no existe. Es aire
        lens_data["material_l2"]=surflist[2]["GLAS"][0]

        lens_data["glass_catalogs"] = gcat 

        lens_data["type"] = "AirSpacedDoublet"
        pyot[key]=lens_data
        return pyot
    elif ns==6:
        lens_data["curvature_s1"]=surflist[0]["CURV"][0]
        lens_data["curvature_s2"]=surflist[1]["CURV"][0]
        lens_data["curvature_s3"]=surflist[2]["CURV"][0]
        lens_data["curvature_s4"]=surflist[3]["CURV"][0]
        lens_data["curvature_s5"]=surflist[4]["CURV"][0]
        lens_data["curvature_s6"]=surflist[5]["CURV"][0]
        
        lens_data["thickness_l1"]=surflist[0]["DISZ"][0]
        lens_data["thickness_l2"]=surflist[1]["DISZ"][0]
        lens_data["air_gap"]=surflist[2]["DISZ"][0]
        lens_data["thickness_l3"]=surflist[3]["DISZ"][0]
        lens_data["thickness_l4"]=surflist[4]["DISZ"][0]
        
        lens_data["material_l1"]=surflist[0]["GLAS"][0]
        lens_data["material_l2"]=surflist[1]["GLAS"][0]
        lens_data["material_l3"]=surflist[3]["GLAS"][0]
        lens_data["material_l4"]=surflist[4]["GLAS"][0]
        
        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        r3=surflist[3]["DIAM"][0]
        r4=surflist[4]["DIAM"][0]
        r5=surflist[5]["DIAM"][0]
    
        assert r0==r1 and r1== r2 and r2 ==r3 and r3==r4 and r4==r5
        
        lens_data["radius"] = r0
        lens_data["glass_catalogs"] = gcat #This is not used for the moment

        lens_data["type"] = "DoubletPair"
        pyot[key]=lens_data
        return pyot
        
    elif ns==7: #Doublet pair with stop in the middle, stop ignored
        lens_data["curvature_s1"]=surflist[0]["CURV"][0]
        lens_data["curvature_s2"]=surflist[1]["CURV"][0]
        lens_data["curvature_s3"]=surflist[2]["CURV"][0]
        lens_data["curvature_s4"]=surflist[4]["CURV"][0]
        lens_data["curvature_s5"]=surflist[5]["CURV"][0]
        lens_data["curvature_s6"]=surflist[6]["CURV"][0]
        
        lens_data["thickness_l1"]=surflist[0]["DISZ"][0]
        lens_data["thickness_l2"]=surflist[1]["DISZ"][0]
        lens_data["air_gap"]=surflist[2]["DISZ"][0]+surflist[3]["DISZ"][0]
        lens_data["thickness_l3"]=surflist[4]["DISZ"][0]
        lens_data["thickness_l4"]=surflist[5]["DISZ"][0]
        
        lens_data["material_l1"]=surflist[0]["GLAS"][0]
        lens_data["material_l2"]=surflist[1]["GLAS"][0]
        lens_data["material_l3"]=surflist[4]["GLAS"][0]
        lens_data["material_l4"]=surflist[5]["GLAS"][0]
        
        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        r3=surflist[4]["DIAM"][0]
        r4=surflist[5]["DIAM"][0]
        r5=surflist[1]["DIAM"][0]

        assert r0==r1 and r1== r2 and r2 ==r3 and r3==r4 and r4==r5
        
        lens_data["radius"] = r0
        lens_data["glass_catalogs"] = gcat #This is not used for the moment

        lens_data["type"] = "DoubletPair"
        pyot[key]=lens_data
        return pyot
  
    else:
        # Print elements present in the library that were not processed or explicitelly excluded
        print("***",key)
        print(description)
        print("Surfaces=",len(surflist))
        for i in surflist:
            print("*", i)

