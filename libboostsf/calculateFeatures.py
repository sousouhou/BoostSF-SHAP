
import os 
from ligandm import Ligand
from receptorm import Receptor


def cal36features(fileReceptor, fileLigand) :
    rec = Receptor(fileReceptor)
    lig = Ligand(fileLigand)
    
    reVal = [0.0 for i in range(0, 36)]
    
    for la in lig.allatoms :
        if la.numType == 1000:
            continue
            
        for ra in rec.allatoms :
            if ra.numType == 1000 :
                continue
        
            d2 = (la.coord_x - ra.coord_x)**2 +  (la.coord_y - ra.coord_y)**2 +  (la.coord_z - ra.coord_z)**2
            
            if d2 >= 12*12 :  # RF-Score cutoff 12A
                continue

            reVal[la.numType*4 + ra.numType]  += 1
    
    return reVal


def calVina5terms(fileReceptor, fileLigand)  :
    
    reVal = [0.0 for i in range(0, 5)]

    dirname, filename = os.path.split(os.path.abspath( __file__))
    vinaPath = dirname + r"\vina.exe"
    
    #r".\vina.exe --receptor 1c5z_protein.pdbqt  --ligand  1c5z_ligand.pdbqt  --score_only "
    cmd = "%s  --receptor  %s  --ligand  %s --score_only > ./tempVinaout.txt"%(vinaPath, fileReceptor, fileLigand)
    # print(cmd)
    os.system(cmd)
    
    pf = open("./tempVinaout.txt", "r")
    strlines = pf.readlines()
    pf.close()

    os.remove("./tempVinaout.txt")
    
    for line in strlines:
        if "gauss 1" in line :      # "  gauss 1     : 85.68172 "
            strlist = line.split(":")
            reVal[0] = float( strlist[1].rstrip() )
            
        if "gauss 2" in line :     
            strlist = line.split(":")
            reVal[1] = float( strlist[1].rstrip() )
            
        if "repulsion" in line :     
            strlist = line.split(":")
            reVal[2] = float( strlist[1].rstrip() )
            
        if "hydrophobic" in line :     
            strlist = line.split(":")
            reVal[3] = float( strlist[1].rstrip() )
            
        if "Hydrogen" in line :     
            strlist = line.split(":")
            reVal[4] = float( strlist[1].rstrip() )
            
    return reVal  
    

def cal41features(fileReceptor, fileLigand)  :
    res1 = cal36features(fileReceptor, fileLigand)
    res2 = calVina5terms(fileReceptor, fileLigand)
    
    return res1 + res2
    

if __name__ == "__main__" :        
    
    # ----- debug 1 -----
    # Feas = cal36features("../test/1a30_protein.pdbqt", "./1a30_ligand.pdbqt")
    # print("1a30: ")
    # for t in Feas :
        # print("%.4f,"%t, end="" )

    # ----- debug 2 -----
    # Feas = calVina5terms("../test/1a30_protein.pdbqt", "./1a30_ligand.pdbqt")
    # print("1a30: ")
    # for t in Feas :
        # print("%.4f,"%t, end="" )

    # ----- debug 3 -----   
    Feas = cal41features("../test/1a30_protein.pdbqt", "../test/1a30_ligand.pdbqt")
    for t in Feas :
        print("%.4f,"%t, end="" )   
    