"""
Created on March 14, 2019 by Andrew Abi-Mansour

The code in this file belongs to James Davidson and is discussed here::

	http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01185.html
	http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01162.html
	http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01900.html  

This is the::

	     _   _   _ _____ ___    __  __    _    ____ _____ ___ _   _ ___ 
	    / \ | | | |_   _/ _ \  |  \/  |  / \  |  _ \_   _|_ _| \ | |_ _|
	   / _ \| | | | | || | | | | |\/| | / _ \ | |_) || |  | ||  \| || | 
	  / ___ \ |_| | | || |_| | | |  | |/ ___ \|  _ < | |  | || |\  || | 
	 /_/   \_\___/  |_| \___/  |_|  |_/_/   \_\_| \_\|_| |___|_| \_|___|                                                            
                                                                 
Tool for automatic MARTINI mapping and parametrization of small organic molecules

Developers::

	Tristan BEREAU (bereau at mpip-mainz.mpg.de)
	Kiran Kanekal (kanekal at mpip-mainz.mpg.de)
	Andrew Abi-Mansour (andrew.gaam at gmail.com)

AUTO_MARTINI is open-source, distributed under the terms of the GNU Public
License, version 2 or later. It is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
received a copy of the GNU General Public License along with PyGran.
If not, see http://www.gnu.org/licenses . See also top-level README
and LICENSE files.
"""


mb = [ """MolPort-000-036-699
  Marvin  05240819062D          
 39 44  0  0  0  0            999 V2000
   -0.9093   -1.0027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1895    0.3923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4328   -1.6783    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6446   -0.2242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9279    1.1739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6970   -1.2486    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9248   -2.3416    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7064   -2.0738    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9991    0.2335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4698    1.7936    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1619   -0.0623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.7811    0.2584    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5409    0.8532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5876    0.4204    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2762    1.6316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3768   -1.8309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1183    1.3358    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1356   -0.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7100   -0.6819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4266    0.7193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5165   -0.5231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2362    0.8781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8523    1.2051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3505    0.6881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8243    2.2513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9452   -0.0342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6619    1.3639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8954    1.3078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6338    2.0925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2068    0.7473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8678   -0.9777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9155   -1.2082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6508   -2.6094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7019    1.1459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9266    2.1454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9635    0.3643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7250   -1.3577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4604   -2.7619    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9991   -2.1361    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  4  1  0  0  0  0
  3  1  4  0  0  0  0
  4  1  1  0  0  0  0
  5  2  4  0  0  0  0
  6  1  4  0  0  0  0
  7  3  4  0  0  0  0
  8  6  4  0  0  0  0
  9  2  4  0  0  0  0
 10  5  4  0  0  0  0
 11  4  1  0  0  0  0
 12 22  1  0  0  0  0
 13  9  4  0  0  0  0
 14 12  1  0  0  0  0
 15 10  4  0  0  0  0
 16  3  1  0  0  0  0
 17  5  2  0  0  0  0
 18 14  4  0  0  0  0
 19 11  1  0  0  0  0
 20 11  1  0  0  0  0
 21 19  1  0  0  0  0
 22 20  1  0  0  0  0
 23 14  4  0  0  0  0
 24 13  4  0  0  0  0
 25 15  4  0  0  0  0
 26 18  4  0  0  0  0
 27 23  4  0  0  0  0
 28 24  4  0  0  0  0
 29 25  4  0  0  0  0
 30 27  4  0  0  0  0
 31 18  1  0  0  0  0
 32 16  1  0  0  0  0
 33 16  1  0  0  0  0
 34 28  1  0  0  0  0
 35 27  1  0  0  0  0
 36 34  1  0  0  0  0
 37 32  1  0  0  0  0
 38 33  1  0  0  0  0
 39 37  1  0  0  0  0
  8  7  4  0  0  0  0
 12 21  1  0  0  0  0
 39 38  1  0  0  0  0
 15 13  4  0  0  0  0
 29 28  4  0  0  0  0
 30 26  4  0  0  0  0
M  END""",
 
"""MolPort-000-005-089
  Marvin  05210810082D          
 19 21  0  0  0  0            999 V2000
    1.4290    3.4338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4290    2.6088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434    2.1963    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8579    2.6088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8579    3.4338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434    3.8463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6425    2.3538    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1274    3.0213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6425    3.6888    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.9525    3.0213    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8974    1.5693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3454    0.9562    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6004    0.1715    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4073    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.9594    0.6131    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7044    1.3977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    3.8463    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    0.7145    4.6714    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.4338    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
  4  7  4  0  0  0  0
  7  8  4  0  0  0  0
  8  9  4  0  0  0  0
  9  5  4  0  0  0  0
  8 10  2  0  0  0  0
  1 17  1  0  0  0  0
  7 11  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  1  0  0  0  0
 13 14  1  0  0  0  0
 14 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 11  1  0  0  0  0
 17 18  2  0  0  0  0
 17 19  1  0  0  0  0
M  CHG  2  17   1  19  -1
M  STY  1   1 DAT
M  SAL   1  1   9
M  SDT   1 MRV_IMPLICIT_H                                        
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 IMPL_H1
M  END
> <PUBCHEM_EXT_DATASOURCE_REGID>
MolPort-000-005-089
> <PUBCHEM_EXT_SUBSTANCE_URL>
http://www.molport.com/buy-chemicals/molecule-link/MolPort-000-005-089
> <PUBCHEM_EXT_DATASOURCE_URL>
http://www.molport.com
$$$$""" ]


from rdkit import Chem
from rdkit.Chem import AllChem

def FragIndicesToMol(oMol,indices):
    em = Chem.EditableMol(Chem.Mol())

    newIndices={}
    for i,idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx]=i

    for i,idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx()==idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx<idx:
                continue
            em.AddBond(newIndices[idx],newIndices[oidx],bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap=newIndices
    return res

def _recursivelyModifyNs(mol,matches,indices=None):
    if indices is None:
        indices=[]
    res=None
    while len(matches) and res is None:
        tIndices=indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol.ToBinary())
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res,indices = _recursivelyModifyNs(nm,matches,indices=tIndices)
        else:
            indices=tIndices
            res=cp
    return res,indices

def AdjustAromaticNs(m,nitrogenPattern='[n&D2&H0;r5,r6]'):
    """
       default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
       to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    plsFix=set()
    for a,b in linkers:
        em.RemoveBond(a,b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at=nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum()==7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [FragIndicesToMol(nm,x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok=True
    for i,frag in enumerate(frags):
        cp = Chem.Mol(frag.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres,indices=_recursivelyModifyNs(frag,matches)
            if not lres:
                #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok=False
                break
            else:
                revMap={}
                for k,v in frag._idxMap.items():
                    revMap[v]=k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m
