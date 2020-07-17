# ----------------- Todo -------------------#
# .PDBファイルから原子座標データだけを取り出す方法
# 1.PDBファイルを読み込む（ファイル名、r）
# 2.スペースごとに区切り、リストに追加していく
# 3.もし配列の何個目の要素がATOMなら、ATOMリストに（原子番号, x, y, z,）を追加
# 注意：PDBデータでは座標の型はstrになってるため、floatに変換して格納する
# ------------------------------------------#

pdb = '1gri.pdb'

#function to store the atom_coord data from the PDB file
#PDBデータから原子座標データを抽出
def extract_atomsCoord(PDB_file):
    atoms_coord = []
    
    #read PDB_file
    with open(PDB_file, 'r') as f: 
        for row in f:
            row_pdb = row.split()

            #separate protein to monomer
            if row_pdb[0] == "TER":     
                break
            
            #extract only atoms data(x, y, z) from PDB_file
            if row_pdb[0] == "ATOM":
                #atom_coord = (atom_number, atomic_name, x, y, z)
                atom_coord = [row_pdb[1], row_pdb[2], float(row_pdb[6]), float(row_pdb[7]), float(row_pdb[8])]   
                atoms_coord.append(atom_coord)
    
    return atoms_coord  #(atom_number, atomic_name,  x, y, z), type(x,y,z) = float