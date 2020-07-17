#====================================================#
#重心からもっとも遠くにある原子上位50個の原子名とその距離                  
#====================================================#
#!! process !!
# 1. extract atom_coord from PDB file
# 2. calc G(center of gravity)
# 2. calc each distance from atom to G
# 3. apply to bubble sort
# 4. print top50 value of distance and its atom_name
#-----------------------------------------------------
#実際に重心-原子間の距離の計算をし、昇順にソートしているのは
#37-99行目になります.
#
#それ以外（原子座標の抽出、重心の計算）は、以前の
#重心および慣性半径算出のプログラムコードと
#同じものを使用しています.
#-----------------------------------------------------

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


#funtion to print top50 value of distance from atom to G(center_of_gravity)
#and its atom name
#TOP50の原子-重心間距離の結果出力用の関数
def print_top50_distance_atomToG(distances_atomToG_bubbleSorted):
    for distance_atomToG in distances_atomToG_bubbleSorted:
        print(distance_atomToG)
    
    return  #(atom_number, atom_name, value of distance(Å))


#function to get top50 value of distance from atom to G, and its atom name
#重心からもっとも遠くにある原子上位50個の原子名とその距離を得るための関数
def top50_distance_atomToG(atoms_coord):
    distances_atomToG = calc_distances_atomToG(atoms_coord)
    distances_atomToG_bubbleSorted = bubblesort(distances_atomToG)

    return distances_atomToG_bubbleSorted[0:49] #(atom_number, atom_name, value of distance(Å))


#function for bubble sort（ascending order）
#得られた配列（距離）を昇順にソートするための関数
def bubblesort(distances_atomToG):
    for i in range(len(distances_atomToG)):
        for j in range(len(distances_atomToG)-1, i, -1):
            if distances_atomToG[j][2] > distances_atomToG[j-1][2]:
                distances_atomToG[j],distances_atomToG[j-1] =  distances_atomToG[j-1], distances_atomToG[j]
    
    distances_atomToG_bubbleSorted = distances_atomToG
    
    return distances_atomToG_bubbleSorted   #(atom_number, atom_name, value of distance(Å))


#function to calc each distance from atom to G
#各原子と重心との間の距離を計算するための関数
def calc_distances_atomToG(atoms_coord):
    G = calc_centerOfGravity(atoms_coord)
    Xg = G[0]   # x_coord of G
    Yg = G[1]   # y_coord of G
    Zg = G[2]   # z_coord of G

    distances_atomToG = []   #distanced beetween G and atom for every atom

    for atom_coord in atoms_coord:

        #atom_coord
        atom_x = atom_coord[2]
        atom_y = atom_coord[3]
        atom_z = atom_coord[4]

        #calc desitance_squared
        distance_squared_x = (atom_x - Xg)**2
        distance_squared_y = (atom_y - Yg)**2
        distance_squared_z = (atom_z - Zg)**2

        distance_squared = (distance_squared_y + distance_squared_y + distance_squared_z)
        distance = (distance_squared)**0.5 

        #append to the list
        distance_atomToG = [atom_coord[0], atom_coord[1], distance]
        distances_atomToG.append(distance_atomToG)
    
    return distances_atomToG    #(atom_number, atom_name, value of distance(Å))


#function for calculating the center_of_gravity(G)
#center_of_gravity(Xg) =  (sum of each atom's x_coord) / (number of atoms + 1)
#G(x, y, z) = (Xg, Yg, Zg) 
#重心計算用の関数
def calc_centerOfGravity(atoms_coord):
    sum_coord = calc_sumCoord(atoms_coord)

    sum_x = sum_coord[0]
    sum_y = sum_coord[1]
    sum_z = sum_coord[2]
    count_atoms = len(atoms_coord)    #number of atoms

    Xg = sum_x/(count_atoms+1)  
    Yg = sum_y/(count_atoms+1)  
    Zg = sum_z/(count_atoms+1)  

    G = [Xg, Yg, Zg] 
    return G    #(Xg, Yg, Zg), type(Xg, Yg, Zg) = float



#declares a function for calculating sums of each atom_coords(x, y, z) 
#重心計算のための、各座標の和の計算用関数
def calc_sumCoord(atoms_coord):
    sum_x = 0   #var for sum of each atom's x_coord
    sum_y = 0   #var for sum of each atom's y_coord
    sum_z = 0   #var for sum of each atom's z_coord

    #calculates sums of each atom_coords
    #重心計算用に、各座標の和の算出
    for atom_coord in atoms_coord:
        #sum of x_coord
        atom_x = atom_coord[2]
        sum_x += atom_x

        #sum of y_coord
        atom_y = atom_coord[3]
        sum_y += atom_y

        #sum of z_coord
        atom_z = atom_coord[4]
        sum_z += atom_z

    return [sum_x, sum_y, sum_z]    #(sum of x_coord、sum of y_coord、sum of z_coord) type = float



# =========================
#          Execution
# =========================

#extract atomsCoord from PDB
pdb_1gri = '1gri.pdb'
atoms_coord_1gri = extract_atomsCoord(pdb_1gri)

#calc distance _atomToG /get bubbleSorted list
top50_distance_atomToG_1gri = top50_distance_atomToG(atoms_coord_1gri)

#print Top50 distace_atomToG_1gri
print_top50_distance_atomToG(top50_distance_atomToG_1gri)

