# ---------------池谷先生からのコメント（done = 修正済み）-------------------

# 1.(done) 関数では、変数は必ず引数で渡す。Global変数は使わない。

# (done)2. コメントは英語で書く。
#    Scienceは必ず世界の人たちと共同研究し、競争するので、
#     廣瀬くんのプログラムが外国の人が見るかもしれません。
#    大変かもしれませんが、これも練習です。英語で苦労して損はないはずです。
#    英語で書きましょう。

#(done)3. 1GRIは、ダイマーになっています。モノマーの値が知りたいです。

#(done)4. すでに書いたように、BioPythonは使わない。

#(done) 5. for loopで、iは暗黙的にindex (integer)なので、座標として使わないほう
# が良いかも。本質ではないですが、なるべく慣習に従うべきです。

# (done)6. 関数は始めにまとめて定義した方が良いと思いますけどね。可読性やオブジ
# ェクト指向という観点から。
# -------------------------------------------------------------------------
# 課題1
# ・GB1のpdbデータから原子座標を読み込み
# ・GB1の重心の計算
# ・GB1の慣性半径の計算
# 
# 課題2
# ・重心からもっとも遠くにある原子上位50個の原子名とその距離の表示（遠い順）
# -------------------------------------------------------------------------
# this program is for calculating the center of gravity(G) and radius of Gyration(Rg).
#
# 1.function "calc_centerOfGravity( PDB_file(ATOM) ) "
# you can get the value of G with executing the function "calc_centerOfGravity( PDB_file(ATOM) ) " .
#
# 2.function "calc_RadiusOfGyration( PDB_file(ATOM) )  "
# you can get the value of Rg with the function "calc_RadiusOfGyration( PDB_file(ATOM) ) " as well.
#
#【Remenber】
# Remenber that you need to execute the function "extract_atomsCoord( PDB_file )" to import atom's coord(x,y,z) 
# of your protein before using those two functions above.
# -------------------------------------------------------------------------- 

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


#function for calculating the center_of_gravity(G)
#center_of_gravity(Xg) =  (sum of each atom's x_coord) / (number of atoms + 1)
#G(x, y, z) = (Xg, Yg, Zg) 
#重心計算用の関数の定義
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


#declares a funciton for calculating the Radius_of_Gyration
#(Radius_of_Gyration)^2 = Σ(each atom's coord - G)^2 / number of atoms
#慣性半径算出用の関数の定義
def calc_RadiusOfGyration(atoms_coord):
    sum_radius_squared = calc_sumRadiusSquared(atoms_coord)
    count_atoms = len(atoms_coord)    #number of atoms

    Rg_squared = sum_radius_squared / count_atoms   #(Radius_of_Gyration)^2
    Rg = Rg_squared**0.5 #Radius_of_Gyration

    return Rg


#declares a function for calculating sums of each atom_coords(x, y, z) 
#重心計算のための、各座標の和の計算用関数の定義
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


#declares a function for calculating the distance between each atom's coord and G(: center of gravity)
#(Radius_of_Gyration)^2 = Σ(each atom's coord - G)^2 / number of atoms
#各原子の座標と、重心の位置との距離の差を求める関数の定義
#慣性半径^2＝Σ（各原子の座標-重心)^2 / 全原子数
def calc_sumRadiusSquared(atoms_coord):
    G = calc_centerOfGravity(atoms_coord)
    Xg = G[0]   # x_coord of G
    Yg = G[1]   # y_coord of G
    Zg = G[2]   # z_coord of G

    sum_radius_squared = 0  #sum of squared radius (: each atom's coord - G) 

    for atom_coord in atoms_coord:
        atoms_x = atom_coord[2]  #x_coord
        atoms_y = atom_coord[3]  #y_coord
        atoms_z = atom_coord[4]  #z_coord

        # calculates the radius 
        radius_squared_x = (atoms_x - Xg)**2    #x_component of radius
        radius_squared_y = (atoms_y - Yg)**2    #y_component of radius
        radius_squared_z = (atoms_z - Zg)**2    #z_component of radius

        radius_squared = (radius_squared_x + radius_squared_y + radius_squared_z)
        
        #add to the var sum_radius_squared
        sum_radius_squared += radius_squared
    
    return sum_radius_squared 



#Execution
pdb = '1gri.pdb'    #pdb data of the your protein

atoms_coord_1gri = extract_atomsCoord(pdb)    #extract atoms_coord from pdb data
center_of_gravity_1gri = calc_centerOfGravity(atoms_coord_1gri)  #calculates the center fo gravity
radius_of_gyration_1gri = calc_RadiusOfGyration(atoms_coord_1gri) #calculates the radius of gyration

# print (G and Rg)
print("Center of Gravty(x, y, z (Å))",":","G", "(", center_of_gravity_1gri[0],",", center_of_gravity_1gri[1],",", center_of_gravity_1gri[2], ")")
print("Radius of Gyration（Å）",":","Rg", radius_of_gyration_1gri)