#!/usr/local/bin/python3.13
#La ligne ci-dessus est nécessaire pour exécuter le fichier sur Linux
#Nom / prénom / groupe
#MECHOUB Gaya G1
#BOUZIDI Massyl G2
#CHERIEF Lamine G5

#******PARTIE ENTREES*********

n = int(input("Quel est le nombre des variables de décision (Xi) ? "))
m = int(input("Quel est le nombre des contraintes ? "))

# Saisie des coefficients de la fonction objectif
FctObjectif = []
print("Entrez les coefficients de la fonction objectif Z : ")
for i in range(n):
    FctObjectif.append(float(input(f"Coefficient de X{i+1} : ")))

A = []  # Matrice des coefficients des contraintes
B = []  # Vecteur des seconds membres
TC = []  # Liste des types de contraintes

for i in range(m):
    ligne = []  #Saisie des coefficients des contraintes
    print(f"Pour la contrainte numéro {i+1}")
    for j in range(n):
        ligne.append(float(input(f"Coefficient de X{j+1} : ")))
    A.append(ligne)
    B.append(float(input(f"Second membre b{i+1} : ")))

    # Saisie des types de contrainte
    typeContrainte = " "    #On estime ici qu'il y a un risque que l'utilisateur saisisse un opérateur non valide ou autre chose que demandé
    while not (typeContrainte == "<=" or typeContrainte == ">=" or typeContrainte == "="):
        try:
            typeContrainte = input(f"Type de contrainte (<=, >= ou =) : ")
            if not (typeContrainte == "<=" or typeContrainte == ">=" or typeContrainte == "="):
                raise ValueError
        except ValueError:
            print("Vous n'avez pas saisi un type valide.")
    TC.append(typeContrainte)

# Saisie du type d'optimisation
typeOptimisation = input("Type d'optimisation (max ou min) : ") #Ici également on s'asssure que l'utilisateur n'entre pas autre chose...
while not (typeOptimisation == "max" or typeOptimisation == "min"):
    print("Vous n'avez pas saisi un type valide.")
    typeOptimisation = input("Type d'optimisation (max ou min) : ")

contraintes_opr = []
for tc in TC:
    if tc == "<=":
        contraintes_opr.append(-1)
    elif tc == "=":
        contraintes_opr.append(0)
    elif tc == ">=":
        contraintes_opr.append(1)

#******PARTIE CORPS DU CODE******

type = typeOptimisation.upper()
obj_func_coeff = FctObjectif
contraintes_coeff = A
contraintes_sm = B

EMPTY_CELL = "*"    #Cellule vide

def get_super(x):  #Fonction pour obtenir les exposants
    return str(x)

def rechercher(vect, maxi=True, positif=True):  #Fonction pour rechercher le pivot
    new_vect = [i for i in vect if i != None and i < 0]
    if positif:
        new_vect = [i for i in vect if i != None and i > 0]
    if new_vect == []:
        return None
    if maxi:
        return vect.index(max(new_vect))
    return vect.index(min(new_vect))

def rechercherCP(vect): #Fonction pour rechercher la colonne de pivot
    if type == "MAX":
        return rechercher(vect, maxi=True, positif=True)
    if degenere:
        return rechercher(vect, maxi=True, positif=False)
    return rechercher(vect, maxi=False, positif=False)

def clone_matrice(matrice): #Fonction pour faire une copie conforme d'une matrice
    return [[i for i in vect] for vect in matrice]

def eliminer_fraction(nbr): #Fonction pour éliminer les fractions
    try :
        return round(nbr) if round(nbr) == nbr else round(nbr, 2)
    except :
        return nbr

def format_equation(vect, op):  #Fonction pour formater les équations c'est-à-dire les contraintes et la fonction objectif
    equa = ""
    is_first = True
    for i in range(len(vect)):
        if vect[i] != 0:
            sign = " - " if vect[i] < 0 else " + "
            coeff = "" if vect[i] in [-1, 1] else str(abs(eliminer_fraction(vect[i]))) + "*"
            sign = (sign.strip() if sign == " - " else "") if is_first else sign
            equa += "{}{}x{}".format(sign, coeff, get_super(i + 1))
            is_first = False
    
    return ('0' if equa == '' else equa) + (" ≥" if op > 0 else " ≤" if op < 0 else " =")

def format_matrice():   
    ADD_SM_CP_COLOMN = True 
    COL = col
    new_mat = clone_matrice(simplex_matrice)
    
    if big_M :
        for i in range(COL) :
            nbr = eliminer_fraction(new_mat[-1][i])
            if nbr == None :
                continue
            Mnbr = eliminer_fraction(simplex_Mobj_func[i])
            sign = '' if Mnbr > 0 and nbr == 0 else '-' if Mnbr < 0 else '+'
            Mnbr = '' if Mnbr == 0 else sign + (str(abs(Mnbr)) if Mnbr not in [1, -1] else '') + 'M'
            new_mat[-1][i] = (str(nbr) if nbr != 0 else '' if Mnbr != '' else '0') + Mnbr

    if ADD_SM_CP_COLOMN :   #Ajout de la colonne SM/CP
        for i in range(row - 1) :
            new_mat[i] += [eliminer_fraction(contraintes_sm_cp[i])]
        new_mat[-1].append(EMPTY_CELL)
        COL += 1
    for i in range(row):        
        for j in range(COL):
            new_mat[i][j] = str(eliminer_fraction(new_mat[i][j])) if new_mat[i][j] != None else EMPTY_CELL
    header = ["x" + get_super(i+1) for i in range(n)] + ["e" + get_super(i+1) for i in range(ecart_nbr)]
    header += ["t" + get_super(i+1) for i in range(art_nbr)]
    header += ["SM", "SM/CP"] if ADD_SM_CP_COLOMN else ["SM"]
    new_mat.insert(0, header)

    # Ajout variables de base
    for i in range(row - 1) :
        new_mat[i + 1].insert(0, header[vbs[i]])
    new_mat[0].insert(0, "Base")
    new_mat[-1].insert(0, "-z")
    return new_mat

def afficherMatrice():
    COL = col + 2
    ROW = row + 1
    matrice = format_matrice()
    col_max_len = [0]*(COL)

    for i in range(COL) :
        col_max_len[i] = max(len(row[i]) for row in matrice)

    print()
    for i in range(ROW) :
        print("|", end='')
        for j in range(COL) :
            elem = matrice[i][j]
            print('', elem, ' '*(col_max_len[j] - len(elem)) + '|', end='')
        print()


n = max([len(vect) for vect in contraintes_coeff + [obj_func_coeff]]) 
####################### Corriger le nombre des variables  ###########################
contraintes_coeff = [constraint + [0]*(n - len(constraint)) for constraint in contraintes_coeff]
obj_func_coeff += [0]*(n - len(obj_func_coeff))

# les seconds membres doivent être positifs
is_contraintes_changed = False
for i in range(m):
    if contraintes_sm[i] < 0:
        is_contraintes_changed = True
        contraintes_sm[i] *= -1
        contraintes_coeff[i] = [-i for i in contraintes_coeff[i]]
        contraintes_opr[i] *= -1
# Vérifier si nous pouvons résoudre le PL avec la méthode du simplex
big_M = any([False if opr == -1 else True for opr in contraintes_opr])

# Calculer combien de variables artificielles et d'écart nous ajouterons
art_nbr = 0
ecart_nbr = m
for opr in contraintes_opr :
    art_nbr += 1 if opr != -1 else 0
    ecart_nbr -= 1 if opr == 0 else 0

# Générer la matrice des variables d'éxcedent et d'écart et des variables artificielles
ecart_coeff = [[0 for i in range(ecart_nbr)] for i in range(m)]
art_coeff = [[0 for i in range(art_nbr)] for i in range(m)]

tmp_e, tmp_a = 0, 0
for i in range(m) :
    if contraintes_opr[i] != 0 :
        ecart_coeff[i][tmp_e] = - contraintes_opr[i]
        tmp_e += 1
    if contraintes_opr[i] != -1: 
        art_coeff[i][tmp_a] = 1
        tmp_a += 1

obj_func_ecoeff = [0]*ecart_nbr
obj_func_acoeff = [0]*art_nbr
# Calculer la nouvelle fonction objectif de Big M
Mobj_func_coeff  = [0]*n
Mobj_func_ecoeff = [0]*ecart_nbr
Mobj_func_acoeff  = [0]*art_nbr
Mobj_func_sm = 0

for i in range(m) :
    if contraintes_opr[i] == -1:
        continue
    for j in range(n):
        Mobj_func_coeff [j] += -contraintes_coeff[i][j]
    
    for j in range(ecart_nbr):
        Mobj_func_ecoeff[j] += -ecart_coeff[i][j]
    
    Mobj_func_sm += contraintes_sm[i]

if type == "MAX":
    Mobj_func_coeff  = [-i for i in Mobj_func_coeff ]
    Mobj_func_ecoeff = [-i for i in Mobj_func_ecoeff]
    Mobj_func_sm *= -1

Mobj_func_sm *= -1
# déclarer un tableau des indices des variables de base pour suivre les changements lors de l'application de l'algorithme
vbs = [None] * m
tmp_e, tmp_t = 0, 0
for i in range(m):
    if contraintes_opr[i] == -1:
        vbs[i] = n + tmp_e
        tmp_e += 1
        continue
    else :
        vbs[i] = n + ecart_nbr + tmp_t
        tmp_t += 1

#Génératio de la matrice du tableau simplexe
simplex_matrice = clone_matrice(contraintes_coeff)
for i in range(m) :
    simplex_matrice[i] += ecart_coeff[i] + art_coeff[i] + [contraintes_sm[i]]
simplex_matrice += [[i for i in obj_func_coeff]]
simplex_matrice[-1] += obj_func_ecoeff + obj_func_acoeff +  [0]

simplex_Mobj_func = Mobj_func_coeff  + Mobj_func_ecoeff + Mobj_func_acoeff  + [Mobj_func_sm]

col = n + ecart_nbr + art_nbr + 1
row = m + 1

degenere = False
k = 0

contraintes_sm_cp  = [0]*m
obj_func_sm_cp = 0


#********PARTIE AFFICHAGE********

# Afficher le tableau initial (k = 0) seulement
print("\nK =", k, " (Tableau simplexe initial du problème) : ") 
afficherMatrice()   #Affichage de la matrice

while True :
    # Trouver la colonne de pivot
    if big_M :
        cp = rechercherCP(simplex_Mobj_func[: n + ecart_nbr])
    if not(big_M) or (big_M and (cp == None and all(not i for i in simplex_Mobj_func[: n + ecart_nbr]))):
        cp = rechercherCP(simplex_matrice[-1][: n + ecart_nbr])
    
    # Test de critere d'arret
    if cp == None:
        break

    # Calculer SM/CP
    for i in range(m):
        if simplex_matrice[i][cp] == 0 :
            contraintes_sm_cp[i] = None
            continue
        contraintes_sm_cp[i] = simplex_matrice[i][-1] / simplex_matrice[i][cp]
    
    # Trouver la ligne de pivot
    lp = rechercher(contraintes_sm_cp, maxi=False, positif=True)

    if lp == None:
        print("(!) La solution est infinie")
        exit(1)
    
    # Trouver le pivot
    pivot = simplex_matrice[lp][cp]

    # Nouvelle variable de base
    vbs[lp] = cp
    k += 1

    # Application de l'algorithme Gausséen
    simplex_matrice[lp] = [nbr / pivot if nbr != None else None for nbr in simplex_matrice[lp]]

    # Calculer la nouvelle matrice
    old_mat = clone_matrice(simplex_matrice)
    for i in range(row):
        if i == lp:
            continue
        for j in range(col):
            try :
                simplex_matrice[i][j] = old_mat[i][j] - old_mat[i][cp] * old_mat[lp][j]
            except :
                pass

    if big_M:
        old_Mobj_func = simplex_Mobj_func.copy()
        for i in range(col):
            try :
                simplex_Mobj_func[i] = old_Mobj_func[i] - old_Mobj_func[cp] * old_mat[lp][i]
            except :
                pass

        # Éliminer la variable sortante
        for i in range(n + ecart_nbr, n + ecart_nbr + art_nbr):
            if i not in vbs:
                for j in range(row):
                    simplex_matrice[j][i] = None
                simplex_Mobj_func[i] = None
    
    # Verifier si la solution est dégénéré
    for i in range(row - 1):
        if (simplex_matrice[i][-1] == 0):
            print("\nLa solution est dégénéré")
            afficherMatrice()

            if degenere:
                print("Nous ne pouvons pas appliquer la règle de Bland pour la deuxième fois ...\n")
                exit(1)

            input("\nCliquez sur ENTRÉE pour appliquer la règle de Bland : ")
            degenere = True
            simplex_matrice = clone_matrice(old_mat)
            if big_M:
                simplex_Mobj_func = old_Mobj_func.copy()
            break

optimal_sol = [0] * n
for i in range(m):
    if vbs[i] < n:
        optimal_sol[vbs[i]] = eliminer_fraction(simplex_matrice[i][-1])

print("\n\nLa solution optimale :\nx* = (", end="")
for i in range(n):
    print(optimal_sol[i], end="" if i == n - 1 else ", ")
print(")\n")

print("L'optimum :\nz* =", -eliminer_fraction(simplex_matrice[-1][-1]), "\n")
