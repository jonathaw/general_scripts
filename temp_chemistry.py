from rdkit.Chem import AllChem
from rdkit.Chem import Draw

smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
mol = AllChem.MolFromSmiles(smiles)

# technically this step isn't required since the drawing code
# will automatically add a 2D conformation to a molecule that has
# no conformation information, I'm including it to show how to
# generate 2D coords with the RDKit:
AllChem.Compute2DCoords(mol)
print(mol)
Draw.MolToFile(mol,"caffeine.png",size=(200,250))