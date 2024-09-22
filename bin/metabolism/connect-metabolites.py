"""Use RDKit to clean reactants and products, and identify shared and similar metabolites between pairs."""

from typing import Iterable, List, Optional, Union
from functools import partial
from itertools import product
import sys

from carabiner import print_err
from carabiner.cast import cast
import pandas as pd
import numpy as np
from rdkit import RDLogger    
from rdkit.Chem import AddHs, rdFingerprintGenerator, MolFromSmiles, MolToInchi, MolToSmiles
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.DataStructs import TanimotoSimilarity
from tqdm.auto import tqdm

RDLogger.DisableLog('rdApp.*')     

removers = (
    lambda x: str.replace(x, ".[H+].", "").replace(".[H+]", "").replace("[H+].", ""),
    SaltRemover(defnData="[H,Na,K,Mg,Ca,Mn,Fe,F,Cl,Br,O,S]").StripMol,
    partial(SaltRemover().StripMol, dontRemoveEverything=True),
    AddHs,
)

def _clean_smiles_list(smiles=str) -> str:
    mol = MolFromSmiles(smiles)
    for remover in removers:
        mol = remover(mol)
    return MolToSmiles(mol)


def clean_metabolites(table: pd.DataFrame, 
                      cols=Iterable[str]) -> pd.DataFrame:
    return table.assign(**{col: table[col].apply(_clean_smiles_list) for col in cols})


def _shared(a_b) -> bool:
    a, b = (map(MolFromSmiles, a_or_b) for a_or_b in a_b)
    b = set(map(MolToInchi, b))
    return any(MolToInchi(item) in b for item in a)


def _have_shared_chemical(df: pd.DataFrame, a: str, b: str) -> pd.Series:
    return df[[a, b]].progress_apply(_shared, axis=1)


def _have_shared_chemical_simple(df: pd.DataFrame, col: str) -> pd.Series:
    return df[col] == 1.

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _fingerprinter(a):
    return mfpgen.GetFingerprint(MolFromSmiles(a))


def _tanimoto(a_b):
    a, b = a_b
    return max(TanimotoSimilarity(*map(_fingerprinter, items)) for items in product(*a_b))


def _max_similarity(df: pd.DataFrame, a: str, b: str):
    return df[[a, b]].progress_apply(_tanimoto, axis=1)


def _smiles_splitter(df: pd.DataFrame, col: str, char: str = ".") -> List[str]:
    return df[col].str.split(char)


def connect_metabolites(table: pd.DataFrame,
                        cols: Iterable[str],
                        id_col: Union[str, Iterable[str]],
                        table2: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    
    if table2 is None:
        table2 = table.copy()
    if isinstance(id_col, str):
        id_col = cast(id_col, to=list)
    cols = cast(cols, to=list)

    all_cols = id_col + cols
    all_cols2 = [f"{col}_2" for col in all_cols]
    cols2 = [f"{col}_2" for col in cols]

    table = table[all_cols].assign(**{f"{col}_split": partial(_smiles_splitter, col=col) for col in cols})
    table2 = (table2[all_cols]
              .rename(columns=dict(zip(all_cols, all_cols2)))
              .assign(**{f"{col}_split": partial(_smiles_splitter, col=col) for col in cols2}))
    table = table.merge(table2, how='cross').query(f"{id_col[0]} < {id_col[0]}_2")

    index_labels = tuple("_".join(items) for items in product(["r1", "p1"], ["r2", "p2"]))
    # print_err(f"By shared chemicals...")
    # table = table.assign(**{f"connection_{index_labels[i]}": partial(_have_shared_chemical, a=f"{a}_split", b=f"{b}_split") 
    #                         for i, (a, b) in enumerate(product(cols, cols2))})
    print_err(f"By maximum similarity...")
    table = table.assign(**{f"max_similarity_{index_labels[i]}": partial(_max_similarity, a=f"{a}_split", b=f"{b}_split") 
                            for i, (a, b) in enumerate(product(cols, cols2))})
    print_err(f"By shared chemicals...")
    table = table.assign(**{f"connection_{index_labels[i]}": partial(_have_shared_chemical_simple, col=f"max_similarity_{index_labels[i]}") 
                            for i, (a, b) in enumerate(product(cols, cols2))})

    return table[[col for col in table if not col.endswith("_split")]]
    
    
def main() -> None:

    tqdm.pandas()

    print_err(f"Reading table from STDIN:")
    table = pd.read_csv(sys.stdin, sep='\t').assign(_id=lambda x: x["rhea_reaction_id"].astype(str).str.cat(x["uniprot_id"], sep=":"))
    print_err(table)
    id_table = table[["_id", "rhea_reaction_id", "uniprot_id"]].copy()

    print_err(f"Cleaning metabolites...")
    table = clean_metabolites(table, cols=("reactants", "products"))
    print_err(f"Connecting metabolites...")
    table = (connect_metabolites(table, cols=("reactants", "products"), id_col=["_id"])
             .merge(id_table, how='left')
             .merge(id_table.rename(columns={col: f"{col}_2" for col in id_table}), how='left'))
    cols_to_rename = ("reactants", "products", "uniprot_id", "rhea_reaction_iduniprot_id")
    table = (table
             .drop(columns=["_id", "_id_2"])
             .rename(columns=dict(zip(cols_to_rename, (f"{col}_1" for col in cols_to_rename)))))
    table[reversed(table.columns)].to_csv(sys.stdout, sep='\t', index=False)

    return None


if __name__ == '__main__':
    main()
