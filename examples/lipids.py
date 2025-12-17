import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from rdkit import Chem

import CFG

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0


def sanitize_smiles(smiles: str) -> str:
    """
    Parse and sanitize a SMILES (Simplified Molecular Input Line Entry System) string.

    Parameters:
        smiles (str): A SMILES string representing a chemical structure.

    Returns:
        str: A sanitized SMILES string.

    Raises:
        ValueError: If the input SMILES string cannot be parsed into a valid molecule.
    """
    # Parse and sanitize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Failed to parse SMILES: " + smiles)
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def build_hydrocarbon(hc_length: int, dbl_pos: int) -> str:
    """
    Build a hydrocarbon chain as a SMILES string.

    Parameters:
        hc_length (int): The number of carbon atoms in the hydrocarbon chain, excluding the carboxyl carbon.
        dbl_pos (int): The position of the double bond in the chain (1-based index).
                       If 0, no double bond is added.

    Returns:
        str: A SMILES string representing the hydrocarbon chain.
    """
    # Create a list of carbon atoms
    atoms = ["C"] * hc_length
    # Initialize a list of bonds between the atoms
    bonds = [""] * (hc_length - 1)
    # Add a double bond at the specified position, if applicable
    if dbl_pos != 0:
        idx = dbl_pos - 2  # Convert 1-based index to 0-based index
        bonds[idx] = "="
    # Stitch the atoms and bonds together into a SMILES string
    s = atoms[0]
    for b, a in zip(bonds, atoms[1:]):
        s += b + a
    return s


def generate_glyceride_smiles(chain_length: int, double_bond_pos: int) -> str:
    """
    Generate a SMILES string for a glyceride molecule.

    Parameters:
        chain_length (int): The number of carbon atoms in the hydrocarbon chain, excluding the carboxyl carbon.
        double_bond_pos (int): The position of the double bond in the chain (1-based index).
                               If 0, no double bond is added.

    Returns:
        str: A sanitized SMILES string representing the glyceride molecule.
    """
    # Build hydrocarbon fragment
    tail = build_hydrocarbon(chain_length, double_bond_pos)

    # The carboxylic acid group:
    smiles = f"OC(=O){tail}"

    # Parse and sanitize
    return sanitize_smiles(smiles)


def generate_triglyceride_smiles(chain_length: int, double_bond_pos: int) -> str:
    """
    Generate a SMILES string for a triglyceride molecule.

    Parameters:
        chain_length (int): The number of carbon atoms in the hydrocarbon chain, excluding the carboxyl carbon.
        double_bond_pos (int): The position of the double bond in the chain (1-based index).
                               If 0, no double bond is added.

    Returns:
        str: A sanitized SMILES string representing the triglyceride molecule.
    """
    # Build hydrocarbon fragment
    chain_s = build_hydrocarbon(chain_length, double_bond_pos)

    # The glycerol backbone with three ester linkages:
    smiles = (f"C(C(OC(=O){chain_s})COC(=O){chain_s})OC(=O){chain_s}")

    # Parse and sanitize
    return sanitize_smiles(smiles)


def generate_phospholipid_smiles(chain_length: int, chain_double_bond_pos: int) -> str:
    """
    Generate a SMILES string for a phospholipid molecule (Phosphatidylcholine).

    Parameters:
        chain_length (int): The number of carbon atoms in the hydrocarbon chain, excluding the carboxyl carbon.
        chain_double_bond_pos (int): The position of the double bond in the chain (1-based index).
                                     If 0, no double bond is added.

    Returns:
        str: A sanitized SMILES string representing the phospholipid molecule.
    """
    # Build hydrocarbon fragment
    frag = build_hydrocarbon(chain_length, chain_double_bond_pos)

    # Assemble the SMILES for the Phosphatidylcholine
    head_group = "OP(=O)([O-])OCC[N+](C)(C)C"
    smiles = f"{frag}C(=O)OC(COC(=O){frag})C{head_group}"

    # Parse and sanitize
    return sanitize_smiles(smiles)


def get_bounds(smiles_list):
    """
    Calculate the assembly index bounds for a list of SMILES strings.

    Parameters:
        smiles_list (list of str): List of SMILES strings representing molecules.

    Returns:
        tuple: A tuple containing:
            - list of Assembly index upper bounds calculated using the CFG DET method.
            - list of Lower bounds on the assembly indices from addition chains.
    """
    # Convert all the SMILES to NetworkX graphs without hydrogen
    graphs = [CFG.smi_to_nx(smile, add_hydrogens=False) for smile in smiles_list]

    # Initialize lists to store the results
    L_det_list = []
    L_ac_list = []

    # Loop over the graphs, calculate AI bounds
    for i, g in enumerate(graphs):
        print(flush=True)
        L_det, _, _ = CFG.calculate_assembly_path_det(g)
        L_det_list.append(L_det)
        L_ac_list.append(l_n[len(g.edges)])

    return L_det_list, L_ac_list


if __name__ == "__main__":
    l_n = [0]
    base_path = os.path.dirname(__file__)
    with open(os.path.join(base_path, 'ln_9999.txt'), 'r') as file: # This file contains precomputed shortest addition chain values
        # Skip the header line
        next(file)
        next(file)

        for line in file:  # Read each line in the file
            # Split the line into columns based on whitespace
            columns = line.split()
            l_n.append(int(columns[3]))

    carbon_lengths = range(2, 18)
    fatty_acid_smiles = [generate_glyceride_smiles(c, 0) for c in carbon_lengths]
    triglyceride_smiles = [generate_triglyceride_smiles(c, 0) for c in carbon_lengths]
    phospholipid_smiles = [generate_phospholipid_smiles(c, 0) for c in carbon_lengths]

    ai_list = [3, 4, 4, 5, 5, 6, 5, 6, 6, 7, 6, 7, 7, 7, 6, 7]
    L_det_list, ai_lower_list = get_bounds(fatty_acid_smiles)
    
    ai_list_trig = [7, 8, 8, 9, 9, 10, 9, 10, 10, 11, 10, 11, 11, 11, 10, 11]
    L_det_list_trig, ai_lower_list_trig = get_bounds(triglyceride_smiles)
    
    ai_list_phos = [14, 15, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 18, 18]
    L_det_list_phos, ai_lower_list_phos = get_bounds(phospholipid_smiles)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))

    xlab = 'Lipid length'
    ylab = 'Joining Operations'
    xs = 16

    # Plot the graphs
    ax1.plot(carbon_lengths, L_det_list, 'o:', label='RePair Upper Bound', lw=2, color='black')
    ax1.plot(carbon_lengths, ai_list, 'o-', label='Assembly Index', lw=2, color='black')
    ax1.plot(carbon_lengths, ai_lower_list, 'o:', label='Integer Chain Lower Bound', lw=2, color='black')
    ax1.set_title("Fatty Acid", fontsize=xs)
    ax1.legend(loc="lower right")
    ax1.set_xticks(carbon_lengths[::2])  # Set ticks to every second value
    ax1.set_xticklabels(carbon_lengths[::2])  # Show only every second label
    ax1.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.set_xlabel(xlab, fontsize=xs)
    ax1.set_ylim(0, max(L_det_list) + 1)
    ax1.set_ylabel(ylab, fontsize=xs)

    ax2.plot(carbon_lengths, L_det_list_trig, 'o:', label='RePair Upper Bound', lw=2, color='blue')
    ax2.plot(carbon_lengths, ai_list_trig, 'o-', label='Assembly Index', lw=2, color='blue')
    ax2.plot(carbon_lengths, ai_lower_list_trig, 'o:', label='Integer Chain Lower Bound', lw=2, color='blue')
    ax2.set_title("Triglyceride", fontsize=xs)
    ax2.legend(loc="lower right")
    ax2.set_xticks(carbon_lengths[::2])  # Set ticks to every second value
    ax2.set_xticklabels(carbon_lengths[::2])  # Show only every second label
    ax2.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
    ax2.tick_params(axis='both', which='both', top=True, right=True)
    ax2.set_xlabel(xlab, fontsize=xs)
    ax2.set_ylim(0, max(L_det_list_trig) + 1)
    ax2.set_ylabel(ylab, fontsize=xs)

    ax3.plot(carbon_lengths, L_det_list_phos, 'o:', label='RePair Upper Bound', lw=2, color='red')
    ax3.plot(carbon_lengths, ai_list_phos, 'o-', label='Assembly Index', lw=2, color='red')
    ax3.plot(carbon_lengths, ai_lower_list_phos, 'o:', label='Integer Chain Lower Bound', lw=2, color='red')
    ax3.set_title("Phospholipid", fontsize=xs)
    ax3.legend(loc="lower right")
    ax3.set_xticks(carbon_lengths[::2])  # Set ticks to every second value
    ax3.set_xticklabels(carbon_lengths[::2])  # Show only every second label
    ax3.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
    ax3.tick_params(axis='both', which='both', top=True, right=True)
    ax3.set_xlabel(xlab, fontsize=xs)
    ax3.set_ylim(0, max(L_det_list_phos) + 1)
    ax3.set_ylabel(ylab, fontsize=xs)

    # Restrict y-axis to integer ticks
    for ax in [ax1, ax2, ax3]:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig('lipids.png', dpi=600)
    plt.savefig('lipids.pdf')
    plt.show()
