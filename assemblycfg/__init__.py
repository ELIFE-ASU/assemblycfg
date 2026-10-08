"""Upper and lower bounds on string and molecular assembly index."""
from ._version import __version__

# Upper bounds on string assembly index
from .string_repair import ai_core, repair_with_pathways

# Upper bounds on molecular assembly index
from .molecule_repair import GraphRepairResult, calculate_assembly_path_graph_repair, graph_repair

# Lower bounds on string and molecular assembly index
from .lz import lz_lower_bound
from .vac import (
    VacError,
    find_vac,
    install_vac,
    mol_unit_counts,
    scalar_chain_length,
    solve_vac,
    string_unit_counts,
    vac_lower_bound,
)

# Molecule conversion
from .molecules import (
    bond_order_rdkit_to_int,
    dict_to_nx,
    get_disconnected_subgraphs,
    mol2graph,
    mol_to_nx,
    molfile_to_mol,
    nx_to_dict,
    print_graph_dict,
    print_virtual_objects,
    remove_hydrogen_from_graph,
    safe_standardize_mol,
    smi_to_mol,
    smi_to_nx,
    standardize_mol,
)
