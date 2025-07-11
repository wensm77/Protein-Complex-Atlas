from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
import gemmi
import os
import concurrent.futures
import argparse
import subprocess
from functools import partial

def get_chain_sequence(chain):
    return ''.join(seq1(residue.resname) for residue in chain if is_aa(residue, standard=True))

def get_poly_seq(structure):
    poly_seq = []
    entity_id = 1
    processed_sequences = {}
    entity_id_map = {}
    one_letter_sequences = {}
    
    for model in structure:
        for chain in model:
            chain_seq = get_chain_sequence(chain)
            if chain_seq in processed_sequences:
                entity_id_map[chain.id] = processed_sequences[chain_seq]
                continue
            
            entity_id_map[chain.id] = entity_id
            one_letter_sequences[chain.id] = chain_seq
            
            for i, residue in enumerate(chain, start=1):
                if is_aa(residue, standard=True):
                    poly_seq.append({
                        "entity_id": entity_id,
                        "num": i,
                        "mon_id": residue.resname,
                        "hetero": "n"
                    })
            processed_sequences.update({chain_seq:entity_id})
            entity_id += 1
    return poly_seq, entity_id_map, one_letter_sequences

def add_entity_poly(doc, poly_seq, one_letter_sequences, entity_id_map):
    block = doc.sole_block()
    entity_poly = block.init_loop('_entity_poly.', ['entity_id', 'pdbx_strand_id', 'type', 'pdbx_seq_one_letter_code_can'])
    
    seen = set()

    for seq in poly_seq:
        entity_id = seq["entity_id"]
        chain_id = [k for k, v in entity_id_map.items() if v == entity_id][0]
        one_letter_seq = one_letter_sequences[chain_id]

        entry = (str(seq["entity_id"]), chain_id, 'polypeptide(L)', one_letter_seq)
        
        if entry not in seen:
            entity_poly.add_row(entry)
            seen.add(entry)

def add_pdbx_poly_seq_scheme(doc, bio_structure, entity_id_map):
    block = doc.sole_block()
    pdbx_poly_seq_scheme = block.init_loop('_pdbx_poly_seq_scheme.', [
        'asym_id', 'seq_id', 'mon_id', 'pdb_seq_num', 'auth_seq_num', 'pdb_strand_id'
    ])
    
    for model in bio_structure:
        for chain in model:
            entity_id = entity_id_map[chain.id]
            for i, residue in enumerate(chain, start=1):
                if is_aa(residue, standard=True):
                    pdbx_poly_seq_scheme.add_row([chain.id, str(i), residue.resname, str(i), str(i), chain.id])

def add_atom_site(cif_doc, structure):
    block = cif_doc.sole_block()
    cif_atom_site = block.init_loop('_atom_site.', [
        'group_PDB', 'label_atom_id', 'type_symbol', 'label_comp_id',
        'label_asym_id', 
        'label_seq_id', 'label_alt_id', 'Cartn_x', 'Cartn_y', 'Cartn_z',
        'occupancy', 'B_iso_or_equiv', 'pdbx_PDB_model_num'
    ])

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_name = atom.name
                    type_symbol = atom.element.name
                    resname = residue.name
                    chain_id = chain.name

                    seq_id = residue.seqid
                    x, y, z = atom.pos
                    occupancy = atom.occupancy() if hasattr(atom, 'occupancy') else 1.0
                    B_factor = atom.bfactor() if hasattr(atom, 'bfactor') else 0.0
                    model_number = model.name

                    cif_atom_site.add_row([
                        'ATOM', atom_name, type_symbol, resname,
                        chain_id,
                        str(seq_id), '.',
                        str(x), str(y), str(z),
                        str(occupancy), str(B_factor), str(model_number)
                    ])

def add_pdbx_struct_mod_residue(doc, structure):
    block = doc.sole_block()
    mod_residue = block.init_loop('_pdbx_struct_mod_residue.', ['label_asym_id', 'label_seq_id', 'label_comp_id', 'details'])

    for model in structure:
        for chain in model:
            chain_name = chain.name
            for residue in chain:
                label_asym_id = chain_name
                label_seq_id = residue.seqid
                label_comp_id = residue.name
                
                details = 'none'

                mod_residue.add_row([label_asym_id, str(label_seq_id), label_comp_id, details])

def add_other_sections(doc):
    block = doc.sole_block()
    
    exptl = block.init_loop('_exptl.', ['method'])
    exptl.add_row(['X-RAY DIFFRACTION'])

    pdbx_database_status = block.init_loop('_pdbx_database_status.', ['status', 'recvd_initial_deposition_date'])
    pdbx_database_status.add_row(['REL', '2000-10-06'])

    em_3d_reconstruction = block.init_loop('_em_3d_reconstruction.', ['details', 'resolution'])
    
    details = "EM reconstruction details go here."
    resolution = 2.0
    em_3d_reconstruction.add_row([details, str(resolution)])

def add_struct_assembly_gen(doc, entity_id_map):
    block = doc.sole_block()
    struct_assembly_gen = block.init_loop('_pdbx_struct_assembly_gen.', ['assembly_id', 'asym_id_list'])

    assembly_id = '1'
    asym_id_list = ','.join(entity_id_map.keys())

    struct_assembly_gen.add_row([assembly_id, asym_id_list])

def add_struct_assembly(doc):
    block = doc.sole_block()
    struct_assembly = block.init_loop('_pdbx_struct_assembly.', ['id', 'details', 'method_details'])
    
    assembly_id = '1'
    details = 'detail_default.'
    method_details = 'X-RAY DIFFRACTION'
    struct_assembly.add_row([assembly_id, details, method_details])

def convert_pdb_to_cif(input_dir, pdb_file, out_path, pt_path):
    parser = PDBParser()
    bio_structure = parser.get_structure('structure', input_dir + pdb_file)
    poly_seq, entity_id_map, one_letter_sequences = get_poly_seq(bio_structure)

    structure = gemmi.read_structure(input_dir + pdb_file)
    structure.setup_entities()
    structure.assign_label_seq_id()
    cif_doc = structure.make_mmcif_document()

    add_entity_poly(cif_doc, poly_seq, one_letter_sequences, entity_id_map)
    add_pdbx_poly_seq_scheme(cif_doc, bio_structure, entity_id_map)
    add_atom_site(cif_doc, structure)
    add_pdbx_struct_mod_residue(cif_doc, structure)
    add_other_sections(cif_doc)
    add_struct_assembly_gen(cif_doc, entity_id_map)
    add_struct_assembly(cif_doc)

    out_put_name = out_path + pdb_file.split('.')[0] + '.cif'
    cif_doc.write_file(out_put_name)
    if os.path.exists(out_put_name):
        print(f"{out_put_name} saved!")
    else:
        print(f"{out_put_name} failed!")
    
    gz_out_name = out_put_name + '.gz'
    subprocess.run(['gzip', '-k', out_put_name], check=True)
    print(f'Compressed to {gz_out_name}')
    
    py_script = 'parse_cif_noX_qxz.py'
    subprocess.run(['python', py_script, gz_out_name, pt_path + pdb_file.split('.')[0]], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDB files to CIF format with multi-threading support")
    parser.add_argument("--input_dir", required=True, help="Directory containing PDB files")
    parser.add_argument("--output_dir", required=True, help="Directory to save CIF files")
    parser.add_argument("--pt_path", required=True, help="Directory to save PT files")
    parser.add_argument("--threads", type=int, default=30, help="Number of threads to use")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    if not os.path.exists(args.pt_path):
        os.makedirs(args.pt_path)
    
    print(f"Start processing PDB files in {args.input_dir} with {args.threads} threads")
    
    pdb_files = [f for f in os.listdir(args.input_dir) if f.endswith('.pdb')]
    count = 0
    for pdb_file in pdb_files:
        count = count + 1
        print(f"{count}: Processing {pdb_file}")
        convert_pdb_to_cif(args.input_dir, pdb_file, args.output_dir, args.pt_path)