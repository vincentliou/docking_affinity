#%%%
# version 3: import the fixed section 2 and 3 here
# version 4: using env 
from dotenv import load_dotenv
import os
import subprocess
import requests
import pandas as pd
import glob
import shutil
from prodigy_prot.modules.prodigy import Prodigy



def run_zdock(receptor_pdb, ligand_pdb,
              zdock_path="./zdock",
              out_file="zdock_out3.txt"):
    """
    Run ZDOCK docking between receptor and ligand.
    """
    cmd = [zdock_path, "-R", receptor_pdb, "-L", ligand_pdb, "-o", out_file]
    subprocess.run(cmd, check=False)
    return out_file

def extract_poses(
                  receptor_pdb,
                  ligand_pdb,
                  create_pl_path,
                  output_dir,
                  zdock_out,
                  num_poses=10):
    
    """
    Extract top poses into PDB format using create.pl (Perl wrapper for create_lig).
    """
    script_dir = os.path.dirname(create_pl_path)
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Change the current working directory to where create.pl and create_lig are located
    original_cwd = os.getcwd()
    os.chdir(script_dir)

    # Use a try...finally block to ensure we always return to the original directory
    try:
        # Construct the command. The Perl script will create files in the current directory.
        cmd = ["perl", create_pl_path, os.path.abspath(zdock_out), os.path.abspath(receptor_pdb), os.path.abspath(ligand_pdb), str(num_poses)]
        
        print(f"Executing command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # After the script runs, the poses are created in script_dir.
        # Now, we need to move them to the desired output_dir.
        poses = []
        for i in range(1, num_poses + 1):
            source_file = f"complex.{i}.pdb"
            dest_file = os.path.join(output_dir, f"pose.{i}.pdb")
            
            if os.path.exists(source_file):
                shutil.move(source_file, dest_file)
                poses.append(dest_file)
            else:
                print(f"Warning: Expected file {source_file} not found. Something may have gone wrong with the Perl script.")
                
    except subprocess.CalledProcessError as e:
        print(f"Error executing Perl script: {e}")
        poses = []
    finally:
        # Always change back to the original working directory
        os.chdir(original_cwd)

    return poses


def prodigy_affinity(pdb_file,selection="J M"):
    """
    Call PRODIGY API to estimate binding affinity (ΔG and Kd).
    """
    """
    NEW: Call the local PRODIGY library to estimate binding affinity (ΔG and Kd).
    This function no longer uses a web API.
    """
    try:
        # The predict_complex_dG_dKd function takes a single PDB file of the complex
        # and returns a tuple: (ΔG, Kd)


        #print(Prodigy(pdb_file))
        cmd = ["prodigy", pdb_file , "-q"]
        # Add the selection arguments if provided
        if selection:
            cmd.append("--selection")
            cmd.extend(selection.split())
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # The output looks like: pose.1_model0 -15.228
        # We need to extract the numerical value
        output_line = res.stdout.strip()
        affinity_dg = float(output_line.split()[-1])
        
        # Calculate Kd from ΔG using the formula: ΔG = RT ln(Kd)
        # R = 1.987 cal/mol·K, T = 298.15 K (25°C)
        # We need to use R in kcal/mol·K, so R = 0.001987
        R = 0.001987
        T = 298.15
        
        kd = 2.71828**(affinity_dg / (R * T))
        print(kd)
        return {"dg": affinity_dg, "kd": kd}
    except Exception as e:
        # Handle potential errors during the prediction
        return {"error": str(e)}   
    return res



def docking_pipeline(receptor_pdb, ligand_pdb,
                     create_pl_path, output_dir,zdock_path="./zdock",num_poses=5):
    """
    Full docking pipeline:
    1. Run ZDOCK
    2. Extract top poses
    3. Send each pose to PRODIGY for affinity estimation
    """
    print("[1] Running ZDOCK...")
    zdock_out = run_zdock(receptor_pdb, ligand_pdb, zdock_path=zdock_path)

    print("[2] Extracting poses...")
    poses = extract_poses(receptor_pdb, ligand_pdb,create_pl_path, output_dir, zdock_out, num_poses=num_poses)

    print("[3] Scoring poses with PRODIGY...")
    results = []
    for i, pose in enumerate(poses, start=1):
        print(f"   → Pose {i}")
        res = prodigy_affinity(pose)
        if "dg" in res:
            results.append({
                "pose": i,
                "pose_file": pose,
                "ΔG (kcal/mol)": res["dg"],
                "Kd (M)": res["kd"]
            })
        else:
            results.append({
                "pose": i,
                "pose_file": pose,
                "error": res.get("error", "Unknown error")
            })

    df = pd.DataFrame(results)
    df.to_csv("docking_results.csv", index=False)
    print("\nDocking complete. Results saved to docking_results.csv")
    return df


if __name__ == "__main__":
    receptor = "5SX4_EGFR_no.pdb"
    ligand = "5SX4_antibody_no.pdb"
    load_dotenv()  # Load environment variables from .env file
    create_pl_path = os.getenv("create_pl_path")
    output_dir = os.getenv("output_dir")
    results = docking_pipeline(receptor, ligand,create_pl_path, output_dir,zdock_path="./zdock", num_poses=5)
    print(results)
#%%%
