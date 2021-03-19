import os
import sys
import argparse
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
import matplotlib.pyplot as plt
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
import numpy as np
import MDAnalysis
from MDAnalysis.analysis import rms, align
from PIL import Image

def file_verification(file_name):
    # verifie l'existence des fichiers
    if os.path.exists(file_name) != 1:
        sys.exit("Error: {} does not exist.".format(file_name))

def get_arguments():
    # recupere les arguments de la ligne de commande
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-pu', dest='PU_name',
                        required=True, help="PU name")
    parser.add_argument('-nbr', dest='nb_rep',
                        required=True, help="Number of replicas")
    return parser.parse_args()
                        
def verify_index(PU_name, replica):
    # verifie qu'un fichier d'index existe dans le dossier de la PU, cree ce fichier le cas echeant
    PU_directory_input = "data/PUs_sele_pdb_Replicat_" + replica + "/" + PU_name + "/"
    if os.path.exists(PU_directory_input + PU_name + '.ndx') != 1:
        bashCommand = 'echo "q" | gmx make_ndx -f ' + PU_directory_input + PU_name + '_NPT_PRO.gro -o ' + PU_directory_input + PU_name + '.ndx'
        print(bashCommand)
        os.system(bashCommand)
        
def generate_index_file(PU_name, replica):
    # genere un fichier index contenant uniquement la PU
    PU_directory_input = "data/PUs_sele_pdb_Replicat_" + replica + "/" + PU_name + "/"
    PU_directory_output = "outputs/" + PU_name + "/Replica_" + replica + "/"
    with open(PU_directory_input + PU_name + '.ndx', mode='r') as in_file, \
         open(PU_directory_output + 'index.ndx', mode='w') as out_file:
        flag_to_write = 0
        for line in in_file:
            if line.startswith('['):
                flag_to_write = 0
            if line.startswith('[ Protein ]'):
                flag_to_write = 1
            if flag_to_write == 1:
                out_file.write(line)

def generate_PU_files_without_water(PU_name, replica):
    # genere des fichiers .gro et .xtc contenant uniquement la PU
    PU_directory_input = "data/PUs_sele_pdb_Replicat_" + replica + "/" + PU_name + "/"
    PU_directory_output = "outputs/" + PU_name + "/Replica_" + replica + "/"
    if os.path.exists(PU_directory_output + PU_name + '_analysis.gro') != 1:
        bashCommand = 'gmx trjconv -f ' + PU_directory_input + PU_name + '_NPT_PRO.gro -s ' + PU_directory_input + PU_name + '_NPT_PRO.tpr -o ' + PU_directory_output + PU_name + '_analysis.gro -n ' + PU_directory_output + 'index.ndx -boxcenter rect -pbc mol -center'
        print(bashCommand)
        os.system(bashCommand)
    if os.path.exists(PU_directory_output + PU_name + '_analysis.xtc') != 1:
        bashCommand = 'gmx trjconv -f ' + PU_directory_input + PU_name + '_NPT_PRO.xtc -s ' + PU_directory_input + PU_name + '_NPT_PRO.tpr -o ' + PU_directory_output + PU_name + '_analysis.xtc -n ' + PU_directory_output + 'index.ndx -boxcenter rect -pbc mol -center -skip 10'
        print(bashCommand)
        os.system(bashCommand)
        
def generate_pdb_visualization(PU_name, replica):
    # genere un image pymol de la PU
    print("Processing visualization...")
    PU_directory_output = "outputs/" + PU_name + "/Replica_" + replica + "/"
    pymol.cmd.load(PU_directory_output + PU_name + '_analysis.gro', PU_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable(PU_name)
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)
    pymol.cmd.bg_color('white')
    pymol.cmd.color('green', 'all')
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png("%s.png"%(PU_directory_output + PU_name))
    print("Done!")
        

#### Radius of gyration ####
def radius_of_gyration(u, PU_directory_output):
    # calcule le rayon de giration de la PU
    print("Processing radius of gyration...")
    Rgyr = []
    protein = u.select_atoms("protein")
    for ts in u.trajectory:
       Rgyr.append((u.trajectory.time, protein.radius_of_gyration()))
    Rgyr = np.array(Rgyr)

    plt.figure()
    ax1 = plt.subplot(111)
    ax1.plot(Rgyr[:,0], Rgyr[:,1], 'r-', label=r"$R_G$")
    ax1.legend(loc="best")
    ax1.set_title("Radius of gyration $R_G$")
    ax1.set_xlabel("Time (ps)")
    ax1.set_ylabel(r"Radius of gyration $R_G$ ($\AA$)")
    ax1.figure.savefig(PU_directory_output + "RGYR.png")
    print("Done!")

#### RMSD ####
def RMSD(u, PU_directory_output):
    # calcule le RMSD de la PU
    print("Processing RMSD...")
    Rmsd = rms.RMSD(u)
    Rmsd.run()
    rmsd = Rmsd.rmsd.T   # utilisation de la transposée pour l'affichage
    time = rmsd[1]
    
    plt.figure()
    ax2 = plt.subplot(111)
    ax2.plot(time, rmsd[2], 'g-',  label="RMSD")
    ax2.legend(loc="best")
    ax2.set_title("RMSD")
    ax2.set_xlabel("Time (ps)")
    ax2.set_ylabel(r"RMSD ($\AA$)")
    ax2.figure.savefig(PU_directory_output + "RMSD.png")
    print("Done!")
    
    
#### RMSF ####
def RMSF(u, PU_directory_output):
    # calcule le RMSF de la PU
    print("Processing RMSF...")
    average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run()
    ref = average.universe
    aligner = align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()
    c_alphas = u.select_atoms('protein and name CA')
    R = rms.RMSF(c_alphas).run()
    
    plt.figure()
    ax3 = plt.subplot(111)
    ax3.plot(c_alphas.resids, R.rmsf, 'b-',  label="RMSF")
    ax3.legend(loc="best")
    ax3.set_title("RMSF")
    ax3.set_xlabel("Residue")
    ax3.set_ylabel(r"RMSF ($\AA$)")
    ax3.figure.savefig(PU_directory_output + "RMSF.png")
    print("Done!")
    
    
#### Contact map ####
def contact_map(PU_name, PU_directory_output):
    # genere la carte de contact de la PU
    print("Processing contact map...")
    traj = md.load(PU_directory_output + PU_name + "_analysis.xtc", top=PU_directory_output + PU_name + "_analysis.gro")
    topology = traj.topology

    frame_contacts = ContactFrequency(traj)
    
    plt.figure()
    fig, ax4 = frame_contacts.residue_contacts.plot()
    plt.title("Contact map")
    _ = plt.xlabel("Residue")
    _ = plt.ylabel("Residue")
    plt.savefig(PU_directory_output + "Contact_map.png")
    print("Done!")
    
    
#### PBxplore ####
def PBxplore(PU_name, PU_directory_output):
    # analyse par Protein Block de la PU
    print("Processing PBxplore...")
    bashCommand = 'PBassign -x ' + PU_directory_output + PU_name + '_analysis.xtc -g ' + PU_directory_output + PU_name + '_analysis.gro -o ' + PU_directory_output + PU_name
    print(bashCommand)
    os.system(bashCommand)
    
    bashCommand = 'PBcount -f ' + PU_directory_output + PU_name + '.PB.fasta -o ' + PU_directory_output + PU_name
    print(bashCommand)
    os.system(bashCommand)
    
    bashCommand = 'PBstat -f ' + PU_directory_output + PU_name + '.PB.count --map --neq --logo -o ' + PU_directory_output + PU_name
    print(bashCommand)
    os.system(bashCommand)
    print("Done!")

def analysis(PU_name, replica):
    PU_directory_output = "outputs/" + PU_name + "/Replica_" + replica + "/"
    
    u = MDAnalysis.Universe(PU_directory_output + PU_name + "_analysis.gro", PU_directory_output + PU_name + "_analysis.xtc")
    
    # analyses
    radius_of_gyration(u, PU_directory_output)
    RMSD(u, PU_directory_output)
    RMSF(u, PU_directory_output)
    contact_map(PU_name, PU_directory_output)
    PBxplore(PU_name, PU_directory_output)

def generate_summary(PU_name, replica):
    # genere une grande image contenant toutes les images obtenues dans la partie d'analyse
    PU_directory_output = "outputs/" + PU_name + "/Replica_" + replica + "/"
    fig_height = 700
    im1 = Image.open(PU_directory_output + PU_name + '.png')
    im1 = im1.resize((int(fig_height * im1.width / im1.height), fig_height))
    im2 = Image.open(PU_directory_output + 'RGYR.png')
    im2 = im2.resize((int(fig_height * im2.width / im2.height), fig_height))
    im3 = Image.open(PU_directory_output + 'RMSD.png')
    im3 = im3.resize((int(fig_height * im3.width / im3.height), fig_height))
    im4 = Image.open(PU_directory_output + 'RMSF.png')
    im4 = im4.resize((int(fig_height * im4.width / im4.height), fig_height))
    im5 = Image.open(PU_directory_output + 'Contact_map.png')
    im5 = im5.resize((int(fig_height * im5.width / im5.height), fig_height))
    im6 = Image.open(PU_directory_output + PU_name + '.PB.map.png')
    im6 = im6.resize((int(fig_height * im6.width / im6.height), fig_height))
    im7 = Image.open(PU_directory_output + PU_name + '.PB.Neq.png')
    im7 = im7.resize((int(fig_height * im7.width / im7.height), fig_height))
    im8 = Image.open(PU_directory_output + PU_name + '.PB.logo.png')
    im8 = im8.resize((int(fig_height * im8.width / im8.height), fig_height))
    
    dst = Image.new('RGB', (im1.width + im2.width + im3.width + im4.width, im1.height + im5.height + im8.height), 'white')
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.paste(im3, (im1.width + im2.width, 0))
    dst.paste(im4, (im1.width + im2.width + im3.width, 0))
    dst.paste(im5, (0, im1.height))
    dst.paste(im6, (im5.width, im1.height))
    dst.paste(im7, (im5.width + im6.width, im1.height))
    dst.paste(im8, (0, im1.height + im5.height))
    dst.save(PU_directory_output + PU_name + '_summary.png')

def main():
    # recuperation des arguments
    args = get_arguments()
    PU_name = args.PU_name
    nb_rep = int(args.nb_rep)
    
    # verification des la presence des fichiers et creation de l'architecture
    for replica in range(1, nb_rep+1):
    	file_verification('data/PUs_sele_pdb_Replicat_' + str(replica) + '/' + PU_name)
    if os.path.exists('outputs/') != 1:
    	os.mkdir('outputs/')
    if os.path.exists('outputs/' + PU_name) != 1:
        os.mkdir('outputs/' + PU_name)
    for replica in range(1, nb_rep+1):
        if os.path.exists("outputs/" + PU_name + "/Replica_" + str(replica)) != 1:
            os.mkdir("outputs/" + PU_name + "/Replica_" + str(replica))
            
    # analyse de la PU
    for replica in range(1, nb_rep+1):
        print("PROCESSING REPLICA #" + str(replica))
        verify_index(PU_name, str(replica))
        generate_index_file(PU_name, str(replica))
        generate_PU_files_without_water(PU_name, str(replica))
        generate_pdb_visualization(PU_name, str(replica))
        analysis(PU_name, str(replica))
        generate_summary(PU_name, str(replica))

if __name__ == '__main__':
    main()

