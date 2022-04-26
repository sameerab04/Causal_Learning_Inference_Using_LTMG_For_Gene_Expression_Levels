
import subprocess 
import sys

LTMG_model = sys.argv[1]
LTMG_script = sys.argv[2]
input_file = sys.argv[3]
gene_name = sys.argv[4]
output_tsv_file = sys.argv[5]
output_gene_name_file = sys.argv[6]
ground_truth_python_file = sys.argv[7]
ref_network = sys.argv[8]
causal_learning_file = sys.argv[9]
method = sys.argv[10]
adjacency_file = sys.argv[11]


ltmg_model_command = "Rscript --vanilla " + LTMG_model
res1 = subprocess.call(ltmg_model_command, shell=True)
print("Done running LTMG model")


ltmg_script_command = "Rscript --vanilla " + LTMG_script + " " + input_file + " " + gene_name + " " + output_tsv_file + " " + output_gene_name_file
res2 = subprocess.call(ltmg_script_command, shell=True)
print("Done running LTMG Script")


ground_truth_command = "python3 " + ground_truth_python_file + " " + ref_network + " " + output_gene_name_file
res3 = subprocess.call(ground_truth_command, shell=True)
print("Done creating ground truth graph")


causal_learning_command = "python3 " + causal_learning_file + " " + output_tsv_file + " " + method + " -geneID " + output_gene_name_file + " -eval " + adjacency_file 
res4 = subprocess.call(causal_learning_command, shell=True)
print(res4)




#### Sample Command 
#python3 causal_learning_script.py "LTMG_function_2.R" "LTMG_script.R" "BEELINE-data/inputs/Curated/HSC/HSC-2000-1-50/ExpressionData.csv" "E1736_220" 
#"HSC-2000-1-50.tsv" "HSC-2000-1-50_gene_names.txt" "Groundtruth.py" "BEELINE-data/inputs/Curated/HSC/HSC-2000-1-50/refNetwork.csv" "causal_learning.py" "PC" 
#"adjacency_matrix2.txt""