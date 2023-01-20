# Executing SigProfilerExtractor to get SBS signatures of all patients.
# Reference : https://github.com/AlexandrovLab/SigProfilerExtractor
import sys
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall

def install_ref(ref_path):
  genInstall.install('GRCh37', offline_files_path=ref_path)

def SBS_NSLC():
    # Main function with default values, except cpu usage and maximum_signatures.
    sig.sigProfilerExtractor("vcf", "NSLC_SBS", "Data/92_patients_NSLC_filtered_VCFS", "GRCh37", cpu=8, minimum_signatures=1, maximum_signatures=8)

if __name__=="__main__":
    #SBS_NSLC()
    ref_path = sys.argv[1]
    install_ref(ref_path)


