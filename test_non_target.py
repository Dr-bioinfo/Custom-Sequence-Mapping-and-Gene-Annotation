import os
import pytest
import subprocess
#from ID_timediff_Info import calculate_time_difference, get_sec

# Tst cases
test_cases = [
    (
        
        "-f1 /home/alphabox0006/Jerry/Unit_testing/off_target/input/resfinder.tsv -f2 /home/alphabox0006/Jerry/Unit_testing/off_target/input/S01_emu_best_hits.csv -f3 /home/alphabox0006/Jerry/Unit_testing/off_target/input/S01_ITS_best_hits.csv -sam /home/alphabox0006/Jerry/Unit_testing/off_target/input/test.fasta -seq /home/alphabox0006/Jerry/Unit_testing/off_target/input/sequence.fasta -gb /home/alphabox0006/Jerry/Unit_testing/off_target/input/sequence.gb -targ /home/alphabox0006/Jerry/Unit_testing/off_target/output/targeted_seq.txt -non /home/alphabox0006/Jerry/Unit_testing/off_target/output/nontargeted_seq.fasta -fin /home/alphabox0006/Jerry/Unit_testing/off_target/output/Batch_Ecoli_OffTargeted.csv",
        "/home/alphabox0006/Jerry/Unit_testing/off_target/output/targeted_seq.txt, /home/alphabox0006/Jerry/Unit_testing/off_target/output/nontargeted_seq.fasta , /home/alphabox0006/Jerry/Unit_testing/off_target/output/Batch_Ecoli_OffTargeted.csv"
    ),
]
# Create a function to run the script with the specified arguments and capture the output
def run_script(args):
    result = subprocess.run(["python", "non_target_auto.py"] + args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Check the return code and raise an exception if non-zero
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, " ".join(result.args), result.stderr)

    stdout, stderr = result.stdout, result.stderr
    print("standard output: ", stdout)
    print("standard error: ", stderr)

# Create a test function using pytest.mark.parametrize
@pytest.mark.parametrize("input_args, expected_output_file", test_cases)
def test_non_target(input_args, expected_output_file):
    # Run the script with the specified arguments
    stdout, stderr = run_script(input_args)
    print("standard output: " ,stdout)
    #print("standard error : ",stderr)
    assert stdout is not None
    assert stderr == ""
    # Check if the expected output file exists
    assert os.path.exists(expected_output_file)
