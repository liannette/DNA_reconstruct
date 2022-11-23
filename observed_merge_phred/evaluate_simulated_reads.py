import os 
import pandas as pd


def count_phred_occurences(fastq_file):
    
    phred_counter = 42 * [0]
    
    cnt = 0
    with open(fastq_file, 'rb') as infile:
        while True:
            # read line
            line = infile.readline().rstrip()
            cnt += 1
            # get quality score lines
            if cnt == 4:
                for char in line:
                    qual = char - 33
                    phred_counter[qual] += 1
                cnt = 0 
                
            # Stop while loop at end of file
            if len(line) == 0:
                break

    return phred_counter
    

def main():
    
    results = {
        "quality_shift": list(),
        "read": list(),
        "quality_score": list(),
        "count": list(),
        }
    
    for q_shift in [0, -10, -20, -30]:
        for read in ["s1", "s2"]:
            infile = os.path.join(
                "output", 
                "simulations", 
                f"gen_n_100000000_l_120_qs_{q_shift}_{read}.fq",
                )
        
            phred_counter = count_phred_occurences(infile)
            for q_score in range(len(phred_counter)):
                results["quality_shift"].append(q_shift)
                results["read"].append(read)
                results["quality_score"].append(q_score)
                results["count"].append(phred_counter[q_score])
                
    
    
    export_path = os.path.join("output", "evaluation", "simulated_reads.csv")
    
    if export_path is not None:
        df = pd.DataFrame.from_dict(results)
        df.to_csv(export_path, na_rep="NA", index=False)
    

if __name__ == "__main__":
    main()