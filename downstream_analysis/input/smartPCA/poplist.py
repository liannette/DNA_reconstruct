
# pop_set = set()

# with open("/home/projects/DNA_reconstruct/downstream_analysis/input/genotypes/data.ind", "r") as infile:
#     for line in infile:
#         population = line.split()[-1]
#         pop_set.add(population)

# with open("poplist.txt", "w") as outfile:
#     for pop in sorted(pop_set):
#         print(pop, file=outfile)

pop_dict = dict()

with open("/home/projects/DNA_reconstruct/downstream_analysis/input/smartPCA/poplist_allmodern.txt", "r") as infile:
    for line in infile:
        population = line.split()[-1]
        pop_dict[population] = 0
        
with open("/home/projects/DNA_reconstruct/downstream_analysis/input/genotypes/data.ind", "r") as infile:
    for line in infile:
        population = line.split()[-1]
        try:
            pop_dict[population] +=1
        except:
            pass

with open("poplist_allmodern_count.txt", "w") as outfile:
    for population in sorted(pop_dict.keys()):
        print(str(pop_dict[population]) + "\t" +  population, file=outfile)
