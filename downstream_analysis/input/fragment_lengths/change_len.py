# with open("Vi33.19", "r") as infile:
#     lines = infile.readlines()
    
# with open("Vi33.19_max150.txt", "w") as outfile:
#     for line in lines:
#         num = line.rstrip()
#         if int(num) > 150:
#             num = 150
#         print(num, file=outfile)

with open("Vi33.19_max150.txt", "r") as infile:
    lines = [int(line.rstrip()) for line in infile.readlines()]
    avg_len = sum(lines)/len(lines)
    print(avg_len)
        