
timing_data_path = "/home/naco3124/snakemake/hypedsearch/src/Timing_data.txt"

with open(timing_data_path, 'r') as t:
    separation_list = []
    for line in t:
        separation = line.split(":\t")
        if float(separation[1].replace("\n", "")) > 40:
            separation_list.append(separation)
            
with open(timing_data_path, 'r') as t: #Looking for the smallest holdup
    covered_spectra = []
    for line in t:
        if "identification" not in line:
            continue
        separation = line.split(":")
        spec_number = separation[1].split(" ")[1]
        covered_spectra.append(int(spec_number))

with open(timing_data_path, 'r') as t: #Looking for the #complete
    count = 0
    for line in t:
        if "identification" not in line:
            continue
        count += 1


separation_list = sorted(separation_list, key = lambda x: x[1], reverse = True)

smallest_x = 1000000
for x in range(0,5000):
    if x not in covered_spectra:
        if x < smallest_x:
            smallest_x = x

print(separation_list[0], smallest_x, count)