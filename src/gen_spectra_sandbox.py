import gen_spectra
from utils import ppm_to_da
from sqlite import database_file

ppm_tolerance = 20
max_len = 10


spectra_id = 
sequence = "DLQTLALEVE"
rival_seq = "EIAEKALVE"
obs_prec = 565.805089
prec_charge = 2
theoretical_prec = gen_spectra.get_precursor(sequence, prec_charge)
rival_theoretical_prec = gen_spectra.get_precursor(sequence, prec_charge)
print(theoretical_prec, rival_theoretical_prec, obs_prec)

# b_mass = gen_spectra.max_mass("DLQTLA", 'b', 1)
# y_mass = gen_spectra.max_mass("WSRM", 'y', 1)

# this_prec = gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, 1, 1, 2)
# print(this_prec)

# conv_b, conv_y = gen_spectra.convert_precursor_to_ion(obs_prec, prec_charge)

# b_mass = gen_spectra.max_mass(sequence, 'b', 2)
# y_mass = gen_spectra.max_mass(sequence, 'y', 2)


# btol = ppm_to_da(conv_b, ppm_tolerance)
# ytol = ppm_to_da(conv_y, ppm_tolerance)

# print(conv_b-btol, b_mass, conv_b+btol)
# print(conv_y-ytol, y_mass, conv_y+ytol)

# #Now to test the queries
# dbf = database_file(max_len, False)
# bqueries,_ = dbf.query_mass(conv_b, btol)
# _,yqueries = dbf.query_mass(y_mass, ytol)
# for i, query in enumerate(bqueries):
#     mass = query[0]
#     if mass == b_mass:
#         print("b_precursor found and entry at", i, query)

# #535.771854285 is a match with y_mass
# for i, query in enumerate(yqueries):
#     mass = query[0]
#     if mass == y_mass:
#         print("y_precursor found and entry at", i, query)